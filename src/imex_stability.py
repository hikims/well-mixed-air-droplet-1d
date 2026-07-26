"""Linear stability of the exact first-order IMEX update used by model_lu_tbb.hpp.

This module is a diagnostic companion to the continuous-time H(lambda)=1
calculation.  It linearizes the actual code path

    explicit reaction -> implicit periodic diffusion -> exact Wbar update
    -> explicit droplet-position update

on the odd field sector.  Odd perturbations do not change Wbar at linear
order, so the reduced tangent map contains the three odd chemical fields and
the droplet displacement only.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from hopf_solver import DropletHopfModel


@dataclass
class IMEXMode:
    multiplier: complex
    rate: complex
    other_rates: np.ndarray
    info: dict[str, Any]


def _indicator_and_center_derivative(model: DropletHopfModel) -> tuple[np.ndarray, np.ndarray]:
    p = model.p
    xi = model.xi
    eps = p.eps_chi
    left = 0.5 * (np.tanh((xi + p.r) / eps) + 1.0)
    right = 0.5 * (np.tanh((p.r - xi) / eps) + 1.0)
    dleft_dxi = 0.5 / eps / np.cosh((xi + p.r) / eps) ** 2
    dright_dxi = -0.5 / eps / np.cosh((p.r - xi) / eps) ** 2
    chi = left * right
    dchi_dxi = dleft_dxi * right + left * dright_dxi
    return chi, -dchi_dxi  # d chi / d x_c


def _reaction_odd_data(model: DropletHopfModel, state: np.ndarray):
    n = model.nx
    p = model.p
    u = state[:n]
    s = state[n : 2 * n]
    c = state[2 * n : 3 * n]
    wbar = float(state[-1])
    _, q0, qs, qc = model.reaction_data(u, s, c)

    a11 = (
        -model.chi_d * p.k0 / p.u0
        - p.ka * wbar / p.u0
        - model.chi_u * p.ku * (1.0 - wbar / p.w0)
        - p.kus * q0
    )
    a12 = -p.kus * qs
    a13 = p.kus * qc
    a21 = -p.lambda_us * q0
    a22 = -p.lambda_us * qs - p.lambda1
    a23 = p.lambda_us * qc
    a31 = p.lambda_us * q0
    a32 = p.lambda_us * qs
    a33 = -p.lambda_us * qc - p.lambda2

    reaction_full = sp.bmat(
        [
            [sp.diags(a11), sp.diags(a12), sp.diags(a13)],
            [sp.diags(a21), sp.diags(a22), sp.diags(a23)],
            [sp.diags(a31), sp.diags(a32), sp.diags(a33)],
        ],
        format="csc",
    )
    b3 = model._odd_basis_3
    reaction_odd = (b3.T @ reaction_full @ b3).tocsc()

    _, dchi_dxc = _indicator_and_center_derivative(model)
    forcing_u = dchi_dxc * (
        p.k0 * (1.0 - u / p.u0)
        + p.ku * (1.0 - wbar / p.w0) * u
    )
    forcing_full = np.concatenate((forcing_u, np.zeros(n), np.zeros(n)))
    forcing_odd = np.asarray(b3.T @ forcing_full).ravel()

    # Velocity derivative with respect to odd field perturbations.  For an odd
    # perturbation, the left interpolation is the negative of the right one.
    j0, j1, weight0, weight1 = model._edge_interpolation(model.center + p.r)
    u_edge = weight0 * u[j0] + weight1 * u[j1]
    s_edge = weight0 * s[j0] + weight1 * s[j1]
    gamma_u, gamma_s = model.gamma_derivatives(float(u_edge), float(s_edge))
    velocity_row_full = np.zeros(3 * n)
    prefactor = 2.0 / p.eta
    velocity_row_full[j0] += prefactor * gamma_u * weight0
    velocity_row_full[j1] += prefactor * gamma_u * weight1
    velocity_row_full[n + j0] += prefactor * gamma_s * weight0
    velocity_row_full[n + j1] += prefactor * gamma_s * weight1
    velocity_row_odd = np.asarray(velocity_row_full @ b3).ravel()

    # Explicit center derivative of the piecewise-linear edge sampling while
    # the grid fields are held fixed, exactly as in interp_periodic_field().
    def edge(target: float):
        coordinate = target % p.L
        scaled = coordinate / model.dx
        left_index = int(np.floor(scaled)) % n
        fraction = float(scaled - np.floor(scaled))
        right_index = (left_index + 1) % n
        return left_index, right_index, 1.0 - fraction, fraction

    jr0, jr1, ar, br = edge(model.center + p.r)
    jl0, jl1, al, bl = edge(model.center - p.r)
    ur = ar * u[jr0] + br * u[jr1]
    sr = ar * s[jr0] + br * s[jr1]
    ul = al * u[jl0] + bl * u[jl1]
    sl = al * s[jl0] + bl * s[jl1]
    gur, gsr = model.gamma_derivatives(float(ur), float(sr))
    gul, gsl = model.gamma_derivatives(float(ul), float(sl))
    ux_r = (u[jr1] - u[jr0]) / model.dx
    sx_r = (s[jr1] - s[jr0]) / model.dx
    ux_l = (u[jl1] - u[jl0]) / model.dx
    sx_l = (s[jl1] - s[jl0]) / model.dx
    velocity_center_derivative = (
        gur * ux_r + gsr * sx_r - gul * ux_l - gsl * sx_l
    ) / p.eta

    return reaction_odd, forcing_odd, velocity_row_odd, float(velocity_center_derivative)


def build_continuous_lab_odd_operator(model: DropletHopfModel, state: np.ndarray) -> sp.csc_matrix:
    """Continuous-time semi-discrete odd operator including droplet displacement."""
    reaction_odd, forcing_odd, velocity_row_odd, velocity_x = _reaction_odd_data(
        model, state
    )
    odd_laplacian = (model._odd_basis.T @ model.laplacian @ model._odd_basis).tocsc()
    diffusion_odd = sp.block_diag(
        (
            model.p.du * odd_laplacian,
            model.p.ds * odd_laplacian,
            model.p.dc * odd_laplacian,
        ),
        format="csc",
    )
    field_operator = diffusion_odd + reaction_odd
    return sp.bmat(
        [
            [field_operator, sp.csc_matrix(forcing_odd[:, None])],
            [sp.csc_matrix(velocity_row_odd[None, :]), sp.csc_matrix([[velocity_x]])],
        ],
        format="csc",
    )


def leading_continuous_mode(model: DropletHopfModel, state: np.ndarray, number: int = 12) -> IMEXMode:
    operator = build_continuous_lab_odd_operator(model, state)
    values = spla.eigs(
        operator,
        k=min(number, operator.shape[0] - 2),
        which="LR",
        return_eigenvectors=False,
        tol=1.0e-10,
        maxiter=60000,
        ncv=min(operator.shape[0], max(4 * number + 8, 48)),
    )
    candidates = [value for value in values if value.imag > 10.0]
    leading = max(candidates, key=lambda value: value.real)
    ordered = np.asarray(sorted(values, key=lambda value: value.real, reverse=True))
    return IMEXMode(
        multiplier=complex(np.nan, np.nan),
        rate=complex(leading),
        other_rates=ordered,
        info={"dimension": int(operator.shape[0])},
    )


def build_imex_odd_map(
    model: DropletHopfModel, state: np.ndarray, dt: float
) -> tuple[spla.LinearOperator, dict[str, Any]]:
    """LinearOperator for one exact code time step on odd fields plus delta x."""
    if dt <= 0.0:
        raise ValueError("dt must be positive")
    reaction_odd, forcing_odd, velocity_row_odd, velocity_x = _reaction_odd_data(
        model, state
    )
    odd_laplacian = (model._odd_basis.T @ model.laplacian @ model._odd_basis).tocsc()
    m = odd_laplacian.shape[0]
    factorizations = [
        spla.splu(
            (
                sp.eye(m, format="csc")
                - dt * diffusion * odd_laplacian
            ).tocsc()
        )
        for diffusion in (model.p.du, model.p.ds, model.p.dc)
    ]
    dimension = 3 * m + 1

    def solve_real_factor(factor: spla.SuperLU, rhs: np.ndarray) -> np.ndarray:
        if np.iscomplexobj(rhs):
            return factor.solve(rhs.real) + 1j * factor.solve(rhs.imag)
        return factor.solve(rhs)

    def matvec(vector: np.ndarray) -> np.ndarray:
        vector = np.asarray(vector)
        fields = vector[:-1]
        displacement = vector[-1]
        rhs = fields + dt * (reaction_odd @ fields + forcing_odd * displacement)
        updated = np.empty_like(fields, dtype=np.result_type(vector, complex))
        for species, factor in enumerate(factorizations):
            sl = slice(species * m, (species + 1) * m)
            updated[sl] = solve_real_factor(factor, rhs[sl])
        new_displacement = displacement + dt * (
            velocity_row_odd @ updated + velocity_x * displacement
        )
        return np.concatenate((updated, [new_displacement]))

    operator = spla.LinearOperator(
        (dimension, dimension), matvec=matvec, dtype=np.complex128
    )
    info = {
        "dimension": dimension,
        "odd_field_dimension": 3 * m,
        "velocity_center_derivative": velocity_x,
        "dt": dt,
    }
    return operator, info


def leading_imex_mode(
    model: DropletHopfModel,
    state: np.ndarray,
    dt: float,
    number: int = 12,
) -> IMEXMode:
    operator, info = build_imex_odd_map(model, state, dt)
    values = spla.eigs(
        operator,
        k=min(number, operator.shape[0] - 2),
        which="LM",
        return_eigenvectors=False,
        tol=1.0e-10,
        maxiter=60000,
        ncv=min(operator.shape[0], max(4 * number + 8, 48)),
    )
    rates = np.log(values) / dt
    candidates = [
        (multiplier, rate)
        for multiplier, rate in zip(values, rates)
        if rate.imag > 10.0
    ]
    multiplier, rate = max(candidates, key=lambda pair: pair[1].real)
    order = np.argsort(rates.real)[::-1]
    return IMEXMode(
        multiplier=complex(multiplier),
        rate=complex(rate),
        other_rates=np.asarray(rates[order]),
        info=info,
    )
