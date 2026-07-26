"""Odd-sector eigenvalue calculation for the well-mixed-air droplet model.

The implementation follows Eqs. (3)--(6) and the uploaded C++ implementation.
It reproduces the periodic second-order finite differences, the exact product
of tanh switches used by sigmoid(), and the periodic linear edge interpolation.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


@dataclass(frozen=True)
class Parameters:
    L: float = 16.0
    r: float = 1.0
    eta: float = 2.0
    du: float = 48.0
    ds: float = 9.6
    dc: float = 1.6
    k0: float = 0.9
    ku: float = 80.0
    ka: float = 1.0
    kw: float = 10000.0
    kus: float = 1.0
    lambda_us: float = 1.0
    lambda1: float = 0.8
    lambda2: float = 5.0
    u0: float = 0.006
    w0: float = 1.0
    c0: float = 0.02
    alpha: float = 1.0
    beta: float = 0.002
    a: float = 0.06
    b: float = 0.06
    gamma0: float = 80000.0
    gamma1: float = 0.0
    eps_chi: float = 0.04

    def to_dict(self) -> dict[str, float]:
        return asdict(self)


def periodic_laplacian(n: int, dx: float) -> sp.csc_matrix:
    main = -2.0 * np.ones(n)
    off = np.ones(n - 1)
    matrix = sp.diags([off, main, off], [-1, 0, 1], shape=(n, n), format="lil")
    matrix[0, n - 1] = 1.0
    matrix[n - 1, 0] = 1.0
    return (matrix.tocsc()) / dx**2


def periodic_derivative(n: int, dx: float) -> sp.csc_matrix:
    off = np.ones(n - 1)
    matrix = sp.diags([-off, off], [-1, 1], shape=(n, n), format="lil")
    matrix[0, n - 1] = -1.0
    matrix[n - 1, 0] = 1.0
    return matrix.tocsc() / (2.0 * dx)


class DropletHopfModel:
    def __init__(
        self,
        nx: int = 401,
        parameters: Parameters | None = None,
        center: float | None = None,
    ):
        if nx < 5 or nx % 2 != 1:
            raise ValueError("nx must be an odd integer >= 5")
        self.nx = int(nx)
        self.p = parameters or Parameters()
        self.dx = self.p.L / self.nx
        self.x = np.arange(self.nx, dtype=float) * self.dx
        # The production code fixes a centered rest droplet at x_c=L/2.
        self.center = (0.5 * self.p.L) if center is None else float(center) % self.p.L
        self.xi = ((self.x - self.center + self.p.L / 2.0) % self.p.L) - self.p.L / 2.0
        # Exact one-droplet form of sigmoid(x, center, rad) in model_lu_tbb.hpp.
        left_switch = 0.5 * (
            np.tanh((self.xi + self.p.r) / self.p.eps_chi) + 1.0
        )
        right_switch = 0.5 * (
            np.tanh((self.p.r - self.xi) / self.p.eps_chi) + 1.0
        )
        self.chi_d = left_switch * right_switch
        self.chi_u = 1.0 - self.chi_d
        self.laplacian = periodic_laplacian(self.nx, self.dx)
        self.derivative = periodic_derivative(self.nx, self.dx)
        self._odd_basis = self._build_odd_basis()
        self._odd_basis_3 = sp.block_diag(
            (self._odd_basis, self._odd_basis, self._odd_basis), format="csc"
        )

    def _build_odd_basis(self) -> sp.csc_matrix:
        """Orthonormal basis of f(-xi)=-f(xi) on an odd periodic grid."""
        m = (self.nx - 1) // 2
        rows: list[int] = []
        cols: list[int] = []
        data: list[float] = []
        coefficient = 1.0 / np.sqrt(2.0)
        for k, j in enumerate(range(1, m + 1)):
            rows.extend((j, self.nx - j))
            cols.extend((k, k))
            data.extend((coefficient, -coefficient))
        return sp.csc_matrix((data, (rows, cols)), shape=(self.nx, m))

    def reaction_data(
        self, u: np.ndarray, s: np.ndarray, c: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        p = self.p
        f = p.alpha * s / (p.beta + s)
        fp = p.alpha * p.beta / (p.beta + s) ** 2
        unsaturated = 1.0 - c / p.c0
        reaction = unsaturated * f * u
        q0 = unsaturated * f
        qs = unsaturated * fp * u
        qc = f * u / p.c0
        return reaction, q0, qs, qc

    def residual_and_jacobian(
        self, state: np.ndarray, s0: float
    ) -> tuple[np.ndarray, sp.csc_matrix]:
        n = self.nx
        p = self.p
        u = state[:n]
        s = state[n : 2 * n]
        c = state[2 * n : 3 * n]
        wbar = float(state[-1])
        reaction, q0, qs, qc = self.reaction_data(u, s, c)

        f_u = (
            p.du * (self.laplacian @ u)
            + self.chi_d * p.k0 * (1.0 - u / p.u0)
            + p.ka * wbar * (1.0 - u / p.u0)
            - self.chi_u * p.ku * (1.0 - wbar / p.w0) * u
            - p.kus * reaction
        )
        f_s = (
            p.ds * (self.laplacian @ s)
            - p.lambda_us * reaction
            + p.lambda1 * (s0 - s)
        )
        f_c = (
            p.dc * (self.laplacian @ c)
            + p.lambda_us * reaction
            - p.lambda2 * c
        )
        phi = np.mean(self.chi_u * u)
        psi = np.mean(1.0 - u / p.u0)
        f_w = (
            -p.kw * wbar
            + p.ku * (1.0 - wbar / p.w0) * phi
            - p.ka * wbar * psi
        )
        residual = np.concatenate((f_u, f_s, f_c, [f_w]))

        a11 = (
            -self.chi_d * p.k0 / p.u0
            - p.ka * wbar / p.u0
            - self.chi_u * p.ku * (1.0 - wbar / p.w0)
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

        j_uu = p.du * self.laplacian + sp.diags(a11, format="csc")
        j_us = sp.diags(a12, format="csc")
        j_uc = sp.diags(a13, format="csc")
        j_su = sp.diags(a21, format="csc")
        j_ss = p.ds * self.laplacian + sp.diags(a22, format="csc")
        j_sc = sp.diags(a23, format="csc")
        j_cu = sp.diags(a31, format="csc")
        j_cs = sp.diags(a32, format="csc")
        j_cc = p.dc * self.laplacian + sp.diags(a33, format="csc")

        d_fu_dw = p.ka * (1.0 - u / p.u0) + self.chi_u * p.ku * u / p.w0
        zero_column = sp.csc_matrix((n, 1))
        u_column = sp.csc_matrix(d_fu_dw[:, None])
        d_fw_du = (
            p.ku * (1.0 - wbar / p.w0) * self.chi_u / n
            + p.ka * wbar / (n * p.u0)
        )
        d_fw_dw = -p.kw - (p.ku / p.w0) * phi - p.ka * psi
        u_row = sp.csc_matrix(d_fw_du[None, :])
        zero_row = sp.csc_matrix((1, n))
        scalar = sp.csc_matrix([[d_fw_dw]])

        jacobian = sp.bmat(
            [
                [j_uu, j_us, j_uc, u_column],
                [j_su, j_ss, j_sc, zero_column],
                [j_cu, j_cs, j_cc, zero_column],
                [u_row, zero_row, zero_row, scalar],
            ],
            format="csc",
        )
        return residual, jacobian

    def solve_steady_state(
        self,
        s0: float,
        initial_state: np.ndarray | None = None,
        tolerance: float = 1.0e-10,
        max_iterations: int = 30,
    ) -> tuple[np.ndarray, dict[str, Any]]:
        n = self.nx
        p = self.p
        if initial_state is None:
            state = np.concatenate(
                (np.zeros(n), np.full(n, s0), np.zeros(n), [0.0])
            )
        else:
            state = np.asarray(initial_state, dtype=float).copy()

        residual_norm = np.inf
        for iteration in range(max_iterations):
            residual, jacobian = self.residual_and_jacobian(state, s0)
            residual_norm = float(np.linalg.norm(residual, ord=np.inf))
            if residual_norm < tolerance:
                break
            correction = spla.spsolve(jacobian, -residual)
            step = 1.0
            accepted = False
            for _ in range(24):
                candidate = state + step * correction
                if (
                    np.min(candidate[n : 2 * n]) <= -0.99 * p.beta
                    or np.min(candidate[2 * n : 3 * n]) < -0.5 * p.c0
                    or candidate[-1] < -0.1 * p.w0
                ):
                    step *= 0.5
                    continue
                candidate_residual, _ = self.residual_and_jacobian(candidate, s0)
                candidate_norm = float(
                    np.linalg.norm(candidate_residual, ord=np.inf)
                )
                if np.isfinite(candidate_norm) and candidate_norm < residual_norm:
                    state = candidate
                    accepted = True
                    break
                step *= 0.5
            if not accepted:
                if residual_norm < 50.0 * tolerance:
                    break
                raise RuntimeError(
                    f"Newton line search failed for nx={n}, s0={s0:.12g}, "
                    f"residual={residual_norm:.3e}"
                )
        else:
            raise RuntimeError(
                f"Newton iteration did not converge for nx={n}, s0={s0:.12g}"
            )

        # Enforce the analytically exact reflection symmetry before differentiating.
        reflected = (-np.arange(n)) % n
        for offset in (0, n, 2 * n):
            values = state[offset : offset + n]
            state[offset : offset + n] = 0.5 * (values + values[reflected])
        residual, _ = self.residual_and_jacobian(state, s0)
        residual_norm = float(np.linalg.norm(residual, ord=np.inf))
        info = {
            "iterations": iteration,
            "residual_inf": residual_norm,
            "wbar": float(state[-1]),
        }
        return state, info

    def _edge_interpolation(self, target: float) -> tuple[int, int, float, float]:
        coordinate = target % self.p.L
        scaled = coordinate / self.dx
        j0 = int(np.floor(scaled)) % self.nx
        fraction = float(scaled - np.floor(scaled))
        j1 = (j0 + 1) % self.nx
        return j0, j1, 1.0 - fraction, fraction

    def gamma_derivatives(self, u: float, s: float) -> tuple[float, float]:
        p = self.p
        scale = p.gamma0 - p.gamma1
        gamma_u = -scale * p.a**2 * u / (p.a**2 + u**2) ** 2
        gamma_s = -0.5 * scale * p.b / (p.b + s) ** 2
        return float(gamma_u), float(gamma_s)

    def build_odd_operators(
        self, state: np.ndarray
    ) -> tuple[sp.csc_matrix, sp.csc_matrix, np.ndarray, np.ndarray, dict[str, float]]:
        n = self.nx
        p = self.p
        u = state[:n]
        s = state[n : 2 * n]
        c = state[2 * n : 3 * n]
        wbar = float(state[-1])
        _, q0, qs, qc = self.reaction_data(u, s, c)

        a11 = (
            -self.chi_d * p.k0 / p.u0
            - p.ka * wbar / p.u0
            - self.chi_u * p.ku * (1.0 - wbar / p.w0)
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

        open_loop_full = sp.bmat(
            [
                [
                    p.du * self.laplacian + sp.diags(a11),
                    sp.diags(a12),
                    sp.diags(a13),
                ],
                [
                    sp.diags(a21),
                    p.ds * self.laplacian + sp.diags(a22),
                    sp.diags(a23),
                ],
                [
                    sp.diags(a31),
                    sp.diags(a32),
                    p.dc * self.laplacian + sp.diags(a33),
                ],
            ],
            format="csc",
        )
        open_loop_odd = (
            self._odd_basis_3.T @ open_loop_full @ self._odd_basis_3
        ).tocsc()

        q_full = np.concatenate(
            (self.derivative @ u, self.derivative @ s, self.derivative @ c)
        )
        q_odd = np.asarray(self._odd_basis_3.T @ q_full).ravel()

        # Match interp_periodic_field(..., XC(0)+rad) in model_lu_tbb.hpp.
        j0, j1, weight0, weight1 = self._edge_interpolation(self.center + p.r)
        u_edge = weight0 * u[j0] + weight1 * u[j1]
        s_edge = weight0 * s[j0] + weight1 * s[j1]
        gamma_u, gamma_s = self.gamma_derivatives(float(u_edge), float(s_edge))
        edge_row_full = np.zeros(3 * n)
        prefactor = 2.0 / p.eta
        edge_row_full[j0] += prefactor * gamma_u * weight0
        edge_row_full[j1] += prefactor * gamma_u * weight1
        edge_row_full[n + j0] += prefactor * gamma_s * weight0
        edge_row_full[n + j1] += prefactor * gamma_s * weight1
        edge_row_odd = np.asarray(edge_row_full @ self._odd_basis_3).ravel()

        feedback = sp.csc_matrix(q_odd[:, None] @ edge_row_odd[None, :])
        closed_loop_odd = open_loop_odd + feedback
        edge_info = {
            "u_edge": float(u_edge),
            "s_edge": float(s_edge),
            "gamma_u_edge": gamma_u,
            "gamma_s_edge": gamma_s,
            "edge_j0": int(j0),
            "edge_j1": int(j1),
            "edge_weight0": float(weight0),
            "edge_weight1": float(weight1),
        }
        return open_loop_odd, closed_loop_odd, q_odd, edge_row_odd, edge_info

    def rightmost_eigenpairs(
        self,
        closed_loop_odd: sp.csc_matrix,
        number: int = 12,
        tolerance: float = 1.0e-10,
        eigenvectors: bool = False,
    ) -> tuple[np.ndarray, np.ndarray | None]:
        dimension = closed_loop_odd.shape[0]
        k = min(number, dimension - 2)
        values, vectors = spla.eigs(
            closed_loop_odd,
            k=k,
            which="LR",
            return_eigenvectors=True,
            tol=tolerance,
            maxiter=60000,
            ncv=min(dimension, max(4 * k + 8, 40)),
        )
        order = np.argsort(values.real)[::-1]
        values = values[order]
        vectors = vectors[:, order]
        return values, vectors if eigenvectors else None

    @staticmethod
    def leading_oscillatory_root(values: np.ndarray) -> complex:
        oscillatory = [value for value in values if value.imag > 1.0]
        if not oscillatory:
            return complex(values[np.argmax(values.real)])
        return complex(max(oscillatory, key=lambda value: value.real))

    @staticmethod
    def response_value(
        lam: complex,
        open_loop_odd: sp.csc_matrix,
        q_odd: np.ndarray,
        edge_row_odd: np.ndarray,
    ) -> complex:
        identity = sp.eye(open_loop_odd.shape[0], dtype=complex, format="csc")
        response = spla.spsolve(
            lam * identity - open_loop_odd.astype(complex), q_odd.astype(complex)
        )
        return complex(edge_row_odd @ response)

    def reconstruct_odd_vector(self, odd_vector: np.ndarray) -> np.ndarray:
        return np.asarray(self._odd_basis_3 @ odd_vector).ravel()
