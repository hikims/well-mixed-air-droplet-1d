#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.optimize import brentq

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(ROOT / "src"))

from hopf_solver import DropletHopfModel, Parameters  # noqa: E402


def calculate_at_s0(model: DropletHopfModel, s0: float, number: int = 12):
    state, steady_info = model.solve_steady_state(s0)
    open_loop, closed_loop, q, edge_row, edge_info = model.build_odd_operators(state)
    values, vectors = model.rightmost_eigenpairs(
        closed_loop, number=number, eigenvectors=True
    )
    leading = model.leading_oscillatory_root(values)
    return {
        "state": state,
        "steady_info": steady_info,
        "open_loop": open_loop,
        "closed_loop": closed_loop,
        "q": q,
        "edge_row": edge_row,
        "edge_info": edge_info,
        "values": values,
        "vectors": vectors,
        "leading": leading,
    }


def write_scan(model: DropletHopfModel, s_values: np.ndarray, output: Path) -> pd.DataFrame:
    rows = []
    for index, s0 in enumerate(s_values, start=1):
        result = calculate_at_s0(model, float(s0), number=10)
        leading = result["leading"]
        values = result["values"]
        real_nonpair = [v.real for v in values if abs(v.imag) < 1.0]
        next_real = max(real_nonpair) if real_nonpair else np.nan
        rows.append(
            {
                "nx": model.nx,
                "s0": float(s0),
                "leading_real": leading.real,
                "leading_imag": leading.imag,
                "leading_frequency_abs": abs(leading.imag),
                "next_real_eigenvalue": next_real,
                "steady_residual_inf": result["steady_info"]["residual_inf"],
                "wbar": result["steady_info"]["wbar"],
            }
        )
        print(
            f"nx={model.nx}: scan {index:02d}/{len(s_values):02d}, "
            f"s0={s0:.6f}, lambda={leading.real:+.8f}{leading.imag:+.8f}i"
        )
    frame = pd.DataFrame(rows)
    frame.to_csv(output, index=False)
    return frame


def critical_point(model: DropletHopfModel, bracket: tuple[float, float]):
    def real_part(s0: float) -> float:
        result = calculate_at_s0(model, s0, number=8)
        value = result["leading"].real
        print(f"nx={model.nx}: root evaluation s0={s0:.12f}, Re(lambda)={value:+.6e}")
        return value

    s_h = brentq(real_part, bracket[0], bracket[1], xtol=1.0e-11, rtol=1.0e-12)
    result = calculate_at_s0(model, s_h, number=24)
    lam = result["leading"]

    difference = 1.0e-4
    below = calculate_at_s0(model, s_h - difference, number=8)["leading"]
    above = calculate_at_s0(model, s_h + difference, number=8)["leading"]
    slope = (above.real - below.real) / (2.0 * difference)
    h_value = model.response_value(
        lam, result["open_loop"], result["q"], result["edge_row"]
    )
    values = result["values"]
    pair_index = int(np.argmin(np.abs(values - lam)))
    vector_odd = result["vectors"][:, pair_index]
    vector_full = model.reconstruct_odd_vector(vector_odd)
    phase_index = int(np.argmax(np.abs(vector_full)))
    vector_full *= np.exp(-1j * np.angle(vector_full[phase_index]))
    vector_full /= np.max(np.abs(vector_full))

    real_other = [v.real for v in values if abs(v - lam) > 1.0e-5 and abs(v - np.conj(lam)) > 1.0e-5]
    next_real = max(real_other)
    summary = {
        "nx": model.nx,
        "dx": model.dx,
        "s_h": float(s_h),
        "lambda_real": float(lam.real),
        "omega_h": float(lam.imag),
        "period_h": float(2.0 * np.pi / abs(lam.imag)),
        "crossing_slope": float(slope),
        "H_real": float(h_value.real),
        "H_imag": float(h_value.imag),
        "H_minus_1_abs": float(abs(h_value - 1.0)),
        "next_rightmost_real_part": float(next_real),
        "steady_residual_inf": float(result["steady_info"]["residual_inf"]),
        "wbar": float(result["steady_info"]["wbar"]),
        **result["edge_info"],
    }
    return summary, result, vector_full


def save_critical_outputs(
    model: DropletHopfModel,
    summary: dict,
    result: dict,
    vector_full: np.ndarray,
    output_dir: Path,
    response_points: int = 161,
):
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "critical_point.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    values = result["values"]
    pd.DataFrame(
        {
            "rank_by_real_part": np.arange(1, len(values) + 1),
            "real": values.real,
            "imag": values.imag,
        }
    ).to_csv(output_dir / "rightmost_spectrum.csv", index=False)

    n = model.nx
    state = result["state"]
    profile = pd.DataFrame(
        {
            "x": model.x,
            "xi": model.xi,
            "chi_d": model.chi_d,
            "u": state[:n],
            "s": state[n : 2 * n],
            "c": state[2 * n : 3 * n],
            "wbar": state[-1],
        }
    ).sort_values("xi")
    profile.to_csv(output_dir / "steady_profile_critical.csv", index=False)

    eigenfunction = pd.DataFrame(
        {
            "x": model.x,
            "xi": model.xi,
            "phi_u_real": vector_full[:n].real,
            "phi_u_imag": vector_full[:n].imag,
            "phi_s_real": vector_full[n : 2 * n].real,
            "phi_s_imag": vector_full[n : 2 * n].imag,
            "phi_c_real": vector_full[2 * n : 3 * n].real,
            "phi_c_imag": vector_full[2 * n : 3 * n].imag,
        }
    ).sort_values("xi")
    eigenfunction.to_csv(output_dir / "critical_eigenfunction.csv", index=False)

    omega_h = summary["omega_h"]
    if response_points <= 1:
        omega_values = np.array([omega_h])
    else:
        omega_values = np.linspace(
            omega_h - 15.0, omega_h + 15.0, response_points
        )
    response_values = np.array(
        [
            model.response_value(
                1j * omega, result["open_loop"], result["q"], result["edge_row"]
            )
            for omega in omega_values
        ]
    )
    pd.DataFrame(
        {
            "omega": omega_values,
            "H_real": response_values.real,
            "H_imag": response_values.imag,
            "H_real_minus_1": response_values.real - 1.0,
            "H_minus_1_abs": np.abs(response_values - 1.0),
        }
    ).to_csv(output_dir / "response_curve.csv", index=False)


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "results")
    parser.add_argument("--grids", type=int, nargs="+", default=[401, 801])
    parser.add_argument("--s-min", type=float, default=0.070)
    parser.add_argument("--s-max", type=float, default=0.105)
    parser.add_argument("--s-step", type=float, default=0.001)
    parser.add_argument("--bracket", type=float, nargs=2, default=[0.077, 0.081])
    parser.add_argument(
        "--force", action="store_true",
        help="recompute outputs even when complete files already exist",
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    args.output.mkdir(parents=True, exist_ok=True)
    parameters = Parameters()
    s_values = np.arange(args.s_min, args.s_max + 0.5 * args.s_step, args.s_step)
    summaries: list[dict] = []

    for nx in args.grids:
        model = DropletHopfModel(nx=nx, parameters=parameters)
        grid_dir = args.output / f"Nx{nx}"
        grid_dir.mkdir(parents=True, exist_ok=True)

        scan_file = grid_dir / "eigenvalue_scan.csv"
        if scan_file.exists() and not args.force:
            scan = pd.read_csv(scan_file)
            expected = np.asarray(s_values)
            complete = (
                len(scan) == len(expected)
                and np.allclose(scan["s0"].to_numpy(), expected)
            )
            if complete:
                print(f"nx={nx}: reusing complete scan {scan_file}")
            else:
                write_scan(model, s_values, scan_file)
        else:
            write_scan(model, s_values, scan_file)

        critical_file = grid_dir / "critical_point.json"
        required = [
            critical_file,
            grid_dir / "rightmost_spectrum.csv",
            grid_dir / "steady_profile_critical.csv",
            grid_dir / "critical_eigenfunction.csv",
            grid_dir / "response_curve.csv",
        ]
        if all(path.exists() for path in required) and not args.force:
            summary = json.loads(critical_file.read_text(encoding="utf-8"))
            print(f"nx={nx}: reusing complete critical-point outputs")
        else:
            summary, result, vector_full = critical_point(
                model, (args.bracket[0], args.bracket[1])
            )
            # The full response sweep is plotted only for the baseline grid.
            # On the refined grid, save the direct H(i omega_H) check only;
            # this avoids many repeated sparse complex factorizations.
            response_points = 161 if nx == min(args.grids) else 1
            save_critical_outputs(
                model, summary, result, vector_full, grid_dir, response_points
            )
        summaries.append(summary)

    summary_frame = pd.DataFrame(summaries)
    summary_frame.to_csv(args.output / "critical_points.csv", index=False)
    if len(summary_frame) >= 2:
        coarse = summary_frame.iloc[-2]
        fine = summary_frame.iloc[-1]
        extrapolated_s = fine.s_h + (fine.s_h - coarse.s_h) / 3.0
        extrapolated_omega = fine.omega_h + (fine.omega_h - coarse.omega_h) / 3.0
    else:
        extrapolated_s = np.nan
        extrapolated_omega = np.nan

    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "grids": args.grids,
        "scan": {"s_min": args.s_min, "s_max": args.s_max, "s_step": args.s_step},
        "root_bracket": args.bracket,
        "indicator": "0.25*(1+tanh((xi+r)/eps_chi))*(1+tanh((r-xi)/eps_chi))",
        "droplet_center": "L/2, matching initialize_state(..., centered=true)",
        "edge_interpolation": "periodic piecewise-linear interpolation at x_c+r",
        "parameters": parameters.to_dict(),
        "second_order_richardson_s_h": float(extrapolated_s),
        "second_order_richardson_omega_h": float(extrapolated_omega),
    }
    with (args.output / "run_metadata.json").open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)

    print("\nCritical-point summary")
    print(summary_frame.to_string(index=False))


if __name__ == "__main__":
    main()
