#!/usr/bin/env python3
"""Compare the continuous Hopf point with the actual first-order IMEX map."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(ROOT / "src"))

from hopf_solver import DropletHopfModel, Parameters  # noqa: E402
from imex_stability import leading_continuous_mode, leading_imex_mode  # noqa: E402


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "time_step_stability.csv")
    parser.add_argument("--nx", type=int, default=401)
    parser.add_argument(
        "--dts",
        type=float,
        nargs="+",
        default=[5.0e-4, 2.5e-4, 1.25e-4, 6.25e-5],
    )
    parser.add_argument("--scan-min", type=float, default=0.075)
    parser.add_argument("--scan-max", type=float, default=0.110)
    parser.add_argument("--scan-step", type=float, default=0.0025)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def locate_bracket(function, values):
    previous_x = float(values[0])
    previous_y = float(function(previous_x))
    for current_x in values[1:]:
        current_x = float(current_x)
        current_y = float(function(current_x))
        if previous_y == 0.0 or previous_y * current_y < 0.0:
            return previous_x, current_x
        previous_x, previous_y = current_x, current_y
    raise RuntimeError("No sign change found in the requested scan interval")


def main():
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.output.exists() and not args.force:
        print(f"reusing {args.output}")
        print(pd.read_csv(args.output).to_string(index=False))
        return
    parameters = Parameters()
    # A grid-aligned center is the stable pinning position selected by the
    # production lab-frame code at finite Nx.
    model = DropletHopfModel(nx=args.nx, parameters=parameters, center=0.0)
    cache = {}

    def state_at(s0):
        key = round(float(s0), 14)
        if key not in cache:
            cache[key] = model.solve_steady_state(float(s0))[0]
        return cache[key]

    scan = np.arange(args.scan_min, args.scan_max + 0.5 * args.scan_step, args.scan_step)
    rows = []

    def continuous_growth(s0):
        mode = leading_continuous_mode(model, state_at(s0))
        print(f"continuous s0={s0:.9f}, growth={mode.rate.real:+.8e}")
        return mode.rate.real

    bracket = locate_bracket(continuous_growth, scan)
    critical = brentq(continuous_growth, *bracket, xtol=1.0e-10, rtol=1.0e-11)
    mode = leading_continuous_mode(model, state_at(critical), number=20)
    rows.append(
        {
            "scheme": "continuous semi-discrete lab-frame",
            "nx": args.nx,
            "dt": 0.0,
            "critical_s0": critical,
            "growth_rate": mode.rate.real,
            "omega": mode.rate.imag,
            "period": 2.0 * np.pi / abs(mode.rate.imag),
            "multiplier_abs": np.nan,
        }
    )

    for dt in args.dts:
        def growth(s0):
            mode = leading_imex_mode(model, state_at(s0), dt)
            print(f"dt={dt:.8g}, s0={s0:.9f}, growth={mode.rate.real:+.8e}")
            return mode.rate.real

        bracket = locate_bracket(growth, scan)
        critical = brentq(growth, *bracket, xtol=1.0e-10, rtol=1.0e-11)
        mode = leading_imex_mode(model, state_at(critical), dt, number=20)
        rows.append(
            {
                "scheme": "production IMEX tangent map",
                "nx": args.nx,
                "dt": dt,
                "critical_s0": critical,
                "growth_rate": mode.rate.real,
                "omega": mode.rate.imag,
                "period": 2.0 * np.pi / abs(mode.rate.imag),
                "multiplier_abs": abs(mode.multiplier),
            }
        )

    frame = pd.DataFrame(rows)
    frame.to_csv(args.output, index=False)
    print("\n", frame.to_string(index=False))


if __name__ == "__main__":
    main()
