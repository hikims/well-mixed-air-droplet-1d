
from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from matplotlib.ticker import ScalarFormatter

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
import numpy as np
import pandas as pd

# ============================================================
# Journal-style font settings matched to make_figures.py
# ============================================================

FONT_FAMILY = "Arial"
BASE_FONT_SIZE = 16
TICK_FONT_SIZE = 14
TITLE_FONT_SIZE = 18


def configure_figure_fonts() -> None:
    """Apply the same global font settings used by make_figures.py."""
    plt.rcParams.update({
        "font.family": FONT_FAMILY,
        "font.sans-serif": [FONT_FAMILY, "Helvetica", "DejaVu Sans"],
        "font.size": BASE_FONT_SIZE,
        "axes.labelsize": BASE_FONT_SIZE,
        "axes.titlesize": TITLE_FONT_SIZE,
        "xtick.labelsize": TICK_FONT_SIZE,
        "ytick.labelsize": TICK_FONT_SIZE,
        "legend.fontsize": BASE_FONT_SIZE,
        "legend.framealpha": 1.0,
        "figure.titlesize": TITLE_FONT_SIZE,
        "mathtext.fontset": "custom",
        "mathtext.rm": FONT_FAMILY,
        "mathtext.it": FONT_FAMILY + ":italic",
        "mathtext.bf": FONT_FAMILY + ":bold",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.unicode_minus": False,
    })


configure_figure_fonts()


def apply_figure_style(fig: plt.Figure) -> None:
    """Enforce the make_figures.py font sizes after artists are created."""
    for ax in fig.axes:
        ax.title.set_fontfamily(FONT_FAMILY)
        ax.title.set_fontsize(TITLE_FONT_SIZE)

        ax.xaxis.label.set_fontfamily(FONT_FAMILY)
        ax.xaxis.label.set_fontsize(BASE_FONT_SIZE)
        ax.yaxis.label.set_fontfamily(FONT_FAMILY)
        ax.yaxis.label.set_fontsize(BASE_FONT_SIZE)

        ax.tick_params(axis="both", which="major", labelsize=TICK_FONT_SIZE)
        ax.tick_params(axis="both", which="minor", labelsize=TICK_FONT_SIZE)

        for tick_label in ax.get_xticklabels() + ax.get_yticklabels():
            tick_label.set_fontfamily(FONT_FAMILY)
            tick_label.set_fontsize(TICK_FONT_SIZE)

        for text in ax.texts:
            text.set_fontfamily(FONT_FAMILY)
            text.set_fontsize(BASE_FONT_SIZE)

        legend = ax.get_legend()
        if legend is not None:
            for text in legend.get_texts():
                text.set_fontfamily(FONT_FAMILY)
                text.set_fontsize(BASE_FONT_SIZE)
            legend.get_title().set_fontfamily(FONT_FAMILY)
            legend.get_title().set_fontsize(BASE_FONT_SIZE)

    if getattr(fig, "_suptitle", None) is not None:
        fig._suptitle.set_fontfamily(FONT_FAMILY)
        fig._suptitle.set_fontsize(TITLE_FONT_SIZE)

    for text in fig.texts:
        if text is not getattr(fig, "_suptitle", None):
            text.set_fontfamily(FONT_FAMILY)
            text.set_fontsize(BASE_FONT_SIZE)


@dataclass(frozen=True)
class ScanSpec:
    label: str
    path: Path


def parse_scan(value: str) -> ScanSpec:
    if "=" not in value:
        raise argparse.ArgumentTypeError("--scan must be LABEL=PATH")
    label, raw_path = value.split("=", 1)
    label = label.strip()
    if not label:
        raise argparse.ArgumentTypeError("scan label cannot be empty")
    return ScanSpec(label=label, path=Path(raw_path).expanduser().resolve())


def estimate_period(t: np.ndarray, v: np.ndarray) -> float:
    if len(t) < 8:
        return math.nan

    centered = v - np.mean(v)
    if np.linalg.norm(centered) <= 1.0e-14:
        return math.nan

    crossings: list[float] = []
    for i in range(len(centered) - 1):
        y0 = float(centered[i])
        y1 = float(centered[i + 1])
        if y0 <= 0.0 < y1:
            dt = float(t[i + 1] - t[i])
            if dt <= 0.0:
                continue
            denom = y1 - y0
            frac = 0.0 if abs(denom) < 1.0e-30 else -y0 / denom
            crossings.append(float(t[i] + frac * dt))

    if len(crossings) >= 4:
        periods = np.diff(np.asarray(crossings, dtype=float))
        periods = periods[np.isfinite(periods) & (periods > 0.0)]
        if periods.size:
            return float(np.median(periods))

    # FFT fallback.
    sample_dt = np.diff(t)
    if sample_dt.size == 0 or np.any(sample_dt <= 0.0):
        return math.nan
    dt_med = float(np.median(sample_dt))
    window = np.hanning(len(centered))
    spectrum = np.abs(np.fft.rfft(centered * window))
    freqs = np.fft.rfftfreq(len(centered), d=dt_med)
    if len(spectrum) <= 1:
        return math.nan
    spectrum[0] = 0.0
    idx = int(np.argmax(spectrum))
    freq = float(freqs[idx])
    return 1.0 / freq if freq > 0.0 else math.nan


def read_hopf_reference(path: Path) -> dict[str, float]:
    frame = pd.read_csv(path)
    if frame.empty:
        raise ValueError(f"empty Hopf reference: {path}")
    row = frame.sort_values("nx").iloc[-1]
    return {
        "nx": int(row["nx"]),
        "s_h": float(row["s_h"]),
        "omega_h": float(row["omega_h"]),
        "period_h": float(row["period_h"]),
    }


def read_imex_reference(path: Path, target_dt: float) -> dict[str, float]:
    frame = pd.read_csv(path)
    frame = frame[frame["scheme"].astype(str).str.contains("IMEX", case=False, na=False)]
    if frame.empty:
        raise ValueError(f"no IMEX rows in {path}")
    idx = int(np.argmin(np.abs(frame["dt"].to_numpy(dtype=float) - target_dt)))
    row = frame.iloc[idx]
    if not math.isclose(float(row["dt"]), target_dt, rel_tol=1.0e-8, abs_tol=1.0e-12):
        raise ValueError(
            f"nearest IMEX dt={float(row['dt']):.12g} does not match target dt={target_dt:.12g}"
        )
    return {
        "nx": int(row["nx"]),
        "dt": float(row["dt"]),
        "critical_s0": float(row["critical_s0"]),
        "omega": float(row["omega"]),
        "period": float(row["period"]),
    }


def scientific_math_body(value: float) -> str:
    if value == 0.0:
        return "0"
    exponent = int(math.floor(math.log10(abs(value))))
    mantissa = value / (10 ** exponent)
    if math.isclose(mantissa, round(mantissa), rel_tol=1.0e-12, abs_tol=1.0e-12):
        mantissa_text = str(int(round(mantissa)))
    else:
        mantissa_text = f"{mantissa:.3g}"
    if mantissa_text == "1":
        return rf"10^{{{exponent}}}"
    return rf"{mantissa_text}\times10^{{{exponent}}}"


def epsilon_legend_label(eps: float) -> str:
    return rf"$\varepsilon={scientific_math_body(float(eps))}$"


def scan_epsilon_label(label: str, fallback_eps: float | None = None) -> str:
    lower = label.lower()
    if "1e-" in lower:
        idx = lower.find("1e-")
        token = []
        for ch in lower[idx:]:
            if ch.isdigit() or ch in {'e', '-', '+', '.'}:
                token.append(ch)
            else:
                break
        token = ''.join(token)
        try:
            return epsilon_legend_label(float(token))
        except ValueError:
            pass
    if fallback_eps is not None:
        return epsilon_legend_label(float(fallback_eps))
    return label


def analyze_scan(
    spec: ScanSpec,
    late_fraction: float,
    active_vrms: float,
    active_amp: float,
) -> pd.DataFrame:
    case_file = spec.path / "perturbation_cases.csv"
    if not case_file.exists():
        raise FileNotFoundError(case_file)
    cases = pd.read_csv(case_file)

    records: list[dict[str, object]] = []
    for row in cases.itertuples(index=False):
        trajectory = spec.path / str(row.folder) / "perturbation.csv"
        if not trajectory.exists():
            raise FileNotFoundError(trajectory)

        data = pd.read_csv(trajectory)
        data = data.replace([np.inf, -np.inf], np.nan)
        data = data.dropna(subset=["t", "delta_x", "velocity"]).sort_values("t")
        if len(data) < 8:
            raise ValueError(f"too few samples in {trajectory}")

        t0 = float(data["t"].iloc[0])
        t1 = float(data["t"].iloc[-1])
        late_start = t0 + (1.0 - late_fraction) * (t1 - t0)
        late = data[data["t"] >= late_start]
        if len(late) < 8:
            raise ValueError(f"too few late samples in {trajectory}")

        delta = late["delta_x"].to_numpy(dtype=float)
        velocity = late["velocity"].to_numpy(dtype=float)
        times = late["t"].to_numpy(dtype=float)

        late_amp = float(np.max(delta) - np.min(delta))
        late_vrms = float(np.sqrt(np.mean(velocity * velocity)))
        late_mean_abs_v = float(np.mean(np.abs(velocity)))
        max_abs_delta = float(np.max(np.abs(data["delta_x"].to_numpy(dtype=float))))
        max_abs_velocity = float(np.max(np.abs(data["velocity"].to_numpy(dtype=float))))
        active = bool(late_amp >= active_amp or late_vrms >= active_vrms)
        period = estimate_period(times, velocity) if active else math.nan

        records.append(
            {
                "scan": spec.label,
                "scan_dir": str(spec.path),
                "s0": float(row.s0),
                "folder": str(row.folder),
                "dt": float(row.dt),
                "eps_shift": float(row.eps_shift),
                "t_relax": float(row.t_relax),
                "t_free": float(row.t_free),
                "save_every": int(row.save_every),
                "center_mode": str(row.center_mode),
                "x_rest": float(row.x_rest),
                "nx": int(row.Nx),
                "dx": float(row.dx),
                "late_start": late_start,
                "late_fraction": late_fraction,
                "late_delta_amp": late_amp,
                "late_vrms": late_vrms,
                "late_mean_abs_velocity": late_mean_abs_v,
                "max_abs_delta": max_abs_delta,
                "max_abs_velocity": max_abs_velocity,
                "estimated_period": period,
                "active_by_threshold": active,
                "trajectory_file": str(trajectory),
            }
        )

    return pd.DataFrame.from_records(records).sort_values("s0").reset_index(drop=True)


def onset_bracket(frame: pd.DataFrame) -> tuple[float, float]:
    ordered = frame.sort_values("s0")
    active = ordered[ordered["active_by_threshold"]]
    if active.empty:
        return math.nan, math.nan
    first_active = float(active["s0"].iloc[0])
    inactive_below = ordered[
        (~ordered["active_by_threshold"]) & (ordered["s0"] < first_active)
    ]
    last_inactive = (
        float(inactive_below["s0"].iloc[-1]) if not inactive_below.empty else math.nan
    )
    return last_inactive, first_active



def setup_axes_grid(ax: plt.Axes) -> None:
    """Apply an EPS-safe grid without transparency."""
    ax.grid(True, color="0.85", linewidth=0.8)


def save_figure(fig: plt.Figure, output_stem: Path, formats: list[str]) -> None:
    """Save one figure in every requested format."""
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    apply_figure_style(fig)
    for fmt in formats:
        outfile = output_stem.with_suffix(f".{fmt}")
        fig.savefig(
            outfile,
            format=fmt,
            dpi=300 if fmt == "png" else None,
            bbox_inches="tight",
            pad_inches=0.10,
            facecolor="white",
            transparent=False,
        )


def plot_activity(
    combined: pd.DataFrame,
    y: str,
    ylabel: str,
    output: Path,
    hopf: dict[str, float],
    imex: dict[str, float],
    formats: list[str],
    show_hopf_reference: bool = True,
) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for label, group in combined.groupby("scan", sort=False):
        group = group.sort_values("s0")
        legend_label = scan_epsilon_label(label, float(group["eps_shift"].iloc[0]))
        ax.plot(group["s0"], group[y], marker="o", label=legend_label)
    if show_hopf_reference:
        ax.axvline(hopf["s_h"], linestyle=":", label=f"continuous Hopf, Nx={hopf['nx']}")
    ax.axvline(
        imex["critical_s0"],
        linestyle="--",
        label="IMEX threshold\n"
        + rf"$\Delta t={scientific_math_body(imex['dt'])}$",
    )
    ax.set_xlabel(r"$s_0$")
    ax.set_ylabel(ylabel)
    ax.set_yscale("log")

    # Paper-specific range and legend placement for the velocity-activity panel.
    if y == "late_vrms":
        # Restrict the right edge to the sampled range so the unused
        # 0.0830 tick is not displayed.
        ax.set_xlim(0.0815, 0.083)
        ax.xaxis.set_major_locator(MultipleLocator(0.0005))
        ax.xaxis.set_major_formatter(FormatStrFormatter("%.4f"))
        # The upper-left part of this panel is empty on the log scale.
        legend_loc = "upper left"
    else:
        legend_loc = "best"

    setup_axes_grid(ax)
    ax.legend(loc=legend_loc, frameon=True)
    fig.tight_layout()
    save_figure(fig, output, formats)
    plt.close(fig)


def plot_period(
    combined: pd.DataFrame,
    output: Path,
    hopf: dict[str, float],
    imex: dict[str, float],
    formats: list[str],
) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    plotted = False
    for label, group in combined.groupby("scan", sort=False):
        active = group[
            group["active_by_threshold"]
            & np.isfinite(group["estimated_period"].to_numpy(dtype=float))
        ].sort_values("s0")
        if active.empty:
            continue
        legend_label = scan_epsilon_label(label, float(group["eps_shift"].iloc[0]))
        ax.plot(active["s0"], active["estimated_period"], marker="o", label=legend_label)
        plotted = True
    ax.axhline(imex["period"], linestyle="--", label="IMEX linear period")
    ax.axvline(imex["critical_s0"], linestyle="--")
    ax.set_xlabel(r"$s_0$")
    ax.set_ylabel("nonlinear period")
    # Restrict the displayed range to the sampled near-onset data.
    # This removes the unused 0.0830 and 0.0832 tick labels.
    ax.set_xlim(0.08220, 0.08280)
    ax.set_xticks([0.0822, 0.0824, 0.0826, 0.0828])
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.4f"))
    ax.set_ylim(0.08160, 0.08200)
    
    setup_axes_grid(ax)
    if plotted:
        # The data occupy the upper part of this narrow window.
        ax.legend(loc="upper right", frameon=True)
    fig.tight_layout()
    save_figure(fig, output, formats)
    plt.close(fig)


def plot_brackets(
    combined: pd.DataFrame,
    output: Path,
    hopf: dict[str, float],
    imex: dict[str, float],
    formats: list[str],
) -> None:
    rows: list[dict[str, float | str]] = []
    for label, group in combined.groupby("scan", sort=False):
        last_inactive, first_active = onset_bracket(group)
        rows.append(
            {
                "scan": label,
                "eps": float(group["eps_shift"].iloc[0]),
                "last_inactive": last_inactive,
                "first_active": first_active,
            }
        )
    bracket = pd.DataFrame(rows).sort_values("eps")

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for row in bracket.itertuples(index=False):
        if math.isnan(row.first_active):
            continue
        if math.isnan(row.last_inactive):
            ax.scatter([row.eps], [row.first_active], label=epsilon_legend_label(float(row.eps)))
        else:
            mid = 0.5 * (row.last_inactive + row.first_active)
            half = 0.5 * (row.first_active - row.last_inactive)
            ax.errorbar(
                [row.eps],
                [mid],
                yerr=[[half], [half]],
                marker="o",
                capsize=5,
                linestyle="none",
                label=epsilon_legend_label(float(row.eps)),
            )
    ax.axhline(hopf["s_h"], linestyle=":", label="continuous Hopf")
    ax.axhline(imex["critical_s0"], linestyle="--", label="IMEX linear threshold")
    ax.set_xscale("log")
    ax.set_xlabel(r"position perturbation $\varepsilon$")
    ax.set_ylabel(r"sampled onset bracket in $s_0$")
    ax.set_ylim(0.078, 0.083)
    setup_axes_grid(ax)
    # The central log-x region is free of the two sampled epsilon points.
    ax.legend(loc="center", frameon=True, ncol=2)
    fig.tight_layout()
    save_figure(fig, output, formats)
    plt.close(fig)


def plot_transition_trajectories(
    frame: pd.DataFrame,
    output: Path,
    zoom_window: float,
    formats: list[str],
) -> None:
    last_inactive, first_active = onset_bracket(frame)
    selected: list[tuple[pd.Series, str]] = []
    if not math.isnan(last_inactive):
        selected.append((frame.iloc[int(np.argmin(np.abs(frame["s0"] - last_inactive)))], "decayed"))
    if not math.isnan(first_active):
        selected.append((frame.iloc[int(np.argmin(np.abs(frame["s0"] - first_active)))], "near-neutral"))
    if not selected:
        return

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    for row, state in selected:
        data = pd.read_csv(str(row["trajectory_file"]))
        t_end = float(data["t"].max())
        zoom = data[data["t"] >= t_end - zoom_window]
        ax.plot(
            zoom["t"] - float(zoom["t"].iloc[0]),
            zoom["delta_x"],
            label=rf"$s_0={float(row['s0']):.6f}$ ({state})",
        )
    ax.set_xlabel("time in final window")
    ax.set_ylabel(r"$\delta x(t)$")

    scan_label = str(frame["scan"].iloc[0])
    ax.set_title(scan_epsilon_label(scan_label, float(frame["eps_shift"].iloc[0])), pad=14)
    fmt = ScalarFormatter(useMathText=True)
    fmt.set_powerlimits((0, 0))
    ax.yaxis.set_major_formatter(fmt)
    label_lower = scan_label.lower()
    if "1e-3" in label_lower:
        ax.set_ylim(-1.5e-4, 2.5e-4)
        ax.set_yticks(np.linspace(-1.5e-4, 1.5e-4, 7))
    elif "1e-5" in label_lower:
        # This produces displayed tick labels from -2.0 to 2.0 with a 1e-6 scale factor.
        ax.set_ylim(-2.0e-6, 3.0e-6)
        ax.set_yticks(np.linspace(-2.0e-6, 2.0e-6, 5)) 

    setup_axes_grid(ax)
    # One-row legend in the unused top band of the prescribed y range.
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, 0.995),
        ncol=1,
        frameon=True,
        fontsize=BASE_FONT_SIZE,
    )
    fig.tight_layout()
    fig.subplots_adjust(top=0.90)
    save_figure(fig, output, formats)
    plt.close(fig)


def write_summary(
    combined: pd.DataFrame,
    output: Path,
    hopf: dict[str, float],
    imex: dict[str, float],
    active_vrms: float,
    active_amp: float,
) -> None:
    lines = [
        "# Grid-aligned near-onset perturbation validation",
        "",
        "## Linear references",
        "",
        f"- continuous odd-sector Hopf: `s_H = {hopf['s_h']:.12g}`, "
        f"`T_H = {hopf['period_h']:.12g}` (`Nx={hopf['nx']}`)",
        f"- production IMEX map at `dt = {imex['dt']:.12g}`: "
        f"`s_crit = {imex['critical_s0']:.12g}`, `T = {imex['period']:.12g}`",
        "",
        "## Diagnostic definition",
        "",
        f"A case is flagged active when `late_vrms >= {active_vrms:g}` or "
        f"`late_delta_amp >= {active_amp:g}`.",
        "The period is estimated only for cases that pass this activity test.",
        "",
        "## Scan results",
        "",
    ]

    for label, group in combined.groupby("scan", sort=False):
        group = group.sort_values("s0")
        last_inactive, first_active = onset_bracket(group)
        first_active_row = (
            group[group["s0"] == first_active].iloc[0]
            if not math.isnan(first_active)
            else None
        )
        lines.extend(
            [
                f"### {label}",
                "",
                f"- `eps = {float(group['eps_shift'].iloc[0]):.12g}`",
                f"- `dt = {float(group['dt'].iloc[0]):.12g}`",
                f"- sampled range: `{float(group['s0'].min()):.12g} <= s0 <= "
                f"{float(group['s0'].max()):.12g}` ({len(group)} cases)",
            ]
        )
        if math.isnan(first_active):
            lines.append("- no active sampled point")
        elif math.isnan(last_inactive):
            lines.append(f"- first active sampled point: `s0 = {first_active:.12g}`")
        else:
            lines.append(
                f"- sampled onset bracket: `{last_inactive:.12g} < s0 <= "
                f"{first_active:.12g}`"
            )
        if first_active_row is not None and math.isfinite(float(first_active_row["estimated_period"])):
            lines.append(
                f"- first-active nonlinear period: "
                f"`{float(first_active_row['estimated_period']):.12g}`"
            )
        contains = (
            not math.isnan(last_inactive)
            and not math.isnan(first_active)
            and last_inactive < imex["critical_s0"] <= first_active
        )
        lines.append(f"- IMEX threshold lies in sampled bracket: `{contains}`")
        lines.append("")

    output.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scan", action="append", type=parse_scan, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--hopf-csv", type=Path, required=True)
    parser.add_argument("--imex-csv", type=Path, required=True)
    parser.add_argument("--late-fraction", type=float, default=0.25)
    parser.add_argument("--active-vrms", type=float, default=1.0e-4)
    parser.add_argument("--active-amp", type=float, default=1.0e-4)
    parser.add_argument("--zoom-window", type=float, default=1.0)
    parser.add_argument(
        "--formats",
        nargs="+",
        default=["png", "eps"],
        choices=["png", "eps", "pdf", "svg"],
        help="Figure output formats (default: png eps)",
    )
    args = parser.parse_args()

    if not 0.0 < args.late_fraction <= 1.0:
        parser.error("--late-fraction must be in (0, 1]")

    out = args.output_dir.expanduser().resolve()
    fig_dir = out / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    frames: list[pd.DataFrame] = []
    for spec in args.scan:
        frames.append(
            analyze_scan(
                spec,
                late_fraction=args.late_fraction,
                active_vrms=args.active_vrms,
                active_amp=args.active_amp,
            )
        )

    combined = pd.concat(frames, ignore_index=True)
    target_dt = float(combined["dt"].iloc[0])
    if not np.allclose(combined["dt"].to_numpy(dtype=float), target_dt):
        raise ValueError("all scans must use the same dt")

    hopf = read_hopf_reference(args.hopf_csv.expanduser().resolve())
    imex = read_imex_reference(args.imex_csv.expanduser().resolve(), target_dt)

    combined.to_csv(out / "grid_aligned_summary.csv", index=False)
    plot_activity(
        combined,
        "late_delta_amp",
        r"late displacement amplitude $A^{\mathrm{late}}_{\delta x}$",
        fig_dir / "late_displacement_amplitude",
        hopf,
        imex,
        args.formats,
        show_hopf_reference=True,
    )
    plot_activity(
        combined,
        "late_vrms",
        r"late velocity activity $V^{\mathrm{late}}_{\mathrm{rms}}$",
        fig_dir / "late_velocity_rms",
        hopf,
        imex,
        args.formats,
        show_hopf_reference=False,
    )
    plot_period(
        combined, fig_dir / "active_periods", hopf, imex, args.formats
    )
    plot_brackets(
        combined,
        fig_dir / "onset_brackets_vs_epsilon",
        hopf,
        imex,
        args.formats,
    )

    for label, group in combined.groupby("scan", sort=False):
        safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in label)
        plot_transition_trajectories(
            group,
            fig_dir / f"transition_trajectories_{safe}",
            zoom_window=args.zoom_window,
            formats=args.formats,
        )

    write_summary(
        combined,
        out / "grid_aligned_summary.md",
        hopf,
        imex,
        active_vrms=args.active_vrms,
        active_amp=args.active_amp,
    )

    metadata = {
        "scans": [{"label": spec.label, "path": str(spec.path)} for spec in args.scan],
        "late_fraction": args.late_fraction,
        "active_vrms": args.active_vrms,
        "active_amp": args.active_amp,
        "zoom_window": args.zoom_window,
        "formats": args.formats,
        "hopf_reference": hopf,
        "imex_reference": imex,
        "n_cases": int(len(combined)),
    }
    (out / "analysis_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )

    print(f"Wrote {out / 'grid_aligned_summary.csv'}")
    print(f"Wrote {out / 'grid_aligned_summary.md'}")
    print(f"Figures: {fig_dir}")
    print(f"Formats: {', '.join(args.formats)}")


if __name__ == "__main__":
    main()
