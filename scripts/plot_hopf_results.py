#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent

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


def setup_grid(ax: plt.Axes) -> None:
    """Use an opaque grid so the EPS backend does not need transparency."""
    ax.grid(True, color="0.85", linewidth=0.8)


def save_figure(fig: plt.Figure, output_stem: Path, formats: list[str]) -> None:
    """Save one Matplotlib figure in all requested formats."""
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    apply_figure_style(fig)
    for fmt in formats:
        fmt = fmt.lower().lstrip(".")
        output = output_stem.with_suffix(f".{fmt}")
        fig.savefig(
            output,
            format=fmt,
            dpi=300 if fmt == "png" else None,
            bbox_inches="tight",
            pad_inches=0.10,
            facecolor="white",
            transparent=False,
        )
        print(f"Wrote {output}")


def save_line_plot(
    data_sets,
    x,
    y,
    xlabel,
    ylabel,
    output_stem,
    formats,
    horizontal_zero=False,
    verticals=None,
):
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for label, frame in data_sets:
        ax.plot(
            frame[x],
            frame[y],
            marker="o",
            markersize=3,
            linewidth=1.2,
            label=label,
        )
    if horizontal_zero:
        ax.axhline(0.0, linewidth=0.8)
    if verticals:
        for value in verticals:
            ax.axvline(value, linestyle="--", linewidth=0.8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    setup_grid(ax)
    if len(data_sets) > 1:
        ax.legend()
    fig.tight_layout()
    save_figure(fig, output_stem, formats)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, default=ROOT / "results")
    parser.add_argument(
        "--formats",
        nargs="+",
        default=["png", "eps"],
        choices=["png", "eps", "pdf", "svg"],
        help="figure formats to write (default: png eps)",
    )
    args = parser.parse_args()

    results = args.results.expanduser().resolve()
    figures = results / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    critical = pd.read_csv(results / "critical_points.csv")
    scans = []
    for nx in critical.nx.astype(int):
        scans.append(
            (f"Nx={nx}", pd.read_csv(results / f"Nx{nx}" / "eigenvalue_scan.csv"))
        )
    verticals = critical.s_h.tolist()

    save_line_plot(
        scans,
        "s0",
        "leading_real",
        r"$s_0$",
        r"$\mathrm{Re}\,\lambda_+$",
        figures / "leading_real_vs_s0",
        args.formats,
        horizontal_zero=True,
        verticals=verticals,
    )
    save_line_plot(
        scans,
        "s0",
        "leading_frequency_abs",
        r"$s_0$",
        r"$|\mathrm{Im}\,\lambda_+|$",
        figures / "leading_imag_vs_s0",
        args.formats,
        verticals=verticals,
    )

    nx = int(critical.iloc[0].nx)
    response = pd.read_csv(results / f"Nx{nx}" / "response_curve.csv")
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.plot(
        response.omega,
        response.H_real_minus_1,
        label=r"$\mathrm{Re}\,H(i\omega)-1$",
    )
    ax.plot(response.omega, response.H_imag, label=r"$\mathrm{Im}\,H(i\omega)$")
    ax.axhline(0.0, linewidth=0.8)
    ax.axvline(float(critical.iloc[0].omega_h), linestyle="--", linewidth=0.8)
    ax.set_xlabel(r"$\omega$")
    ax.set_ylabel("response residual")
    setup_grid(ax)
    ax.legend()
    fig.tight_layout()
    save_figure(fig, figures / "response_residual_vs_omega", args.formats)
    plt.close(fig)

    spectrum = pd.read_csv(results / f"Nx{nx}" / "rightmost_spectrum.csv")
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    # Spectrum points are discrete; do not connect them with lines.
    ax.scatter(spectrum.real, spectrum.imag, s=28)
    ax.axvline(0.0, linewidth=0.8)
    ax.axhline(0.0, linewidth=0.8)
    ax.set_xlabel(r"$\mathrm{Re}\,\lambda$")
    ax.set_ylabel(r"$\mathrm{Im}\,\lambda$")
    setup_grid(ax)
    fig.tight_layout()
    save_figure(fig, figures / "critical_spectrum", args.formats)
    plt.close(fig)

    profile = pd.read_csv(results / f"Nx{nx}" / "steady_profile_critical.csv")
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.plot(profile.xi, profile.u / 0.006, label=r"$U_0/u_0$")
    ax.plot(profile.xi, profile.s / float(critical.iloc[0].s_h), label=r"$S_0/s_H$")
    ax.plot(profile.xi, profile.c / 0.02, label=r"$C_0/c_0$")
    ax.set_xlim(-4.0, 4.0)
    ax.set_xlabel(r"$\xi$")
    ax.set_ylabel("normalized concentration")
    setup_grid(ax)
    ax.legend()
    fig.tight_layout()
    save_figure(fig, figures / "critical_steady_profile", args.formats)
    plt.close(fig)

    time_step_file = results / "time_step_stability.csv"
    metadata_file = results / "run_metadata.json"
    if time_step_file.exists() and metadata_file.exists():
        time_step = pd.read_csv(time_step_file)
        metadata = json.loads(metadata_file.read_text(encoding="utf-8"))
        imex = time_step[time_step.dt > 0.0].sort_values("dt")
        continuous = time_step[time_step.dt == 0.0]
        fig, ax = plt.subplots(figsize=(6.4, 4.2))
        ax.plot(imex.dt, imex.critical_s0, marker="o", label="production IMEX map")
        if not continuous.empty:
            ax.axhline(
                float(continuous.iloc[0].critical_s0),
                linestyle="--",
                label="continuous semi-discrete",
            )
        ax.axhline(
            float(metadata["second_order_richardson_s_h"]),
            linestyle=":",
            label=r"$H(\lambda)=1$ spatial extrapolation",
        )
        ax.set_xscale("log")
        ax.set_xlabel(r"time step $\Delta t$")
        ax.set_ylabel(r"critical $s_0$")
        setup_grid(ax)
        ax.legend()
        fig.tight_layout()
        save_figure(fig, figures / "time_step_shift_of_critical_s0", args.formats)
        plt.close(fig)

    trajectory_file = results / "nonlinear_dt_check" / "trajectories_downsampled.csv"
    if trajectory_file.exists():
        trajectories = pd.read_csv(trajectory_file)
        fig, ax = plt.subplots(figsize=(6.4, 4.2))
        for (label, dt), frame in trajectories.groupby(["label", "dt"]):
            ax.plot(frame.t, frame.delta_x, label=f"{label}: dt={dt:g}")
        ax.set_xlabel(r"$t$")
        ax.set_ylabel(r"$\delta x(t)$")
        setup_grid(ax)
        ax.legend()
        fig.tight_layout()
        save_figure(fig, figures / "nonlinear_dt_check_displacement", args.formats)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(6.4, 4.2))
        for (label, dt), frame in trajectories.groupby(["label", "dt"]):
            ax.plot(frame.t, frame.velocity, label=f"{label}: dt={dt:g}")
        ax.set_xlabel(r"$t$")
        ax.set_ylabel("droplet velocity")
        setup_grid(ax)
        ax.legend()
        fig.tight_layout()
        save_figure(fig, figures / "nonlinear_dt_check_velocity", args.formats)
        plt.close(fig)


if __name__ == "__main__":
    main()
