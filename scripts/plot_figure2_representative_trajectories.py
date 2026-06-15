#!/usr/bin/env python3
import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# Journal-style font and saving settings
# ============================================================

FONT_FAMILY = "Arial"
BASE_FONT_SIZE = 16
TICK_FONT_SIZE = 14
TITLE_FONT_SIZE = 18


def configure_figure_fonts():
    """
    Apply global matplotlib font settings.

    Main text: Arial, 16 pt
    Tick labels: Arial, 14 pt
    Titles: Arial, 18 pt
    """
    plt.rcParams.update({
        "font.family": FONT_FAMILY,
        "font.sans-serif": [FONT_FAMILY, "Helvetica", "DejaVu Sans"],
        "font.size": BASE_FONT_SIZE,
        "axes.labelsize": BASE_FONT_SIZE,
        "axes.titlesize": TITLE_FONT_SIZE,
        "xtick.labelsize": TICK_FONT_SIZE,
        "ytick.labelsize": TICK_FONT_SIZE,
        "legend.fontsize": BASE_FONT_SIZE,
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


def apply_figure_style(fig):
    """
    Enforce font sizes after all labels/titles/legends have been created.
    """
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


def save_figure(fig, out, dpi=250):
    """
    Save both PNG and EPS using the same base filename.

    Example:
        out = figure/foo.png
    writes:
        figure/foo.png
        figure/foo.eps
    """
    from pathlib import Path

    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)

    apply_figure_style(fig)

    fig.savefig(out, dpi=dpi, bbox_inches="tight")

    eps_out = out.with_suffix(".eps")
    fig.savefig(eps_out, format="eps", dpi=dpi, bbox_inches="tight")

    print(f"saved: {out}")
    print(f"saved: {eps_out}")



# ============================================================
# Plot settings
# ============================================================

STATIONARY_WINDOW = 50.0
OSCILLATORY_WINDOW = 0.5
ROTATING_WINDOW = 0.5

LINE_WIDTH = 1.3
DPI = 250


# ============================================================
# Data readers
# ============================================================

def read_summary(path: Path):
    cases = []

    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)

        for row in reader:
            cases.append({
                "s0": float(row["s0"]),
                "label": row["label"],
                "folder": row["folder"],
            })

    return cases


def read_traj(path: Path):
    data = np.genfromtxt(path, delimiter=",", names=True)

    if data.ndim == 0:
        data = np.array([data])

    return data


# ============================================================
# State inference and time-window selection
# ============================================================

def infer_state(cas):
    """
    Infer the dynamical state from the label or from s0.

    If the C++ code labels cases as "case 1", "case 2", ...
    this function uses the current phase diagram thresholds.
    """
    label = str(cas["label"]).lower()
    s0 = float(cas["s0"])

    if "stationary" in label:
        return "stationary"

    if "oscillatory" in label:
        return "oscillatory"

    if "rotating" in label:
        return "rotating"

    # Based on the current phase diagram.
    if s0 < 0.10:
        return "stationary"
    elif s0 < 0.155:
        return "oscillatory"
    else:
        return "rotating"


def display_label(cas):
    return infer_state(cas)


def window_for_case(cas):
    state = infer_state(cas)

    if state == "stationary":
        return STATIONARY_WINDOW

    if state == "oscillatory":
        return OSCILLATORY_WINDOW

    if state == "rotating":
        return ROTATING_WINDOW

    return 5.0


def crop_tail(t, y, cas):
    """
    Keep only the final time window and shift time so that
    the plotted interval starts at t=0.
    """
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)

    window = window_for_case(cas)

    tmax = np.max(t)
    mask = t >= (tmax - window)

    tc = t[mask].copy()
    yc = y[mask].copy()

    if len(tc) > 0:
        tc -= tc[0]

    return tc, yc


# ============================================================
# Plot helpers
# ============================================================

def break_periodic_jumps(x, y, jump_threshold):
    """
    Insert NaN at discontinuities caused by periodic wrapping.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if len(y) <= 1:
        return x, y

    dy = np.abs(np.diff(y))
    jumps = np.where(dy > jump_threshold)[0]

    if len(jumps) == 0:
        return x, y

    xs = []
    ys = []
    last = 0

    for j in jumps:
        xs.extend(list(x[last:j + 1]))
        ys.extend(list(y[last:j + 1]))

        xs.append(np.nan)
        ys.append(np.nan)

        last = j + 1

    xs.extend(list(x[last:]))
    ys.extend(list(y[last:]))

    return np.array(xs), np.array(ys)


def set_theta_axis(ax):
    ax.set_ylim(-0.1, 2.1 * np.pi)
    ax.set_yticks([0.0, np.pi, 2.0 * np.pi])
    ax.set_yticklabels([r"$0$", r"$\pi$", r"$2\pi$"])
    ax.grid(True, alpha=0.25)


# ============================================================
# compact wrapped trajectories, final windows only
# ============================================================

def plot_compact(data_root: Path, fig_dir: Path, cases):
    fig, axes = plt.subplots(
        1,
        len(cases),
        figsize=(5.0 * len(cases), 3.7),
        sharey=True
    )

    if len(cases) == 1:
        axes = [axes]

    for ax, cas in zip(axes, cases):
        traj_path = data_root / cas["folder"] / "trajectory.csv"
        tr = read_traj(traj_path)

        t = np.asarray(tr["t"], dtype=float)
        theta = np.asarray(tr["theta_wrapped"], dtype=float)

        # Use only final representative behavior.
        t, theta = crop_tail(t, theta, cas)

        # Break periodic jumps.
        t_plot, theta_plot = break_periodic_jumps(
            t,
            theta,
            jump_threshold=np.pi
        )

        ax.plot(t_plot, theta_plot, color="black", lw=1.2)

        label = display_label(cas)

        ax.set_title(
            rf"{label}" + "\n" + rf"$s_0={cas['s0']:.3f}$"
        )

        ax.set_xlabel(r"$t-t_0$")
        set_theta_axis(ax)

    axes[0].set_ylabel(r"$\theta(t)$")

    fig.tight_layout()

    out = fig_dir / "figure2_representative_trajectories_compact.png"
    save_figure(fig, out, dpi=DPI)
    plt.close(fig)

    print(f"saved: {out}")


# ============================================================
# compact unwrapped view for rotating state
# ============================================================

def plot_rotating_unwrapped(data_root: Path, fig_dir: Path, cases):
    rotating_cases = [
        cas for cas in cases
        if infer_state(cas) == "rotating"
    ]

    if not rotating_cases:
        return

    fig, ax = plt.subplots(figsize=(6.2, 4.0))

    for cas in rotating_cases:
        traj_path = data_root / cas["folder"] / "trajectory.csv"
        tr = read_traj(traj_path)

        t = np.asarray(tr["t"], dtype=float)
        theta = np.asarray(tr["theta_unwrapped"], dtype=float)

        # For rotating state, use the last part but do not wrap.
        t, theta = crop_tail(t, theta, cas)

        # Shift theta for visual clarity.
        theta = theta - theta[0]

        ax.plot(
            t,
            theta,
            color="black",
            lw=1.5,
            label=rf"$s_0={cas['s0']:.3f}$"
        )

    ax.set_xlabel(r"$t-t_0$")
    ax.set_ylabel(r"$\theta(t)-\theta(t_0)$")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")

    fig.tight_layout()

    out = fig_dir / "figure2_rotating_unwrapped_final_window.png"
    save_figure(fig, out, dpi=DPI)
    plt.close(fig)

    print(f"saved: {out}")


# ============================================================
# Main
# ============================================================

def main():
    data_root = Path("figure2_data")
    fig_dir = Path("figure")

    if len(sys.argv) >= 2:
        data_root = Path(sys.argv[1])

    if len(sys.argv) >= 3:
        fig_dir = Path(sys.argv[2])

    fig_dir.mkdir(parents=True, exist_ok=True)

    summary_path = data_root / "representative_cases.csv"

    cases = read_summary(summary_path)

    plot_compact(data_root, fig_dir, cases)
    plot_rotating_unwrapped(data_root, fig_dir, cases)


if __name__ == "__main__":
    main()