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
# Settings
# ============================================================

DPI = 250

STATIONARY_WINDOW = 50.0
OSCILLATORY_WINDOW = 0.5
ROTATING_WINDOW = 0.5

# Common y-range for the bottom row of the compact figure.
# Adjust here if needed.
TENSION_YMIN = -0.002
TENSION_YMAX = 0.0025


# ============================================================
# Readers
# ============================================================

def read_cases(path: Path):
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


def read_csv_named(path: Path):
    data = np.genfromtxt(path, delimiter=",", names=True)

    if data.ndim == 0:
        data = np.array([data])

    return data


# ============================================================
# State and time-window helpers
# ============================================================

def infer_state(cas):
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


def window_for_case(cas):
    state = infer_state(cas)

    if state == "stationary":
        return STATIONARY_WINDOW

    if state == "oscillatory":
        return OSCILLATORY_WINDOW

    if state == "rotating":
        return ROTATING_WINDOW

    return 1.0


def crop_tail(t, arrays, cas):
    """
    Keep only the final time window and shift time so that
    the plotted interval starts at t=0.
    """
    t = np.asarray(t, dtype=float)

    tmax = np.max(t)
    window = window_for_case(cas)
    mask = t >= (tmax - window)

    tc = t[mask].copy()

    if len(tc) > 0:
        tc -= tc[0]

    cropped = [
        np.asarray(a, dtype=float)[mask].copy()
        for a in arrays
    ]

    return tc, cropped


def break_periodic_jumps(x, y, jump_threshold=np.pi):
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
    ax.set_yticks([0, np.pi, 2 * np.pi])
    ax.set_yticklabels([r"$0$", r"$\pi$", r"$2\pi$"])
    ax.grid(True, alpha=0.25)


# ============================================================
# Main compact figure: theta and surface-tension components
# ============================================================

def plot_compact_decomposition(data_root: Path, fig_dir: Path, cases):
    fig, axes = plt.subplots(
        2,
        len(cases),
        figsize=(5.2 * len(cases), 6.2),
        sharex=False
    )

    if len(cases) == 1:
        axes = np.array([[axes[0]], [axes[1]]])

    for j, cas in enumerate(cases):
        data_path = data_root / cas["folder"] / "surface_tension.csv"
        data = read_csv_named(data_path)

        t = np.asarray(data["t"], dtype=float)
        theta_wrapped = np.asarray(data["theta_wrapped"], dtype=float)

        dg = np.asarray(data["delta_g_norm"], dtype=float)
        dh = np.asarray(data["delta_h_norm"], dtype=float)
        dgamma = np.asarray(data["delta_gamma_norm"], dtype=float)

        t_crop, (theta_crop, dg_crop, dh_crop, dgamma_crop) = crop_tail(
            t,
            [theta_wrapped, dg, dh, dgamma],
            cas
        )

        t_theta, theta_plot = break_periodic_jumps(
            t_crop,
            theta_crop,
            jump_threshold=np.pi
        )

        state = infer_state(cas)

        ax0 = axes[0, j]
        ax1 = axes[1, j]

        # Top row: trajectory
        ax0.plot(t_theta, theta_plot, color="black", lw=1.2)

        ax0.set_title(
            rf"{state}" + "\n" + rf"$s_0={cas['s0']:.3f}$",
            fontsize=13
        )

        ax0.set_xlabel(r"$t-t_0$")
        ax0.set_ylabel(r"$\theta(t)$")
        set_theta_axis(ax0)

        # Bottom row: surface-tension imbalance
        ax1.axhline(0.0, color="black", lw=0.8, alpha=0.5)

        ax1.plot(
            t_crop,
            dg_crop,
            lw=1.4,
            label=r"$\Delta g/G_0$"
        )

        ax1.plot(
            t_crop,
            dh_crop,
            lw=1.4,
            label=r"$\Delta h/G_0$"
        )

        ax1.plot(
            t_crop,
            dgamma_crop,
            lw=1.8,
            label=r"$\Delta\gamma/G_0$"
        )

        ax1.set_xlabel(r"$t-t_0$")
        ax1.set_ylabel(r"normalized tension diff.")

        # Common y-axis range for all bottom-row panels.
        ax1.set_ylim(TENSION_YMIN, TENSION_YMAX)

        ax1.grid(True, alpha=0.25)

        if j == len(cases) - 1:
            ax1.legend(loc="best", fontsize=8)

    fig.tight_layout()

    out = fig_dir / "figure4_surface_tension_decomposition_compact.png"
    save_figure(fig, out, dpi=DPI)
    plt.close(fig)

    print(f"saved: {out}")


# ============================================================
# Oscillatory detail figure
# ============================================================

def plot_oscillatory_detail(data_root: Path, fig_dir: Path, cases):
    osc_cases = [
        cas for cas in cases
        if infer_state(cas) == "oscillatory"
    ]

    if not osc_cases:
        return

    cas = osc_cases[0]

    data_path = data_root / cas["folder"] / "surface_tension.csv"
    data = read_csv_named(data_path)

    t = np.asarray(data["t"], dtype=float)

    theta = np.asarray(data["theta_wrapped"], dtype=float)
    velocity = np.asarray(data["velocity"], dtype=float)

    dg = np.asarray(data["delta_g_norm"], dtype=float)
    dh = np.asarray(data["delta_h_norm"], dtype=float)
    dgamma = np.asarray(data["delta_gamma_norm"], dtype=float)

    t_crop, (theta, velocity, dg, dh, dgamma) = crop_tail(
        t,
        [theta, velocity, dg, dh, dgamma],
        cas
    )

    t_theta, theta_plot = break_periodic_jumps(
        t_crop,
        theta,
        jump_threshold=np.pi
    )

    fig, axes = plt.subplots(
        3,
        1,
        figsize=(8.5, 7.0),
        sharex=True
    )

    axes[0].plot(t_theta, theta_plot, color="black", lw=1.3)
    axes[0].set_ylabel(r"$\theta(t)$")
    axes[0].set_yticks([0, np.pi, 2 * np.pi])
    axes[0].set_yticklabels([r"$0$", r"$\pi$", r"$2\pi$"])
    axes[0].grid(True, alpha=0.25)

    axes[1].axhline(0.0, color="black", lw=0.8, alpha=0.5)
    axes[1].plot(t_crop, velocity, color="black", lw=1.3)
    axes[1].set_ylabel(r"$\dot{x}_c(t)$")
    axes[1].grid(True, alpha=0.25)

    axes[2].axhline(0.0, color="black", lw=0.8, alpha=0.5)

    axes[2].plot(
        t_crop,
        dg,
        lw=1.4,
        label=r"$\Delta g/G_0$"
    )

    axes[2].plot(
        t_crop,
        dh,
        lw=1.4,
        label=r"$\Delta h/G_0$"
    )

    axes[2].plot(
        t_crop,
        dgamma,
        lw=1.8,
        label=r"$\Delta\gamma/G_0$"
    )

    axes[2].set_ylabel(r"normalized tension diff.")
    axes[2].set_xlabel(r"$t-t_0$")
    axes[2].grid(True, alpha=0.25)
    axes[2].legend(loc="best", fontsize=9)

    fig.tight_layout()

    out = fig_dir / "figure4_oscillatory_surface_tension_detail.png"
    save_figure(fig, out, dpi=DPI)
    plt.close(fig)

    print(f"saved: {out}")


# ============================================================
# Rotating detail figure
# ============================================================

def plot_rotating_detail(data_root: Path, fig_dir: Path, cases):
    rot_cases = [
        cas for cas in cases
        if infer_state(cas) == "rotating"
    ]

    if not rot_cases:
        return

    cas = rot_cases[0]

    data_path = data_root / cas["folder"] / "surface_tension.csv"
    data = read_csv_named(data_path)

    t = np.asarray(data["t"], dtype=float)

    theta_unwrapped = np.asarray(data["theta_unwrapped"], dtype=float)

    dg = np.asarray(data["delta_g_norm"], dtype=float)
    dh = np.asarray(data["delta_h_norm"], dtype=float)
    dgamma = np.asarray(data["delta_gamma_norm"], dtype=float)

    t_crop, (theta_unwrapped, dg, dh, dgamma) = crop_tail(
        t,
        [theta_unwrapped, dg, dh, dgamma],
        cas
    )

    theta_unwrapped = theta_unwrapped - theta_unwrapped[0]

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(8.5, 5.5),
        sharex=True
    )

    axes[0].plot(t_crop, theta_unwrapped, color="black", lw=1.4)
    axes[0].set_ylabel(r"$\theta(t)-\theta(t_0)$")
    axes[0].grid(True, alpha=0.25)

    axes[1].axhline(0.0, color="black", lw=0.8, alpha=0.5)

    axes[1].plot(
        t_crop,
        dg,
        lw=1.4,
        label=r"$\Delta g/G_0$"
    )

    axes[1].plot(
        t_crop,
        dh,
        lw=1.4,
        label=r"$\Delta h/G_0$"
    )

    axes[1].plot(
        t_crop,
        dgamma,
        lw=1.8,
        label=r"$\Delta\gamma/G_0$"
    )

    axes[1].set_ylabel(r"normalized tension diff.")
    axes[1].set_xlabel(r"$t-t_0$")
    axes[1].grid(True, alpha=0.25)
    axes[1].legend(loc="best", fontsize=9)

    fig.tight_layout()

    out = fig_dir / "figure4_rotating_surface_tension_detail.png"
    save_figure(fig, out, dpi=DPI)
    plt.close(fig)

    print(f"saved: {out}")


# ============================================================
# Main
# ============================================================

def main():
    data_root = Path("figure4_data")
    fig_dir = Path("figure")

    if len(sys.argv) >= 2:
        data_root = Path(sys.argv[1])

    if len(sys.argv) >= 3:
        fig_dir = Path(sys.argv[2])

    fig_dir.mkdir(parents=True, exist_ok=True)

    cases_path = data_root / "surface_tension_cases.csv"
    cases = read_cases(cases_path)

    plot_compact_decomposition(data_root, fig_dir, cases)
    plot_oscillatory_detail(data_root, fig_dir, cases)
    plot_rotating_detail(data_root, fig_dir, cases)


if __name__ == "__main__":
    main()