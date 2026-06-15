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

YMIN = 0.0
YMAX = 1.5
DPI = 250

# Droplet radius used in the simulations.
# You can override from command line:
#   python3 scripts/plot_figure3_representative_profiles.py figure2_data figure 1.0
DEFAULT_RADIUS = 1.0


# ============================================================
# Robust readers
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


def read_two_column(path: Path):
    """
    Read files of the form:
        x value

    Broken lines are skipped.
    """
    if not path.exists():
        raise FileNotFoundError(path)

    xs = []
    ys = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            parts = line.split()

            if len(parts) < 2:
                continue

            try:
                x = float(parts[0])
                y = float(parts[1])
            except Exception:
                # Try to repair broken scientific notation, e.g. "2.3 e-06"
                if len(parts) >= 3:
                    try:
                        x = float(parts[0])
                        y = float(parts[1] + parts[2])
                    except Exception:
                        continue
                else:
                    continue

            xs.append(x)
            ys.append(y)

    if len(xs) == 0:
        raise RuntimeError(f"No valid two-column data in {path}")

    return np.column_stack([np.asarray(xs), np.asarray(ys)])


def read_droplet_position(path: Path):
    """
    x file generated by save_snapshot has:
        xc 0
    For N_BS=1, the first row gives droplet center.
    """
    arr = read_two_column(path)
    return float(arr[0, 0])


# ============================================================
# State labels
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


# ============================================================
# Coordinate conversion
# ============================================================

def estimate_domain_length(x):
    """
    Profile files contain x = ix * dx, ix=0,...,Nx-1.
    Since the endpoint L is not included, estimate L as max(x)+dx.
    """
    x = np.asarray(x, dtype=float)

    if len(x) < 2:
        return float(np.max(x))

    dx = np.median(np.diff(np.sort(x)))
    return float(np.max(x) + dx)


def periodic_relative_coordinate(x, xc, Lx):
    """
    Convert x to xi=x-xc on the periodic domain, mapped to [-Lx/2,Lx/2).
    """
    xi = x - xc

    xi = (xi + 0.5 * Lx) % Lx - 0.5 * Lx

    return xi


def load_profiles_in_comoving_frame(case_dir: Path):
    u = read_two_column(case_dir / "u")
    s = read_two_column(case_dir / "s")
    c = read_two_column(case_dir / "c")
    w = read_two_column(case_dir / "w")
    xc = read_droplet_position(case_dir / "x")

    x = u[:, 0]
    Lx = estimate_domain_length(x)
    xi = periodic_relative_coordinate(x, xc, Lx)

    order = np.argsort(xi)

    profiles = {
        "xi": xi[order],
        "u": u[:, 1][order],
        "s": s[:, 1][order],
        "c": c[:, 1][order],
        "w": w[:, 1][order],
        "xc": xc,
        "Lx": Lx,
    }

    return profiles


# ============================================================
# Plotting
# ============================================================

def plot_case_profiles(ax, cas, profiles, radius):
    xi = profiles["xi"]

    ax.axvspan(-radius, radius, color="0.90", zorder=0)
    ax.axvline(-radius, color="black", linestyle="--", linewidth=0.9, alpha=0.6)
    ax.axvline(radius, color="black", linestyle="--", linewidth=0.9, alpha=0.6)

    ax.plot(xi, profiles["u"], color="C0", lw=2.0, label=r"$u/u_0$")
    ax.plot(xi, profiles["s"], color="C1", lw=2.0, label=r"$s/s_0$")
    ax.plot(xi, profiles["c"], color="C2", lw=2.0, label=r"$c/c_0$")
    ax.plot(xi, profiles["w"], color="C3", lw=2.0, label=r"$\bar W/w_0$")

    state = infer_state(cas)

    ax.set_title(
        rf"{state}" + "\n" + rf"$s_0={cas['s0']:.3f}$",
        fontsize=13
    )

    ax.set_xlabel(r"$\xi=x-x_c$")
    ax.set_ylim(YMIN, YMAX)
    ax.grid(True, alpha=0.25)


def plot_profiles_three_columns(data_root: Path, fig_dir: Path, cases, radius):
    fig, axes = plt.subplots(
        1,
        len(cases),
        figsize=(5.2 * len(cases), 4.0),
        sharey=True
    )

    if len(cases) == 1:
        axes = [axes]

    for ax, cas in zip(axes, cases):
        case_dir = data_root / cas["folder"]
        profiles = load_profiles_in_comoving_frame(case_dir)
        plot_case_profiles(ax, cas, profiles, radius)

    axes[0].set_ylabel("normalized concentration")

    # Put legend only once.
    handles, labels = axes[-1].get_legend_handles_labels()
    axes[-1].legend(handles, labels, loc="upper right", fontsize=9)

    fig.tight_layout()

    out = fig_dir / "figure3_representative_profiles_compact.png"
    save_figure(fig, out, dpi=DPI)
    plt.close(fig)

    print(f"saved: {out}")


def plot_s_only_zoom(data_root: Path, fig_dir: Path, cases, radius):
    """
    Optional zoomed plot of SDS profile.
    This is often useful because s/s0 varies only slightly near 1.
    """
    fig, axes = plt.subplots(
        1,
        len(cases),
        figsize=(5.2 * len(cases), 3.7),
        sharey=False
    )

    if len(cases) == 1:
        axes = [axes]

    for ax, cas in zip(axes, cases):
        case_dir = data_root / cas["folder"]
        profiles = load_profiles_in_comoving_frame(case_dir)

        xi = profiles["xi"]
        s = profiles["s"]

        ax.axvspan(-radius, radius, color="0.90", zorder=0)
        ax.axvline(-radius, color="black", linestyle="--", linewidth=0.9, alpha=0.6)
        ax.axvline(radius, color="black", linestyle="--", linewidth=0.9, alpha=0.6)

        ax.plot(xi, s, color="C1", lw=2.0)

        state = infer_state(cas)
        ax.set_title(
            rf"{state}" + "\n" + rf"$s_0={cas['s0']:.3f}$",
            fontsize=13
        )

        ax.set_xlabel(r"$\xi=x-x_c$")
        ax.set_ylabel(r"$s/s_0$")
        ax.grid(True, alpha=0.25)

        ymin = np.min(s)
        ymax = np.max(s)
        margin = 0.10 * max(ymax - ymin, 1e-6)
        ax.set_ylim(ymin - margin, ymax + margin)

    fig.tight_layout()

    out = fig_dir / "figure3_sds_profiles_zoom.png"
    save_figure(fig, out, dpi=DPI)
    plt.close(fig)

    print(f"saved: {out}")


# ============================================================
# Main
# ============================================================

def main():
    data_root = Path("figure2_data")
    fig_dir = Path("figure")
    radius = DEFAULT_RADIUS

    if len(sys.argv) >= 2:
        data_root = Path(sys.argv[1])

    if len(sys.argv) >= 3:
        fig_dir = Path(sys.argv[2])

    if len(sys.argv) >= 4:
        radius = float(sys.argv[3])

    fig_dir.mkdir(parents=True, exist_ok=True)

    summary_path = data_root / "representative_cases.csv"
    cases = read_summary(summary_path)

    plot_profiles_three_columns(data_root, fig_dir, cases, radius)
    plot_s_only_zoom(data_root, fig_dir, cases, radius)


if __name__ == "__main__":
    main()
