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
DPI = 250


def configure_figure_fonts():
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


def save_figure(fig, out, dpi=DPI):
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)

    apply_figure_style(fig)

    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    eps_out = out.with_suffix(".eps")
    fig.savefig(eps_out, format="eps", dpi=dpi, bbox_inches="tight")

    print(f"saved: {out}")
    print(f"saved: {eps_out}")


# ============================================================
# Data loading
# ============================================================


def read_phase_csv(path: Path):
    rows = []

    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)

    if not rows:
        raise RuntimeError(f"No rows found in {path}")

    s0 = np.array([float(r["s0"]) for r in rows], dtype=float)
    amp = np.array([float(r["amplitude_theta"]) for r in rows], dtype=float)

    if "period" in rows[0]:
        period = np.array([float(r["period"]) for r in rows], dtype=float)
    else:
        period = np.full_like(s0, np.nan, dtype=float)

    if "state" in rows[0]:
        state = np.array([str(r["state"]) for r in rows], dtype=str)
    else:
        state = infer_state_from_amplitude(amp)

    order = np.argsort(s0)

    return s0[order], amp[order], period[order], state[order]


def infer_state_from_amplitude(amp):
    state = []
    for a in amp:
        if a < 0.05:
            state.append("stationary")
        elif a > 1.95 * np.pi:
            state.append("rotating")
        else:
            state.append("oscillatory")
    return np.array(state, dtype=str)


# ============================================================
# Plot helpers
# ============================================================


def state_regions(s0, state):
    regions = []
    start = 0

    for i in range(1, len(state)):
        if state[i] != state[i - 1]:
            regions.append((start, i - 1, state[i - 1]))
            start = i

    regions.append((start, len(state) - 1, state[-1]))
    return regions


def state_boundaries(s0, state):
    out = []
    for i in range(1, len(state)):
        if state[i] != state[i - 1]:
            out.append(0.5 * (s0[i] + s0[i - 1]))
    return out


def add_region_labels(ax, s0, state):
    y_text = 1.72 * np.pi

    for i0, i1, st in state_regions(s0, state):
        xm = 0.5 * (s0[i0] + s0[i1])
        ax.text(xm, y_text, st, ha="center", va="center")

    for xb in state_boundaries(s0, state):
        ax.axvline(xb, color="black", linestyle=":", linewidth=1.0)


# ============================================================
# Main plot
# ============================================================


def plot_phase_diagram(infile: Path, outfile: Path):
    s0, amp, period, state = read_phase_csv(infile)

    fig, ax = plt.subplots(figsize=(8.5, 4.6))

    ax.plot(
        s0,
        amp,
        "o",
        color="black",
        markersize=4.5,
        linewidth=1.2,
    )

    ax.axhline(np.pi, color="black", linewidth=0.8, alpha=0.35)
    ax.axhline(2.0 * np.pi, color="black", linewidth=0.8, alpha=0.35)

    add_region_labels(ax, s0, state)

    ax.set_xlabel(r"$s_0$")
    ax.set_ylabel(r"amplitude $\Delta\theta$")
    ax.set_ylim(-0.10, 2.12 * np.pi)

    ax.set_yticks([0.0, np.pi, 2.0 * np.pi])
    ax.set_yticklabels([r"$0$", r"$\pi$", r"$2\pi$"])

    ax.grid(True, alpha=0.25)

    fig.tight_layout()
    save_figure(fig, outfile, dpi=DPI)
    plt.close(fig)


# ============================================================
# Main
# ============================================================


def main():
    infile = Path("results/phase_diagram_long.csv")
    outfile = Path("figure/figure1_phase_diagram.png")

    if len(sys.argv) >= 2:
        infile = Path(sys.argv[1])

    if len(sys.argv) >= 3:
        outarg = Path(sys.argv[2])
        if outarg.suffix.lower() in [".png", ".eps", ".pdf", ".svg"]:
            outfile = outarg.with_suffix(".png")
        else:
            outfile = outarg / "figure1_phase_diagram.png"

    plot_phase_diagram(infile, outfile)


if __name__ == "__main__":
    main()
