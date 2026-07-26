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


def save_figure(fig, out, dpi=250, aliases=None):
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

    targets = [out]
    if aliases:
        targets.extend(Path(a) for a in aliases)

    for target in targets:
        target.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(target, dpi=dpi, bbox_inches="tight")
        eps_out = target.with_suffix(".eps")
        fig.savefig(eps_out, format="eps", dpi=dpi, bbox_inches="tight")
        print(f"saved: {target}")
        print(f"saved: {eps_out}")



# ============================================================
# Settings
# ============================================================

DPI = 250
DEFAULT_RADIUS = 1.0

# Activity diagnostics
ACTIVITY_WINDOW = 4.0
VRMS_WINDOW = 2.0
SMOOTH_PTS = 11


# ============================================================
# Basic readers
# ============================================================

def read_info(path: Path):
    info = {}

    if not path.exists():
        return info

    with open(path, "r") as f:
        for line in f:
            parts = line.strip().split()

            if len(parts) >= 2:
                key = parts[0]
                val = parts[1]

                try:
                    info[key] = float(val)
                except Exception:
                    info[key] = val

    return info


def read_cover_protocol(path: Path):
    """
    Read figure5_data/cover_protocol.csv.

    Expected columns:
        t_global,t_local,segment,x_wrapped,x_unwrapped,
        theta_wrapped,theta_unwrapped,velocity,Wbar
    """
    rows = []

    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)

        for row in reader:
            rows.append(row)

    if not rows:
        raise RuntimeError(f"No data rows in {path}")

    data = {
        "t_global": np.array([float(r["t_global"]) for r in rows]),
        "t_local": np.array([float(r["t_local"]) for r in rows]),
        "segment": np.array([r["segment"] for r in rows]),
        "x_wrapped": np.array([float(r["x_wrapped"]) for r in rows]),
        "x_unwrapped": np.array([float(r["x_unwrapped"]) for r in rows]),
        "theta_wrapped": np.array([float(r["theta_wrapped"]) for r in rows]),
        "theta_unwrapped": np.array([float(r["theta_unwrapped"]) for r in rows]),
        "velocity": np.array([float(r["velocity"]) for r in rows]),
        "Wbar": np.array([float(r["Wbar"]) for r in rows]),
    }

    return data


def read_two_column(path: Path):
    """
    Robust reader for files with two columns:
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
                # Try to repair broken scientific notation:
                # e.g. "2.64159 e-06"
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

    if not xs:
        raise RuntimeError(f"No valid data in {path}")

    return np.column_stack([np.asarray(xs), np.asarray(ys)])


# ============================================================
# Numerical diagnostics
# ============================================================

def moving_average(y, k):
    if k <= 1:
        return np.asarray(y, dtype=float).copy()

    k = int(k)
    pad = k // 2

    y = np.asarray(y, dtype=float)
    ypad = np.pad(y, (pad, pad), mode="edge")
    kernel = np.ones(k) / float(k)

    out = np.convolve(ypad, kernel, mode="valid")

    # If k is even, convolve length may differ by one.
    if len(out) > len(y):
        out = out[:len(y)]
    elif len(out) < len(y):
        out = np.pad(out, (0, len(y) - len(out)), mode="edge")

    return out


def sliding_amplitude(t, y, window):
    """
    A(t_i) = max(y)-min(y) over a sliding window centered at t_i.
    """
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)

    n = len(t)
    out = np.zeros(n)

    half = 0.5 * window
    j0 = 0
    j1 = 0

    for i in range(n):
        left = t[i] - half
        right = t[i] + half

        while j0 < n and t[j0] < left:
            j0 += 1

        while j1 < n and t[j1] <= right:
            j1 += 1

        seg = y[j0:j1]

        if len(seg) >= 2:
            out[i] = np.max(seg) - np.min(seg)
        else:
            out[i] = 0.0

    return out


def sliding_rms(t, y, window):
    """
    RMS over a sliding window centered at t_i.
    """
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)

    n = len(t)
    out = np.zeros(n)

    half = 0.5 * window
    j0 = 0
    j1 = 0

    for i in range(n):
        left = t[i] - half
        right = t[i] + half

        while j0 < n and t[j0] < left:
            j0 += 1

        while j1 < n and t[j1] <= right:
            j1 += 1

        seg = y[j0:j1]

        if len(seg) >= 2:
            out[i] = np.sqrt(np.mean(seg * seg))
        else:
            out[i] = 0.0

    return out


# ============================================================
# Cover segment helpers
# ============================================================

def segment_bounds(t, segment):
    """
    Return list of intervals:
        [(t_start, t_end, segment_name), ...]
    """
    bounds = []

    start = 0

    for i in range(1, len(segment)):
        if segment[i] != segment[i - 1]:
            bounds.append((t[start], t[i - 1], segment[i - 1]))
            start = i

    bounds.append((t[start], t[-1], segment[-1]))

    return bounds


def shade_protocol_regions(ax, bounds):
    colors = {
        "open": "0.94",
        "closed": "0.86",
        "reopened": "0.94",
    }

    for t0, t1, name in bounds:
        ax.axvspan(
            t0,
            t1,
            color=colors.get(name, "0.92"),
            zorder=0
        )

    for _, t1, _ in bounds[:-1]:
        ax.axvline(
            t1,
            color="black",
            linestyle=":",
            linewidth=1.2,
            zorder=3
        )


def add_protocol_labels(ax, bounds, y_frac=0.88):
    ymin, ymax = ax.get_ylim()
    y = ymin + y_frac * (ymax - ymin)

    for t0, t1, name in bounds:
        ax.text(
            0.5 * (t0 + t1),
            y,
            name,
            ha="center",
            va="center",
            fontsize=13
        )


# ============================================================
# Main clean timeseries figure
# ============================================================

def plot_cover_protocol_clean(data_root: Path, fig_dir: Path):
    csv_path = data_root / "cover_protocol.csv"
    info_path = data_root / "cover_protocol_info.txt"

    data = read_cover_protocol(csv_path)
    info = read_info(info_path)

    t = data["t_global"]
    segment = data["segment"]

    theta = np.mod(data["theta_wrapped"], 2.0 * np.pi)
    wbar = data["Wbar"]
    velocity = data["velocity"]

    bounds = segment_bounds(t, segment)

    # Activity diagnostics
    dtheta = sliding_amplitude(t, theta, ACTIVITY_WINDOW)
    dtheta = moving_average(dtheta, SMOOTH_PTS)

    vrms = sliding_rms(t, velocity, VRMS_WINDOW)
    vrms = moving_average(vrms, SMOOTH_PTS)

    fig, axes = plt.subplots(
        3,
        1,
        figsize=(10.5, 8.0),
        sharex=True
    )

    # --------------------------------------------------------
    # Panel 1: angular activity
    # --------------------------------------------------------
    ax = axes[0]
    shade_protocol_regions(ax, bounds)
    ax.plot(t, dtheta, color="black", lw=2.0, zorder=4)

    ax.set_ylabel(r"$\Delta\theta(t)$")
    ax.set_ylim(-0.05, 2.1 * np.pi)
    ax.set_yticks([0, np.pi, 2 * np.pi])
    ax.set_yticklabels([r"$0$", r"$\pi$", r"$2\pi$"])
    ax.grid(True, alpha=0.25)

    # --------------------------------------------------------
    # Panel 2: vapor
    # --------------------------------------------------------
    ax = axes[1]
    shade_protocol_regions(ax, bounds)
    ax.plot(t, wbar, color="tab:red", lw=2.2, zorder=4)

    ax.set_ylabel(r"$\bar{W}(t)/w_0$")
    ax.grid(True, alpha=0.25)

    # --------------------------------------------------------
    # Panel 3: velocity RMS
    # --------------------------------------------------------
    ax = axes[2]
    shade_protocol_regions(ax, bounds)
    ax.plot(t, vrms, color="black", lw=2.0, zorder=4)

    ax.set_ylabel(r"$V_{\rm rms}(t)$")
    ax.set_xlabel(r"$t$")
    ax.grid(True, alpha=0.25)

    for ax in axes:
        add_protocol_labels(ax, bounds, y_frac=0.88)

    fig.tight_layout()

    out = fig_dir / "figure5_cover_protocol_clean.png"
    save_figure(fig, out, dpi=DPI)
    plt.close(fig)

    print(f"saved: {out}")


# ============================================================
# Profile comparison at segment ends
# ============================================================

def estimate_domain_length(x):
    x = np.asarray(x, dtype=float)

    if len(x) < 2:
        return float(np.max(x))

    dx = np.median(np.diff(np.sort(x)))

    return float(np.max(x) + dx)


def periodic_relative_coordinate(x, xc, Lx):
    xi = x - xc
    xi = (xi + 0.5 * Lx) % Lx - 0.5 * Lx
    return xi


def load_profiles_in_comoving_frame(folder: Path):
    u = read_two_column(folder / "u")
    s = read_two_column(folder / "s")
    c = read_two_column(folder / "c")
    w = read_two_column(folder / "w")
    xfile = read_two_column(folder / "x")

    xc = float(xfile[0, 0])

    x = u[:, 0]
    Lx = estimate_domain_length(x)
    xi = periodic_relative_coordinate(x, xc, Lx)

    order = np.argsort(xi)

    return {
        "xi": xi[order],
        "u": u[:, 1][order],
        "s": s[:, 1][order],
        "c": c[:, 1][order],
        "w": w[:, 1][order],
        "xc": xc,
        "Lx": Lx,
    }


def plot_case_profile(ax, profiles, title, radius):
    xi = profiles["xi"]

    ax.axvspan(
        -radius,
        radius,
        color="0.90",
        zorder=0
    )

    ax.axvline(
        -radius,
        color="black",
        linestyle="--",
        linewidth=0.9,
        alpha=0.6
    )

    ax.axvline(
        radius,
        color="black",
        linestyle="--",
        linewidth=0.9,
        alpha=0.6
    )

    ax.plot(
        xi,
        profiles["u"],
        color="C0",
        lw=2.0,
        label=r"$u/u_0$"
    )

    ax.plot(
        xi,
        profiles["s"],
        color="C1",
        lw=2.0,
        label=r"$s/s_0$"
    )

    ax.plot(
        xi,
        profiles["c"],
        color="C2",
        lw=2.0,
        label=r"$c/c_0$"
    )

    ax.plot(
        xi,
        profiles["w"],
        color="C3",
        lw=2.0,
        label=r"$\bar W/w_0$"
    )

    ax.set_title(title)
    ax.set_xlabel(r"$\xi=x-x_c$")
    ax.set_ylim(0.0, 1.5)
    ax.grid(True, alpha=0.25)


def plot_cover_profiles(data_root: Path, fig_dir: Path, supp_dir: Path, radius: float):
    profile_folders = [
        ("open", data_root / "profile_open_end"),
        ("closed", data_root / "profile_closed_end"),
        ("reopened", data_root / "profile_reopened_end"),
    ]

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(15.5, 4.0),
        sharey=True
    )

    for ax, (name, folder) in zip(axes, profile_folders):
        profiles = load_profiles_in_comoving_frame(folder)
        plot_case_profile(ax, profiles, name, radius)

    axes[0].set_ylabel("normalized concentration")

    handles, labels = axes[-1].get_legend_handles_labels()
    axes[-1].legend(handles, labels, loc="upper right", fontsize=9)

    fig.tight_layout()

    # Main text uses only figure5_cover_protocol_clean.
    # This profile figure is saved with the Supplementary filename used in supplementary_material.tex,
    # and also with the legacy filename for compatibility.
    out = supp_dir / "figS1_cover_protocol_profiles.png"
    legacy = fig_dir / "figure5_cover_protocol_profiles.png"
    save_figure(fig, out, dpi=DPI, aliases=[legacy])
    plt.close(fig)

    print(f"saved: {out}")


# ============================================================
# Main
# ============================================================

def main():
    """
    Expected Makefile call:
        python3 scripts/plot_figure6_cover_protocol.py figure5_data figure 1.0

    Arguments:
        argv[1] = data root directory, e.g. figure5_data
        argv[2] = figure directory, e.g. figure
        argv[3] = droplet radius, e.g. 1.0
    """
    data_root = Path("figure5_data")
    fig_dir = Path("figure")
    radius = DEFAULT_RADIUS

    supp_dir = Path("supplement_figures")

    if len(sys.argv) >= 2:
        data_root = Path(sys.argv[1])

    if len(sys.argv) >= 3:
        fig_dir = Path(sys.argv[2])

    if len(sys.argv) >= 4:
        # Backward compatible parsing:
        #   plot_figure6_cover_protocol.py data_root fig_dir radius
        # New optional form:
        #   plot_figure6_cover_protocol.py data_root fig_dir supp_dir radius
        try:
            radius = float(sys.argv[3])
        except ValueError:
            supp_dir = Path(sys.argv[3])
            if len(sys.argv) >= 5:
                radius = float(sys.argv[4])

    fig_dir.mkdir(parents=True, exist_ok=True)
    supp_dir.mkdir(parents=True, exist_ok=True)

    plot_cover_protocol_clean(data_root, fig_dir)
    plot_cover_profiles(data_root, fig_dir, supp_dir, radius)


if __name__ == "__main__":
    main()