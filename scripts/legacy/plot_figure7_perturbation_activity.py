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

# Use the final part of each perturbation simulation to measure activity.
# If t_free = 80, this uses the last 20 time units.
LATE_WINDOW = 20.0

# For representative perturbation plots.
EXAMPLE_TMAX = 20.0

# Window for smoothed velocity RMS in example plots.
VRMS_WINDOW = 0.5

# Envelope for representative displacement response plots.
ENVELOPE_WINDOW_EXAMPLE = 0.25

# Plot style for raw oscillatory signal.
RAW_ALPHA = 0.18
RAW_LINEWIDTH = 0.35
ENVELOPE_LINEWIDTH = 1.8

# Target cases shown in the representative response figure.
# The script chooses the nearest available s0 values.
EXAMPLE_TARGETS = [0.095, 0.100, 0.105]


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
                "folder": row["folder"],
                "eps_shift": float(row["eps_shift"]),
                "t_relax": float(row["t_relax"]),
                "t_free": float(row["t_free"]),
                "dt": float(row["dt"]),
                "save_every": int(float(row["save_every"])),
            })

    return cases


def read_csv_named(path: Path):
    data = np.genfromtxt(path, delimiter=",", names=True)

    if data.ndim == 0:
        data = np.array([data])

    return data


# ============================================================
# Diagnostics
# ============================================================

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


def sliding_max_abs(t, y, window):
    """
    Sliding envelope:
        E(t_i) = max |y| over a window centered at t_i.
    """
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)

    n = len(t)
    out = np.zeros(n)

    half = 0.5 * window
    j0 = 0
    j1 = 0

    ay = np.abs(y)

    for i in range(n):
        left = t[i] - half
        right = t[i] + half

        while j0 < n and t[j0] < left:
            j0 += 1

        while j1 < n and t[j1] <= right:
            j1 += 1

        seg = ay[j0:j1]

        if len(seg) >= 1:
            out[i] = np.max(seg)
        else:
            out[i] = 0.0

    return out


def moving_average(y, k):
    """
    Simple moving average.
    """
    y = np.asarray(y, dtype=float)

    if k <= 1:
        return y.copy()

    k = int(k)
    pad = k // 2

    yy = np.pad(y, (pad, pad), mode="edge")
    kernel = np.ones(k) / float(k)

    out = np.convolve(yy, kernel, mode="valid")

    if len(out) > len(y):
        out = out[:len(y)]
    elif len(out) < len(y):
        out = np.pad(out, (0, len(y) - len(out)), mode="edge")

    return out


def late_window_mask(t, window=LATE_WINDOW):
    """
    Mask for the final time window.
    """
    t = np.asarray(t, dtype=float)

    tmax = np.max(t)
    return t >= (tmax - window)


def displacement_amplitude(y):
    """
    max(y)-min(y).
    This removes constant offsets from the translational neutral mode.
    """
    y = np.asarray(y, dtype=float)

    if len(y) == 0:
        return np.nan

    return float(np.max(y) - np.min(y))


def compute_metrics_for_case(data_root: Path, cas):
    folder = cas["folder"]
    path = data_root / folder / "perturbation.csv"

    data = read_csv_named(path)

    t = np.asarray(data["t"], dtype=float)
    delta = np.asarray(data["delta_x"], dtype=float)

    if "velocity" in data.dtype.names:
        velocity = np.asarray(data["velocity"], dtype=float)
    else:
        xc = np.asarray(data["xc"], dtype=float)
        velocity = np.gradient(xc, t)

    mask = late_window_mask(t, LATE_WINDOW)

    late_delta = delta[mask]
    late_velocity = velocity[mask]

    late_delta_amp = displacement_amplitude(late_delta)

    if len(late_velocity) >= 2:
        late_vrms = float(np.sqrt(np.mean(late_velocity * late_velocity)))
    else:
        late_vrms = np.nan

    max_abs_delta = float(np.max(np.abs(delta)))
    max_abs_velocity = float(np.max(np.abs(velocity)))

    return {
        "s0": cas["s0"],
        "folder": folder,
        "late_delta_amp": late_delta_amp,
        "late_vrms": late_vrms,
        "max_abs_delta": max_abs_delta,
        "max_abs_velocity": max_abs_velocity,
    }


def estimate_activity_onset_s0(s0, activity):
    """
    Estimate transition point using a simple 5% activity threshold.
    This is not an eigenvalue calculation; it is a visual activity-onset guide.
    """
    s0 = np.asarray(s0, dtype=float)
    activity = np.asarray(activity, dtype=float)

    order = np.argsort(s0)
    s0 = s0[order]
    activity = activity[order]

    finite = np.isfinite(activity)

    if np.count_nonzero(finite) < 3:
        return None

    max_act = np.max(activity[finite])

    if max_act <= 0:
        return None

    threshold = 0.05 * max_act

    for i in range(1, len(s0)):
        if not np.isfinite(activity[i - 1]) or not np.isfinite(activity[i]):
            continue

        if activity[i - 1] < threshold <= activity[i]:
            x0 = s0[i - 1]
            x1 = s0[i]
            y0 = activity[i - 1]
            y1 = activity[i]

            if y1 == y0:
                return x1

            return x0 + (threshold - y0) * (x1 - x0) / (y1 - y0)

    return None


# ============================================================
# Plot 1: Activity onset near transition
# ============================================================

def plot_activity_onset(fig_dir: Path, metrics):
    s0 = np.array([m["s0"] for m in metrics], dtype=float)
    amp = np.array([m["late_delta_amp"] for m in metrics], dtype=float)
    vrms = np.array([m["late_vrms"] for m in metrics], dtype=float)

    order = np.argsort(s0)
    s0 = s0[order]
    amp = amp[order]
    vrms = vrms[order]

    onset_amp = estimate_activity_onset_s0(s0, amp)
    onset_vrms = estimate_activity_onset_s0(s0, vrms)

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(7.4, 6.2),
        sharex=True
    )

    # --------------------------------------------------------
    # Late displacement amplitude
    # --------------------------------------------------------
    axes[0].plot(
        s0,
        amp,
        "o-",
        color="black",
        lw=1.3,
        markersize=4
    )

    axes[0].set_ylabel(r"late amplitude $A_{\delta x}$")
    axes[0].grid(True, alpha=0.25)

    if onset_amp is not None:
        axes[0].axvline(onset_amp, color="black", linestyle=":", lw=1.0)

    # --------------------------------------------------------
    # Late velocity RMS
    # --------------------------------------------------------
    axes[1].plot(
        s0,
        vrms,
        "o-",
        color="black",
        lw=1.3,
        markersize=4
    )

    axes[1].set_ylabel(r"late activity $V_{\rm rms}$")
    axes[1].set_xlabel(r"$s_0$")
    axes[1].grid(True, alpha=0.25)

    if onset_vrms is not None:
        axes[1].axvline(onset_vrms, color="black", linestyle=":", lw=1.0)

    fig.tight_layout()

    # New recommended filename.
    out = fig_dir / "figure6_activity_onset.png"
    save_figure(fig, out, dpi=DPI)
    print(f"saved: {out}")

    # Compatibility with old workflow / Makefile expectations.
    out_old = fig_dir / "figure6_growth_rate.png"
    save_figure(fig, out_old, dpi=DPI)
    print(f"saved: {out_old}")

    plt.close(fig)


# ============================================================
# Plot 2: Representative perturbation responses
# ============================================================

def choose_examples(metrics):
    """
    Choose examples closest to EXAMPLE_TARGETS.
    """
    examples = []
    used = set()

    for target in EXAMPLE_TARGETS:
        best = min(
            metrics,
            key=lambda m: abs(m["s0"] - target)
        )

        if best["folder"] not in used:
            examples.append(best)
            used.add(best["folder"])

    # If duplicates reduced the count, add low/mid/high activity cases.
    if len(examples) < 3:
        sorted_by_s0 = sorted(metrics, key=lambda m: m["s0"])
        fallback = [
            sorted_by_s0[0],
            sorted_by_s0[len(sorted_by_s0) // 2],
            sorted_by_s0[-1],
        ]

        for m in fallback:
            if m["folder"] not in used:
                examples.append(m)
                used.add(m["folder"])

            if len(examples) >= 3:
                break

    return examples


def plot_perturbation_examples(data_root: Path, fig_dir: Path, supp_dir: Path, metrics):
    examples = choose_examples(metrics)

    if len(examples) == 0:
        return

    n = len(examples)

    fig, axes = plt.subplots(
        2,
        n,
        figsize=(5.0 * n, 5.8),
        sharex=False
    )

    if n == 1:
        axes = np.array([[axes[0]], [axes[1]]])

    for j, m in enumerate(examples):
        data = read_csv_named(data_root / m["folder"] / "perturbation.csv")

        t = np.asarray(data["t"], dtype=float)
        delta = np.asarray(data["delta_x"], dtype=float)

        if "velocity" in data.dtype.names:
            velocity = np.asarray(data["velocity"], dtype=float)
        else:
            xc = np.asarray(data["xc"], dtype=float)
            velocity = np.gradient(xc, t)

        # Show early-to-intermediate response only.
        mask = t <= min(EXAMPLE_TMAX, np.max(t))

        t_show = t[mask]
        delta_show = delta[mask]
        velocity_show = velocity[mask]

        # Envelope of displacement perturbation.
        delta_env = sliding_max_abs(
            t_show,
            delta_show,
            ENVELOPE_WINDOW_EXAMPLE
        )
        delta_env = moving_average(delta_env, 5)

        # Velocity RMS.
        vrms_show = sliding_rms(
            t_show,
            velocity_show,
            VRMS_WINDOW
        )
        vrms_show = moving_average(vrms_show, 5)

        # ----------------------------------------------------
        # Top row: displacement perturbation with envelope
        # ----------------------------------------------------
        ax = axes[0, j]

        ax.axhline(0.0, color="black", lw=0.8, alpha=0.45)

        # Raw oscillatory response, shown faintly.
        ax.plot(
            t_show,
            delta_show,
            color="0.45",
            lw=RAW_LINEWIDTH,
            alpha=RAW_ALPHA,
            rasterized=True
        )

        # Envelope band.
        ax.fill_between(
            t_show,
            -delta_env,
            delta_env,
            color="0.80",
            alpha=0.45,
            linewidth=0
        )

        # Envelope curves.
        ax.plot(
            t_show,
            delta_env,
            color="black",
            lw=ENVELOPE_LINEWIDTH,
            label=r"$E_{\delta x}(t)$"
        )

        ax.plot(
            t_show,
            -delta_env,
            color="black",
            lw=ENVELOPE_LINEWIDTH
        )

        ax.set_title(rf"$s_0={m['s0']:.3f}$")
        ax.set_ylabel(r"$\delta x(t)$")
        ax.grid(True, alpha=0.25)

        ymax = np.nanmax(delta_env)

        if not np.isfinite(ymax) or ymax <= 0:
            ymax = np.nanmax(np.abs(delta_show))

        if np.isfinite(ymax) and ymax > 0:
            ax.set_ylim(-1.15 * ymax, 1.15 * ymax)

        ax.legend(loc="upper right", fontsize=8)

        # ----------------------------------------------------
        # Bottom row: velocity activity
        # ----------------------------------------------------
        ax = axes[1, j]

        ax.plot(
            t_show,
            vrms_show,
            color="C3",
            lw=1.8
        )

        ax.set_ylabel(r"$V_{\rm rms}(t)$")
        ax.set_xlabel(r"$t$")
        ax.grid(True, alpha=0.25)

    fig.tight_layout()

    # Main text uses figure6_activity_onset. Representative responses are Supplementary Fig. S2.
    out = supp_dir / "figS2_perturbation_examples.png"
    legacy = fig_dir / "figure6_perturbation_examples.png"
    save_figure(fig, out, dpi=DPI, aliases=[legacy])
    plt.close(fig)

    print(f"saved: {out}")


# ============================================================
# Save metrics
# ============================================================

def save_metrics(data_root: Path, metrics):
    out_csv = data_root / "activity_metrics.csv"

    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "s0",
                "folder",
                "late_delta_amp",
                "late_vrms",
                "max_abs_delta",
                "max_abs_velocity",
            ]
        )

        writer.writeheader()

        for m in metrics:
            writer.writerow(m)

    print(f"saved: {out_csv}")


# ============================================================
# Main
# ============================================================

def main():
    data_root = Path("figure6_data")
    fig_dir = Path("figure")

    supp_dir = Path("supplement_figures")

    if len(sys.argv) >= 2:
        data_root = Path(sys.argv[1])

    if len(sys.argv) >= 3:
        fig_dir = Path(sys.argv[2])

    if len(sys.argv) >= 4:
        supp_dir = Path(sys.argv[3])

    fig_dir.mkdir(parents=True, exist_ok=True)
    supp_dir.mkdir(parents=True, exist_ok=True)

    cases_path = data_root / "perturbation_cases.csv"
    cases = read_cases(cases_path)

    metrics = []

    for cas in cases:
        m = compute_metrics_for_case(data_root, cas)
        metrics.append(m)

    save_metrics(data_root, metrics)

    plot_activity_onset(fig_dir, metrics)
    plot_perturbation_examples(data_root, fig_dir, supp_dir, metrics)


if __name__ == "__main__":
    main()