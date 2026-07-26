#!/usr/bin/env python3
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DPI = 250
LEGACY_FIG_DIR = None
FONT_FAMILY = "Arial"
BASE_FONT_SIZE = 16
TICK_FONT_SIZE = 14
TITLE_FONT_SIZE = 18


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


def safe_read(path):
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path)


def ensure_dir(path: Path):
    path.mkdir(parents=True, exist_ok=True)


def savefig(fig, path: Path, aliases=None):
    targets = [Path(path)]
    if aliases:
        targets.extend(Path(a) for a in aliases)
    apply_figure_style(fig)
    for target in targets:
        target.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(target, dpi=DPI, bbox_inches="tight")
        fig.savefig(target.with_suffix(".eps"), format="eps", dpi=DPI, bbox_inches="tight")
        print(f"saved: {target}")
        print(f"saved: {target.with_suffix('.eps')}")
    plt.close(fig)


def plot_grid(df, fig_dir):
    g = df[df["suite"].str.contains("grid_Nx", na=False)].copy()
    if g.empty:
        return
    # Representative cases only for readability.
    reps = g[g["suite"].str.contains("representative", na=False)].copy()
    if reps.empty:
        reps = g

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for nx, sub in reps.groupby("Nx"):
        sub = sub.sort_values("s0")
        ax.plot(sub["s0"], sub["amp_theta"], marker="o", label=f"Nx={nx}")
    ax.set_xlabel(r"$s_0$")
    ax.set_ylabel(r"late amplitude $\Delta\theta$")
    ax.grid(True, alpha=0.3)
    ax.legend()
    savefig(fig, fig_dir / "figS3a_grid_refinement_amplitude.png", aliases=[LEGACY_FIG_DIR / "grid_refinement_amplitude.png"] if LEGACY_FIG_DIR is not None else None)

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for nx, sub in reps.groupby("Nx"):
        sub = sub.sort_values("s0")
        ax.plot(sub["s0"], sub["vrms_late"], marker="o", label=f"Nx={nx}")
    ax.set_xlabel(r"$s_0$")
    ax.set_ylabel(r"late activity $V_{\rm rms}$")
    ax.grid(True, alpha=0.3)
    ax.legend()
    savefig(fig, fig_dir / "figS3b_grid_refinement_vrms.png", aliases=[LEGACY_FIG_DIR / "grid_refinement_vrms.png"] if LEGACY_FIG_DIR is not None else None)


def plot_time(df, fig_dir):
    t = df[df["suite"].str.startswith("dt_", na=False)].copy()
    if t.empty:
        return
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for s0, sub in t.groupby("s0"):
        sub = sub.sort_values("dt")
        ax.plot(sub["dt"], sub["amp_theta"], marker="o", label=fr"$s_0={s0:.3f}$")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\Delta t$")
    ax.set_ylabel(r"late amplitude $\Delta\theta$")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    savefig(fig, fig_dir / "figS4a_timestep_amplitude.png", aliases=[LEGACY_FIG_DIR / "time_step_refinement_amplitude.png"] if LEGACY_FIG_DIR is not None else None)

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for s0, sub in t.groupby("s0"):
        sub = sub.sort_values("dt")
        ax.plot(sub["dt"], sub["vrms_late"], marker="o", label=fr"$s_0={s0:.3f}$")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\Delta t$")
    ax.set_ylabel(r"late activity $V_{\rm rms}$")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    savefig(fig, fig_dir / "figS4b_timestep_vrms.png", aliases=[LEGACY_FIG_DIR / "time_step_refinement_vrms.png"] if LEGACY_FIG_DIR is not None else None)


def plot_transition(df, fig_dir):
    tr = df[df["suite"].str.startswith("transition_", na=False)].copy()
    if tr.empty:
        return
    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    for suite, sub in tr.groupby("suite"):
        sub = sub.sort_values("s0")
        ax.plot(sub["s0"], sub["amp_theta"], marker="o", label=suite)
    ax.axhline(2*np.pi, color="black", lw=0.9, alpha=0.5)
    ax.set_xlabel(r"$s_0$")
    ax.set_ylabel(r"late amplitude $\Delta\theta$")
    ax.grid(True, alpha=0.3)
    ax.legend()
    savefig(fig, fig_dir / "figS5a_transition_refinement.png", aliases=[LEGACY_FIG_DIR / "transition_refinement_amplitude.png"] if LEGACY_FIG_DIR is not None else None)


def plot_long(df, fig_dir):
    lg = df[df["suite"].str.startswith("longtime_", na=False)].copy()
    if lg.empty:
        return
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for suite, sub in lg.groupby("suite"):
        sub = sub.sort_values("s0")
        ax.plot(sub["s0"], sub["amp_theta"], marker="o", label=suite)
    ax.set_xlabel(r"$s_0$")
    ax.set_ylabel(r"late amplitude $\Delta\theta$")
    ax.grid(True, alpha=0.3)
    ax.legend()
    savefig(fig, fig_dir / "figS5b_longtime_check.png", aliases=[LEGACY_FIG_DIR / "longtime_check_amplitude.png"] if LEGACY_FIG_DIR is not None else None)


def plot_perturb(df, fig_dir):
    if df.empty:
        return
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for eps, sub in df.groupby("eps"):
        sub = sub.sort_values("s0")
        ax.plot(sub["s0"], sub["late_delta_amp"], marker="o", label=fr"$\epsilon={eps:g}$")
    ax.set_xlabel(r"$s_0$")
    ax.set_ylabel(r"late amplitude $A_{\delta x}$")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    savefig(fig, fig_dir / "figS6a_perturbation_epsilon_amplitude.png", aliases=[LEGACY_FIG_DIR / "perturbation_epsilon_amplitude.png"] if LEGACY_FIG_DIR is not None else None)

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for eps, sub in df.groupby("eps"):
        sub = sub.sort_values("s0")
        ax.plot(sub["s0"], sub["late_vrms"], marker="o", label=fr"$\epsilon={eps:g}$")
    ax.set_xlabel(r"$s_0$")
    ax.set_ylabel(r"late activity $V_{\rm rms}$")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    savefig(fig, fig_dir / "figS6b_perturbation_epsilon_vrms.png", aliases=[LEGACY_FIG_DIR / "perturbation_epsilon_vrms.png"] if LEGACY_FIG_DIR is not None else None)


def plot_cover(df, fig_dir):
    if df.empty:
        return
    df = df.sort_values("kw_open")
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.plot(df["kw_open"], df["open_vrms"], marker="o", label="open")
    ax.plot(df["kw_open"], df["closed_vrms"], marker="o", label="closed")
    ax.plot(df["kw_open"], df["reopen_vrms"], marker="o", label="reopened")
    ax.set_xscale("log")
    ax.set_xlabel(r"open-cover leakage $k_w^+$")
    ax.set_ylabel(r"final-window activity $V_{\rm rms}$")
    ax.grid(True, alpha=0.3)
    ax.legend()
    savefig(fig, fig_dir / "figS7a_cover_kw_vrms.png", aliases=[LEGACY_FIG_DIR / "cover_kw_sensitivity_vrms.png"] if LEGACY_FIG_DIR is not None else None)

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.plot(df["kw_open"], df["open_wbar_end"], marker="o", label="open")
    ax.plot(df["kw_open"], df["closed_wbar_end"], marker="o", label="closed")
    ax.plot(df["kw_open"], df["reopen_wbar_end"], marker="o", label="reopened")
    ax.set_xscale("log")
    ax.set_xlabel(r"open-cover leakage $k_w^+$")
    ax.set_ylabel(r"end value $\bar W/w_0$")
    ax.grid(True, alpha=0.3)
    ax.legend()
    savefig(fig, fig_dir / "figS7b_cover_kw_wbar.png", aliases=[LEGACY_FIG_DIR / "cover_kw_sensitivity_wbar.png"] if LEGACY_FIG_DIR is not None else None)


def write_latex_snippets(out_dir, open_df, pert_df, cover_df):
    path = out_dir / "latex_snippets.tex"
    lines = []
    lines.append("% Auto-generated by scripts/plot_supplement_robustness_checks.py")
    lines.append("% Copy selected snippets into the manuscript after checking the values.")
    lines.append("")

    # Representative grid table.
    g = open_df[open_df["suite"].str.contains("grid_Nx", na=False) & open_df["suite"].str.contains("representative", na=False)].copy()
    if not g.empty:
        lines.append(r"\begin{table}[t]")
        lines.append(r"\centering")
        lines.append(r"\caption{Grid-refinement check for representative single-droplet states.}")
        lines.append(r"\label{tab:grid-refinement-check}")
        lines.append(r"\begin{tabular}{cccccc}")
        lines.append(r"\toprule")
        lines.append(r"$s_0$ & $N_x$ & $\Delta x$ & state & $\Delta\theta$ & $V_{\rm rms}$\\")
        lines.append(r"\midrule")
        for _, row in g.sort_values(["s0", "Nx"]).iterrows():
            amp = row["amp_theta"]
            vrms = row["vrms_late"]
            lines.append(f"{row['s0']:.3f} & {int(row['Nx'])} & {row['dx']:.5g} & {row['state']} & {amp:.4g} & {vrms:.4g}\\")
        lines.append(r"\bottomrule")
        lines.append(r"\end{tabular}")
        lines.append(r"\end{table}")
        lines.append("")

    # Perturbation epsilon summary.
    if not pert_df.empty:
        lines.append("% Perturbation-amplitude robustness summary")
        pivot = pert_df.pivot_table(index="s0", columns="eps", values="late_vrms", aggfunc="mean")
        lines.append("% late V_rms values:")
        for s0, row in pivot.iterrows():
            vals = ", ".join([f"eps={eps:g}: {val:.4g}" for eps, val in row.dropna().items()])
            lines.append(f"% s0={s0:.3f}: {vals}")
        lines.append("")

    if not cover_df.empty:
        lines.append("% Cover leakage sensitivity summary")
        for _, row in cover_df.sort_values("kw_open").iterrows():
            lines.append(f"% kw_open={row['kw_open']:.4g}: open V={row['open_vrms']:.4g}, closed V={row['closed_vrms']:.4g}, reopened V={row['reopen_vrms']:.4g}; closed W={row['closed_wbar_end']:.4g}")

    path.write_text("\n".join(lines))
    print(f"saved: {path}")


def main():
    out_dir = Path(sys.argv[1]) if len(sys.argv) >= 2 else Path("defense_data")
    global LEGACY_FIG_DIR
    fig_dir = Path("supplement_figures")
    LEGACY_FIG_DIR = out_dir / "figures"
    ensure_dir(fig_dir)
    ensure_dir(LEGACY_FIG_DIR)

    open_df = safe_read(out_dir / "open_case_metrics.csv")
    pert_df = safe_read(out_dir / "perturbation_metrics.csv")
    cover_df = safe_read(out_dir / "cover_sensitivity.csv")

    if not open_df.empty:
        open_df.to_csv(out_dir / "open_case_metrics_sorted.csv", index=False)
        plot_grid(open_df, fig_dir)
        plot_time(open_df, fig_dir)
        plot_transition(open_df, fig_dir)
        plot_long(open_df, fig_dir)

    if not pert_df.empty:
        pert_df.to_csv(out_dir / "perturbation_metrics_sorted.csv", index=False)
        plot_perturb(pert_df, fig_dir)

    if not cover_df.empty:
        cover_df.to_csv(out_dir / "cover_sensitivity_sorted.csv", index=False)
        plot_cover(cover_df, fig_dir)

    write_latex_snippets(out_dir, open_df, pert_df, cover_df)

    print("Done.")


if __name__ == "__main__":
    main()
