#!/usr/bin/env python3
"""
LDPC Code Optimization Process — Dynamic Visualization
=======================================================
Reads H_SV_RES.mat from the parent 'example' folder and animates the
LDPC optimization trajectory across all 200 optimization rounds.

Layout (top → bottom):
  • Shift-value matrix  (H, 13 × 35)  — full width
  • Decoding threshold  |  Normalized complexity   (time-series)
  • FER curves          |  Average-iteration curves

The animation runs for exactly TOTAL_ROUNDS frames (one per optimization
round). Plots are refreshed only when the round counter reaches a saved
checkpoint; the round counter in the title updates every frame.

Run:
    python ldpc_visualization.py

Configuration constants are near the top of this file.
"""

import os
import sys
import bisect
import warnings

import numpy as np
import scipy.io
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.animation import FuncAnimation
from matplotlib.colors import ListedColormap
import matplotlib.ticker as ticker
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────────────────────────────────────
# USER CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────
UPDATE_INTERVAL_MS: int = 100    # ms between animation frames (change freely)
TOTAL_ROUNDS: int = 200          # total number of optimization rounds

# ─────────────────────────────────────────────────────────────────────────────
# PATHS  (H_SV_RES.mat lives one directory above this script)
# ─────────────────────────────────────────────────────────────────────────────
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
MAT_FILE = os.path.join(_THIS_DIR, "..", "H_SV_RES.mat")

# ─────────────────────────────────────────────────────────────────────────────
# PROFESSIONAL COLOR PALETTE  (all sub-plots share the same family)
# ─────────────────────────────────────────────────────────────────────────────
PAL = dict(
    fig_bg      = "#F0F4F8",      # outer figure background
    ax_bg       = "#FFFFFF",      # individual axes background
    grid        = "#DDE3EC",      # major grid lines
    text_dark   = "#1E293B",      # titles, bold labels
    text_mid    = "#475569",      # axis labels, tick labels
    # curves
    initial     = "#E53E3E",      # initial-state FER / AvgIT  (red)
    current     = "#2B6CB0",      # current-state curve         (royal blue)
    history     = "#90CDF4",      # faded history lines         (light blue)
    hist_alpha  = 0.25,           # alpha for history lines
    # scalar plots
    threshold   = "#276749",      # emerald green
    threshold_pt= "#38A169",
    complexity  = "#553C9A",      # violet
    complexity_pt= "#7C3AED",
    # H-matrix colormap  (sequential, perceptually uniform)
    hmat_cmap   = "YlOrBr",       # yellow → orange → brown  (non-zero shifts)
    hmat_zero   = "#F0F4F8",      # background colour for zero blocks (-1)
    hmat_grid   = "#C8D0DC",      # thin cell-separator grid
)

# ─────────────────────────────────────────────────────────────────────────────
# FONT / STYLE SETUP  (enlarged for readability)
# ─────────────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family"          : "DejaVu Sans",
    "font.size"            : 11,
    "axes.titlesize"       : 13,
    "axes.titleweight"     : "bold",
    "axes.labelsize"       : 11,
    "axes.labelcolor"      : PAL["text_mid"],
    "axes.edgecolor"       : "#B0BAC8",
    "axes.linewidth"       : 1.5,
    "xtick.color"          : PAL["text_mid"],
    "ytick.color"          : PAL["text_mid"],
    "xtick.labelsize"      : 10,
    "ytick.labelsize"      : 10,
    "xtick.major.width"    : 1.5,
    "ytick.major.width"    : 1.5,
    "xtick.major.size"     : 5,
    "ytick.major.size"     : 5,
    "xtick.minor.width"    : 1.0,
    "ytick.minor.width"    : 1.0,
    "grid.color"           : PAL["grid"],
    "grid.linewidth"       : 0.8,
    "legend.fontsize"      : 10,
    "legend.framealpha"    : 0.9,
    "figure.facecolor"     : PAL["fig_bg"],
    "axes.facecolor"       : PAL["ax_bg"],
    "savefig.facecolor"    : PAL["fig_bg"],
})

# ─────────────────────────────────────────────────────────────────────────────
# LOAD DATA
# ─────────────────────────────────────────────────────────────────────────────
def load_data(path: str):
    mat = scipy.io.loadmat(path, squeeze_me=False)
    res = mat["H_SV_RES"]              # (N_states, 6) object array

    N = res.shape[0]
    H_mats      = [res[i, 0].astype(np.int16)        for i in range(N)]
    iter_nums   = [int(res[i, 1].flat[0])             for i in range(N)]
    thresholds  = [float(res[i, 2].flat[0])           for i in range(N)]
    complexities= [float(res[i, 3].flat[0])           for i in range(N)]
    fer_curves  = [res[i, 4].flatten().astype(float)  for i in range(N)]
    avgit_curves= [res[i, 5].flatten().astype(float)  for i in range(N)]

    ebn0 = np.array([2.8, 3.2, 3.6, 4.0, 4.4, 4.8])
    return N, H_mats, iter_nums, thresholds, complexities, fer_curves, avgit_curves, ebn0

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE / AXES SETUP
# ─────────────────────────────────────────────────────────────────────────────
def build_figure():
    # Tall figure: H matrix spans full width on top; 2×2 grid of plots below.
    fig = plt.figure(figsize=(18, 14))
    fig.patch.set_facecolor(PAL["fig_bg"])

    # Outer GridSpec: 2 rows — row 0 = H matrix, row 1 = 2×2 scalar/curve plots.
    # height_ratios chosen so the H matrix (13×35 cells, aspect="equal") gets
    # enough vertical room at full figure width.
    outer = gridspec.GridSpec(
        2, 1,
        figure=fig,
        left=0.06, right=0.97,
        top=0.93,  bottom=0.07,
        hspace=0.38,
        height_ratios=[1.2, 1],
    )

    ax_H = fig.add_subplot(outer[0])

    inner = gridspec.GridSpecFromSubplotSpec(
        2, 2,
        subplot_spec=outer[1],
        hspace=0.48, wspace=0.38,
    )
    ax_thr  = fig.add_subplot(inner[0, 0])
    ax_cmp  = fig.add_subplot(inner[0, 1])
    ax_fer  = fig.add_subplot(inner[1, 0])
    ax_avg  = fig.add_subplot(inner[1, 1])

    for ax in (ax_H, ax_thr, ax_cmp, ax_fer, ax_avg):
        ax.set_facecolor(PAL["ax_bg"])
        ax.tick_params(which="major", direction="in", length=5, width=1.5)
        ax.tick_params(which="minor", direction="in", length=3, width=1.0)

    return fig, ax_H, ax_thr, ax_cmp, ax_fer, ax_avg

# ─────────────────────────────────────────────────────────────────────────────
# H MATRIX HELPER
# ─────────────────────────────────────────────────────────────────────────────
def draw_H_matrix(ax, H: np.ndarray, Z: int = 64):
    """Render the 13 × 35 shift-value matrix.

    -1  → background colour (zero block)
    ≥ 0 → single accent colour (cyclic-shift block present)
    """
    ax.cla()

    rows, cols = H.shape          # 13, 35
    # Binary: 0 = zero block (-1), 1 = non-zero block (≥ 0)
    binary = (H >= 0).astype(float)

    cmap_two = ListedColormap([PAL["hmat_zero"], "#4A90D9"])

    img = ax.imshow(
        binary,
        cmap=cmap_two,
        vmin=0, vmax=1,
        aspect="equal",
        interpolation="nearest",
        origin="upper",
    )

    # thin cell grid
    for x in np.arange(-0.5, cols, 1):
        ax.axvline(x, color=PAL["hmat_grid"], linewidth=0.45, zorder=3)
    for y in np.arange(-0.5, rows, 1):
        ax.axhline(y, color=PAL["hmat_grid"], linewidth=0.45, zorder=3)

    ax.set_xlim(-0.5, cols - 0.5)
    ax.set_ylim(rows - 0.5, -0.5)
    ax.set_xticks(np.arange(0, cols, 5))
    ax.set_yticks(np.arange(rows))
    ax.set_yticklabels([str(i) for i in range(rows)], fontsize=9)
    ax.set_xticklabels([str(i) for i in range(0, cols, 5)], fontsize=9)
    ax.set_xlabel("Column index", labelpad=4)
    ax.set_ylabel("Row index", labelpad=4)
    ax.set_title("Shift-Value Matrix  H  (13 × 35)", color=PAL["text_dark"], pad=6)

    return img

# ─────────────────────────────────────────────────────────────────────────────
# SCALAR TIME-SERIES HELPERS
# ─────────────────────────────────────────────────────────────────────────────
def setup_scalar_ax(ax, title: str, ylabel: str, color: str):
    ax.set_title(title, color=PAL["text_dark"], pad=5)
    ax.set_xlabel("Optimization round", labelpad=4)
    ax.set_ylabel(ylabel, labelpad=4, color=PAL["text_mid"])
    ax.set_xlim(0, TOTAL_ROUNDS)
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.3g"))
    ax.grid(True, axis="both")
    ax.set_axisbelow(True)


def update_scalar(ax, xs, ys, color: str, marker_color: str, ylabel: str, title: str):
    ax.cla()
    setup_scalar_ax(ax, title, ylabel, color)
    if len(xs) == 0:
        return
    ax.plot(xs, ys,
            color=color, linewidth=1.4, alpha=0.85,
            marker="o", markersize=4, markerfacecolor=marker_color,
            markeredgewidth=0.6, markeredgecolor="white", zorder=4)
    # fill under
    ax.fill_between(xs, ys, min(ys) - 0.01 * abs(min(ys) + 1e-9),
                    color=color, alpha=0.10)
    # highlight last point
    if len(xs) >= 1:
        ax.scatter([xs[-1]], [ys[-1]], color=marker_color,
                   s=40, zorder=5, edgecolors="white", linewidths=0.8)
    ax.margins(x=0.04, y=0.12)
    ax.grid(True)

# ─────────────────────────────────────────────────────────────────────────────
# FER / AVGIT CURVE HELPERS
# ─────────────────────────────────────────────────────────────────────────────
def setup_curve_ax(ax, title: str, ylabel: str, log_scale: bool = True):
    ax.set_title(title, color=PAL["text_dark"], pad=5)
    ax.set_xlabel("$E_b / N_0$  (dB)", labelpad=4)
    ax.set_ylabel(ylabel, labelpad=4)
    if log_scale:
        ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.5)
    ax.set_axisbelow(True)
    ax.set_xlim(2.6, 5.0)


def update_curves(ax, ebn0, initial_curve, history_curves, current_curve,
                  title: str, ylabel: str, log_scale: bool = True):
    ax.cla()
    setup_curve_ax(ax, title, ylabel, log_scale=log_scale)

    def _plot(curve):
        if log_scale:
            return np.clip(curve, 1e-6, 1)
        return curve

    def _draw(curve, **kwargs):
        if log_scale:
            ax.semilogy(ebn0, _plot(curve), **kwargs)
        else:
            ax.plot(ebn0, curve, **kwargs)

    # ── history curves (all previous optimization states, faded) ──
    for curve in history_curves:
        _draw(curve,
              color=PAL["history"], linewidth=0.9,
              alpha=PAL["hist_alpha"], zorder=2)

    # ── current optimization curve ──
    if current_curve is not None:
        _draw(current_curve,
              color=PAL["current"], linewidth=2.0,
              marker="D", markersize=4.5,
              markerfacecolor=PAL["current"],
              markeredgecolor="white", markeredgewidth=0.7,
              zorder=4, label="Current state")

    # ── initial curve (always on top, in red) ──
    if initial_curve is not None:
        _draw(initial_curve,
              color=PAL["initial"], linewidth=2.0,
              linestyle="--",
              marker="s", markersize=4.5,
              markerfacecolor=PAL["initial"],
              markeredgecolor="white", markeredgewidth=0.7,
              zorder=5, label="Initial state")

    # legend
    legend_handles = []
    if initial_curve is not None:
        legend_handles.append(
            Line2D([0], [0], color=PAL["initial"], linewidth=1.8,
                   linestyle="--", marker="s", markersize=4, label="Initial state"))
    if current_curve is not None:
        legend_handles.append(
            Line2D([0], [0], color=PAL["current"], linewidth=1.8,
                   marker="D", markersize=4, label="Current state"))
    if history_curves:
        legend_handles.append(
            Line2D([0], [0], color=PAL["history"], linewidth=1.5,
                   alpha=0.6, label=f"Previous ({len(history_curves)})"))
    if legend_handles:
        ax.legend(handles=legend_handles, loc="upper right",
                  framealpha=0.9, edgecolor="#B0BAC8", fontsize=7.5)

# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────
def main():
    if not os.path.isfile(MAT_FILE):
        sys.exit(f"[ERROR] Cannot find MAT file:\n  {MAT_FILE}")

    (N, H_mats, iter_nums, thresholds,
     complexities, fer_curves, avgit_curves, ebn0) = load_data(MAT_FILE)

    # iter_nums must be sorted (they represent checkpoint round indices)
    sorted_rounds = sorted(iter_nums)

    fig, ax_H, ax_thr, ax_cmp, ax_fer, ax_avg = build_figure()

    # ── static super-title with round counter ──
    title_text = fig.suptitle(
        f"LDPC Code Optimization  —  Round 1 / {TOTAL_ROUNDS}",
        fontsize=15, fontweight="bold",
        color=PAL["text_dark"], y=0.975,
    )

    # ── state tracker: remember which checkpoint was last drawn ──
    state = dict(prev_state_idx=-2)   # -2 = never drawn

    # ─── animation update function ───────────────────────────────────────────
    def update(frame: int):
        # frame ∈ [0, TOTAL_ROUNDS) → display round = frame + 1
        round_num = frame + 1

        # Always refresh the round counter in the title
        title_text.set_text(
            f"LDPC Code Optimization  —  Round {round_num} / {TOTAL_ROUNDS}"
        )

        # Find the most recent checkpoint at or before the current round
        # bisect_right gives insertion point; subtract 1 for the last ≤ round_num
        pos = bisect.bisect_right(sorted_rounds, round_num) - 1
        if pos < 0:
            # No checkpoint reached yet — nothing to draw except the title
            fig.canvas.draw_idle()
            return

        state_idx = iter_nums.index(sorted_rounds[pos])

        if state_idx == state["prev_state_idx"]:
            # Checkpoint unchanged — only the title needs redrawing
            fig.canvas.draw_idle()
            return

        state["prev_state_idx"] = state_idx

        # ── H matrix ──
        draw_H_matrix(ax_H, H_mats[state_idx])

        # ── scalar time-series (all checkpoints up to and including state_idx) ──
        xs_shown  = [iter_nums[i]    for i in range(state_idx + 1)]
        thr_shown = [thresholds[i]   for i in range(state_idx + 1)]
        cmp_shown = [complexities[i] for i in range(state_idx + 1)]

        update_scalar(ax_thr, xs_shown, thr_shown,
                      PAL["threshold"], PAL["threshold_pt"],
                      "Threshold (dB)", "Decoding Threshold")
        update_scalar(ax_cmp, xs_shown, cmp_shown,
                      PAL["complexity"], PAL["complexity_pt"],
                      "Complexity index", "Normalised Complexity")

        # ── FER / AvgIT curves ──
        initial_fer = fer_curves[0]
        initial_avg = avgit_curves[0]

        if state_idx == 0:
            update_curves(ax_fer, ebn0, initial_fer, [], None,
                          "Frame Error Rate (FER)", "FER")
            update_curves(ax_avg, ebn0, initial_avg, [], None,
                          "Average Iterations (AvgIT)", "AvgIT", log_scale=False)
        else:
            history_fer = fer_curves[1:state_idx]
            history_avg = avgit_curves[1:state_idx]
            current_fer = fer_curves[state_idx]
            current_avg = avgit_curves[state_idx]

            update_curves(ax_fer, ebn0, initial_fer, history_fer, current_fer,
                          "Frame Error Rate (FER)", "FER")
            update_curves(ax_avg, ebn0, initial_avg, history_avg, current_avg,
                          "Average Iterations (AvgIT)", "AvgIT", log_scale=False)

        fig.canvas.draw_idle()

    anim = FuncAnimation(
        fig,
        update,
        frames=TOTAL_ROUNDS,        # exactly 200 frames = 200 optimization rounds
        interval=UPDATE_INTERVAL_MS,
        repeat=False,
    )

    plt.show()


if __name__ == "__main__":
    main()
