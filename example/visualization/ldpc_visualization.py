#!/usr/bin/env python3
"""
LDPC Code Optimization Process — Dynamic Visualization
=======================================================
Reads H_SV_RES.mat from the sibling 'example' folder and animates the
LDPC optimization trajectory across up to 50 saved states.

Six sub-plots on one figure:
  • Shift-value matrix  (H, 13 × 35)
  • Decoding threshold  (dB, time-series)
  • Normalized complexity (time-series)
  • FER curves          (Eb/N0 sweep, accumulated)
  • Average-iteration curves (same convention as FER)

Run:
    python ldpc_visualization.py

Configuration constants are near the top of this file.
"""

import os
import sys
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
UPDATE_INTERVAL_MS: int = 1000   # ms between animation frames (change freely)

# ─────────────────────────────────────────────────────────────────────────────
# PATHS
# ─────────────────────────────────────────────────────────────────────────────
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
MAT_FILE = os.path.join(_THIS_DIR, "H_SV_RES.mat")

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
# FONT / STYLE SETUP
# ─────────────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family"        : "DejaVu Sans",
    "font.size"          : 9,
    "axes.titlesize"     : 10,
    "axes.titleweight"   : "bold",
    "axes.labelsize"     : 9,
    "axes.labelcolor"    : PAL["text_mid"],
    "axes.edgecolor"     : "#B0BAC8",
    "axes.linewidth"     : 0.8,
    "xtick.color"        : PAL["text_mid"],
    "ytick.color"        : PAL["text_mid"],
    "xtick.labelsize"    : 8,
    "ytick.labelsize"    : 8,
    "grid.color"         : PAL["grid"],
    "grid.linewidth"     : 0.6,
    "legend.fontsize"    : 8,
    "legend.framealpha"  : 0.9,
    "figure.facecolor"   : PAL["fig_bg"],
    "axes.facecolor"     : PAL["ax_bg"],
    "savefig.facecolor"  : PAL["fig_bg"],
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
    fig = plt.figure(figsize=(18, 9.5))
    fig.patch.set_facecolor(PAL["fig_bg"])

    # GridSpec: left column = H matrix (full height); right = 2 × 2 grid
    outer = gridspec.GridSpec(
        1, 2,
        figure=fig,
        left=0.05, right=0.97,
        top=0.92,  bottom=0.08,
        wspace=0.32,
        width_ratios=[1, 1.55],
    )

    ax_H = fig.add_subplot(outer[0, 0])

    inner = gridspec.GridSpecFromSubplotSpec(
        2, 2,
        subplot_spec=outer[0, 1],
        hspace=0.42, wspace=0.35,
    )
    ax_thr  = fig.add_subplot(inner[0, 0])
    ax_cmp  = fig.add_subplot(inner[0, 1])
    ax_fer  = fig.add_subplot(inner[1, 0])
    ax_avg  = fig.add_subplot(inner[1, 1])

    for ax in (ax_H, ax_thr, ax_cmp, ax_fer, ax_avg):
        ax.set_facecolor(PAL["ax_bg"])
        ax.tick_params(which="both", direction="in", length=3)

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
    ax.set_yticklabels([str(i) for i in range(rows)], fontsize=7)
    ax.set_xticklabels([str(i) for i in range(0, cols, 5)], fontsize=7)
    ax.set_xlabel("Column index", labelpad=3)
    ax.set_ylabel("Row index", labelpad=3)
    ax.set_title("Shift-Value Matrix  H  (13 × 35)", color=PAL["text_dark"], pad=5)

    return img

# ─────────────────────────────────────────────────────────────────────────────
# SCALAR TIME-SERIES HELPERS
# ─────────────────────────────────────────────────────────────────────────────
def setup_scalar_ax(ax, title: str, ylabel: str, color: str):
    ax.set_title(title, color=PAL["text_dark"], pad=4)
    ax.set_xlabel("Optimization round", labelpad=3)
    ax.set_ylabel(ylabel, labelpad=3, color=PAL["text_mid"])
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
    ax.set_title(title, color=PAL["text_dark"], pad=4)
    ax.set_xlabel("$E_b / N_0$  (dB)", labelpad=3)
    ax.set_ylabel(ylabel, labelpad=3)
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

    fig, ax_H, ax_thr, ax_cmp, ax_fer, ax_avg = build_figure()

    # ── static super-title with round counter ──
    title_text = fig.suptitle(
        "LDPC Code Optimization  —  Round 0",
        fontsize=13, fontweight="bold",
        color=PAL["text_dark"], y=0.975,
    )

    # ── pre-build running lists for scalar plots ──
    xs_shown:     list = []
    thr_shown:    list = []
    cmp_shown:    list = []

    # ── state for curve plots ──
    state = dict(frame=-1)

    # ─── animation update function ───────────────────────────────────────────
    def update(frame: int):
        if frame >= N:
            return

        state["frame"] = frame

        # header
        title_text.set_text(
            "LDPC Code Optimization  —  Round {:d}".format(iter_nums[frame])
        )

        # ── H matrix ──
        draw_H_matrix(ax_H, H_mats[frame])

        # ── scalar plots ──
        xs_shown.append(iter_nums[frame])
        thr_shown.append(thresholds[frame])
        cmp_shown.append(complexities[frame])

        update_scalar(ax_thr, xs_shown, thr_shown,
                      PAL["threshold"], PAL["threshold_pt"],
                      "Threshold (dB)", "Decoding Threshold")
        update_scalar(ax_cmp, xs_shown, cmp_shown,
                      PAL["complexity"], PAL["complexity_pt"],
                      "Complexity index", "Normalised Complexity")

        # ── curve plots ──
        initial_fer  = fer_curves[0]
        initial_avg  = avgit_curves[0]

        if frame == 0:
            # only initial state
            update_curves(ax_fer, ebn0, initial_fer, [], None,
                          "Frame Error Rate (FER)", "FER")
            update_curves(ax_avg, ebn0, initial_avg, [], None,
                          "Average Iterations (AvgIT)", "AvgIT", log_scale=False)
        else:
            # history = states 1 … frame-1 (exclude initial and current)
            history_fer  = fer_curves[1:frame]
            history_avg  = avgit_curves[1:frame]
            current_fer  = fer_curves[frame]
            current_avg  = avgit_curves[frame]

            update_curves(ax_fer, ebn0, initial_fer, history_fer, current_fer,
                          "Frame Error Rate (FER)", "FER")
            update_curves(ax_avg, ebn0, initial_avg, history_avg, current_avg,
                          "Average Iterations (AvgIT)", "AvgIT", log_scale=False)

        fig.canvas.draw_idle()

    anim = FuncAnimation(
        fig,
        update,
        frames=N,
        interval=UPDATE_INTERVAL_MS,
        repeat=False,
    )

    # ── add a progress bar / instruction strip at the bottom ──
    fig.text(
        0.5, 0.012,
        "Each frame represents one saved optimisation state  |  "
        f"Interval: {UPDATE_INTERVAL_MS / 1000:.1f} s  |  "
        "Optimisation rounds: 0 → 200",
        ha="center", va="bottom",
        fontsize=7.5, color=PAL["text_mid"], style="italic",
    )

    plt.show()


if __name__ == "__main__":
    main()
