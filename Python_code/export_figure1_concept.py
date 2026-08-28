from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


CODEBASE_ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT_DIR = CODEBASE_ROOT / "outputs" / "manuscript_figures"

ACCENT = "#176B87"
INK = "#27343A"
MID_GREY = "#6F7A80"
BOX_EDGE = "#879096"
BOX_FILL = "#F5F7F8"
BAND_FILL = "#DCECEF"


def add_box(ax, xy, width, height, text, fontsize=9):
    box = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.025,rounding_size=0.025",
        linewidth=1.2,
        facecolor=BOX_FILL,
        edgecolor=BOX_EDGE,
    )
    ax.add_patch(box)
    ax.text(
        xy[0] + width / 2,
        xy[1] + height / 2,
        text,
        ha="center",
        va="center",
        fontsize=fontsize,
        color=INK,
    )
    return box


def add_arrow(ax, start, end):
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=12,
        linewidth=1.35,
        color=INK,
    )
    ax.add_patch(arrow)


def add_question_box(ax, xy, width, height, number, heading, question):
    box = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.025,rounding_size=0.025",
        linewidth=1.2,
        facecolor=BOX_FILL,
        edgecolor=BOX_EDGE,
    )
    ax.add_patch(box)
    ax.text(
        xy[0] + width * 0.16,
        xy[1] + height / 2,
        f"{number}\n{heading}",
        ha="center",
        va="center",
        fontsize=8.5,
        fontweight="bold",
        color=INK,
        linespacing=1.15,
    )
    ax.plot(
        [xy[0] + width * 0.31, xy[0] + width * 0.31],
        [xy[1] + height * 0.18, xy[1] + height * 0.82],
        color=BOX_EDGE,
        linewidth=0.8,
    )
    ax.text(
        xy[0] + width * 0.655,
        xy[1] + height / 2,
        question,
        ha="center",
        va="center",
        fontsize=8.2,
        color=INK,
        linespacing=1.2,
    )
    return box


def panel_label(ax, label, title):
    ax.text(0.0, 1.03, label, transform=ax.transAxes, fontsize=13, fontweight="bold")
    ax.text(0.12, 1.03, title, transform=ax.transAxes, fontsize=11, fontweight="bold")


def make_figure():
    fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.6))
    fig.patch.set_facecolor("white")

    ax = axes[0, 0]
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    panel_label(ax, "(a)", "Learning-action feedback")
    add_box(ax, (0.04, 0.60), 0.34, 0.18, "Policy selection\nengage or withdraw")
    add_box(ax, (0.62, 0.60), 0.34, 0.18, "Action-dependent\nobservations")
    add_box(ax, (0.62, 0.18), 0.34, 0.18, "Evidence accumulation\nupdates the likelihood\nmatrix", 8.5)
    add_box(ax, (0.04, 0.18), 0.34, 0.18, "Updated expected\nfree energy")
    add_arrow(ax, (0.38, 0.69), (0.62, 0.69))
    add_arrow(ax, (0.79, 0.60), (0.79, 0.36))
    add_arrow(ax, (0.62, 0.27), (0.38, 0.27))
    add_arrow(ax, (0.21, 0.36), (0.21, 0.60))
    ax.text(0.5, 0.03, "Behaviour changes evidence; evidence changes future behaviour", ha="center", fontsize=9, color=INK)

    ax = axes[0, 1]
    panel_label(ax, "(b)", "Emergent transition structure")
    z = np.linspace(-1.25, 1.25, 700)
    field = z**3 - z
    fold_z = 1 / np.sqrt(3)
    stable_low = z <= -fold_z
    unstable = np.abs(z) < fold_z
    stable_high = z >= fold_z
    ax.plot(field[stable_high], z[stable_high], color=ACCENT, linewidth=2.8, label="Engaged branch")
    ax.plot(
        field[stable_low],
        z[stable_low],
        color=INK,
        linewidth=2.8,
        linestyle=(0, (7, 2.2, 1.5, 2.2)),
        label="Withdrawn branch",
    )
    ax.plot(field[unstable], z[unstable], color=MID_GREY, linestyle=":", linewidth=1.5)
    fold_field = 2 / (3 * np.sqrt(3))
    ax.axvspan(-fold_field, fold_field, color=BAND_FILL, alpha=0.65)
    ax.text(
        0.18,
        0.42,
        "Multistable\nregion",
        ha="center",
        va="center",
        fontsize=8.5,
        color=INK,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.92, "pad": 1.8},
    )
    ax.set_xlabel("Preference field Δc")
    ax.set_ylabel("Order parameter z")
    ax.set_xlim(-0.75, 0.75)
    ax.set_ylim(-1.25, 1.25)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.text(
        0.03,
        0.96,
        "Folds permit abrupt shifts\nand hysteresis",
        transform=ax.transAxes,
        va="top",
        fontsize=8.0,
        color=INK,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9, "pad": 1.5},
    )

    ax = axes[1, 0]
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    panel_label(ax, "(c)", "What the model tests")
    add_question_box(
        ax,
        (0.04, 0.68),
        0.92,
        0.19,
        "1",
        "ONSET",
        "Can gradual adversity\ntrigger abrupt withdrawal?",
    )
    add_question_box(
        ax,
        (0.04, 0.405),
        0.92,
        0.19,
        "2",
        "RECOVERY",
        "Which model control most\nstrongly supports re-engagement?",
    )
    add_question_box(
        ax,
        (0.04, 0.13),
        0.92,
        0.19,
        "3",
        "EMPIRICAL",
        "Can dense clinical data recover\na previously identified regime shift?",
    )

    ax = axes[1, 1]
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    panel_label(ax, "(d)", "Evidence hierarchy")
    add_box(ax, (0.06, 0.74), 0.88, 0.15, "MODEL\nΔc, γ, z and transition dynamics", 8.7)
    add_box(ax, (0.06, 0.52), 0.88, 0.15, "PRIMARY TRAJECTORY\nWichers: published day-127 transition recovered", 8.3)
    add_box(ax, (0.06, 0.30), 0.88, 0.15, "SECONDARY VALIDATION\nAction history and sparse trajectories", 8.3)
    add_box(ax, (0.06, 0.08), 0.88, 0.15, "DIRECT TEST\nManipulated actions + individual model fitting", 8.3)

    for panel in axes.ravel():
        panel.set_facecolor("white")

    plt.tight_layout(rect=(0.02, 0.02, 0.98, 0.98), h_pad=2.0, w_pad=1.8)
    return fig


def main():
    MANUSCRIPT_DIR.mkdir(parents=True, exist_ok=True)
    fig = make_figure()
    fig.savefig(MANUSCRIPT_DIR / "Fig_1.pdf", bbox_inches="tight")
    fig.savefig(MANUSCRIPT_DIR / "Fig_1.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Figure 1 regenerated as a theory-led conceptual diagram.")


if __name__ == "__main__":
    main()
