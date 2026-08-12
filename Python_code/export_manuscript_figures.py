import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import sim_psychopathology as step1
import step2_hysteresis as step2
import step4_decay as step4
import step5_asymmetric_memory as step5
import export_figure1_concept as concept


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT_DIR = ROOT / "Manuscript"


def _copy_figure(source_stem, destination_stem):
    for suffix in [".pdf", ".png"]:
        shutil.copy2(
            Path(f"{source_stem}{suffix}"),
            MANUSCRIPT_DIR / f"{destination_stem}{suffix}",
        )


def export_fig1_to_fig4():
    concept.main()
    step1.plot_onset_with_theory()
    step2.fig_clinical_prediction()
    step2.fig_recovery_boundary()

    _copy_figure(Path(step1.OUT) / "fig_P2_onset_revised", "Fig_2")
    _copy_figure(Path(step2.OUT) / "fig_S2_clinical_prediction", "Fig_3")
    _copy_figure(Path(step2.OUT) / "fig_S2_recovery_boundary", "Fig_4")


def export_fig5():
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.7))

    p = step4.REGIME["p"]
    a0 = step4.REGIME["alpha_0"]
    gh = step4.REGIME["gamma_healthy"]
    Nh = step4.REGIME["N_healthy"]
    gr = step4.REGIME["gamma_rate"]
    gf = step4.REGIME["gamma_floor"]
    na = 80
    T_ills = [0, 50, 100, 200, 300, 400, 600]
    etas = [1.0, 0.998, 0.995, 0.990, 0.980]
    eta_colours = ["#7a7a7a", "#4c78a8", "#2ca25f", "#f28e2b", "#c44e52"]

    ax = axes[0]
    dc = 0.10
    dcr = dc / 200
    dcf = -dc
    for eta, col in zip(etas, eta_colours):
        meds = []
        for ti, t_ill in enumerate(T_ills):
            vals = []
            for i in range(na):
                vals.append(
                    step4.run_agent_decay(
                        p, a0, seed=80000 + ti * 100 + i,
                        gh=gh, dch=dc, Nh=Nh, gr=gr, dcr=dcr, gf=gf, dcf=dcf,
                        T_ill=t_ill, rg=gh, rdc=dc, eta=eta
                    )
                )
            meds.append(np.median(vals))
        ax.plot(T_ills, meds, "o-", color=col, ms=4.5, lw=1.6, label=f"η={eta:.3f}" if eta < 1 else "η=1.000")
    ax.axhline(y=0.5, color="#9e9e9e", ls="--", lw=0.8)
    ax.set_xlabel("Adverse-schedule exposure after healthy phase")
    ax.set_ylabel("Median late P(engaged) after full reversal")
    ax.set_title(r"(a) Full reversal, $\Delta c_0 = 0.10$")
    ax.set_ylim(0.45, 1.0)
    ax.legend(
        frameon=False, fontsize=7, loc="lower left",
        bbox_to_anchor=(0.02, 0.12), borderaxespad=0.0
    )

    ax = axes[1]
    threshold_curves = [
        ("Symmetric decay (η=0.995)", 0.995, 0.995, "#2ca25f"),
        (r"Asymmetric mild ($\eta_{eng}$=0.990, $\eta_{wdr}$=0.997)", 0.990, 0.997, "#f28e2b"),
        (r"Asymmetric strong ($\eta_{eng}$=0.985, $\eta_{wdr}$=0.997)", 0.985, 0.997, "#7b3294"),
    ]
    t_ills = [0, 50, 100, 200, 300, 400, 500, 600]
    for label, eta_eng, eta_wdr, color in threshold_curves:
        dc_needed = [step5.find_dc_threshold(t_ill, eta_eng, eta_wdr) for t_ill in t_ills]
        ax.plot(t_ills, dc_needed, marker="o", color=color, ms=4.5, lw=1.8, label=label)
    ax.set_xlabel("Adverse-schedule exposure after healthy phase")
    ax.set_ylabel(r"Min. restored $\Delta c$ needed for recovery")
    ax.set_title(r"(b) Recovery thresholds with asymmetric memory")
    ax.set_ylim(-0.005, 0.105)
    ax.legend(frameon=False, fontsize=6.7, loc="upper left")
    ax.text(
        0.98, 0.97, r"Healthy $\Delta c$ = 0.15",
        transform=ax.transAxes, ha="right", va="top",
        fontsize=6.8, color="#6e6e6e"
    )

    plt.tight_layout()
    fig.savefig(MANUSCRIPT_DIR / "Fig_5.pdf")
    fig.savefig(MANUSCRIPT_DIR / "Fig_5.png", dpi=300)
    plt.close(fig)


def main():
    MANUSCRIPT_DIR.mkdir(parents=True, exist_ok=True)
    export_fig1_to_fig4()
    export_fig5()


if __name__ == "__main__":
    main()
