"""
Step 5: Asymmetric Memory Extension for Path Asymmetry
======================================================
Introduces strategy-specific evidence retention: engagement-related
Dirichlet counts can decay faster than withdrawal-related counts (η_eng <
η_wdr).  Unlike symmetric decay (Step 4 null result), this structural
asymmetry can create genuine illness-duration hysteresis under full
same-variable reversal.

Clinical interpretation
-----------------------
Withdrawal strategies may become more habit-like or self-reinforcing,
retaining their evidential weight longer, while engagement evidence becomes
less accessible after prolonged illness.  The model gives this a
mechanistic reading: asymmetric η produces a ratchet in the Dirichlet
belief state that resists reversal as T_ill grows.

Non-monotonicity note
---------------------
The "Asymmetric mild" (η_eng=0.990, η_wdr=0.997) condition shows weak and
occasionally non-monotonic T_ill dependence.  This is expected: at mild
asymmetry levels, stochastic variation in Dirichlet trajectories can
dominate the small hysteresis signal.  The effect is clear and monotone
only for the stronger asymmetry condition (η_eng=0.985).  The mild
condition is retained for completeness but its T_ill slope should not be
over-interpreted.
"""

import csv
import os

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from core_functions import compute_efe_difference, policy_posterior_pi1
from psychopathology_regimes import recovery_regime


OUT = os.path.join(os.path.dirname(__file__), "Figs_psychopathology")
TABLE_OUT = os.path.join(os.path.dirname(__file__), "..", "outputs", "tables")
os.makedirs(OUT, exist_ok=True)
os.makedirs(TABLE_OUT, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "axes.titleweight": "bold",
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 8,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

REGIME = recovery_regime()

### DT ---> Sample sizes increased from original (na=60/40) to reduce noise
### DT ---> sufficiently to distinguish the mild asymmetry condition.
_NA_MAIN = 150    # main T_ill curves
_NA_HEAT = 100    # heatmap cells
_NA_THRESH = 80   # threshold sweep


def save_figure(fig, stem):
    fig.savefig(os.path.join(OUT, f"{stem}.png"))
    fig.savefig(os.path.join(OUT, f"{stem}.pdf"))


def run_agent_asymmetric_decay(p, a0, seed, gh, dch, Nh, gr, dcr, gf, dcf,
                               T_ill, rg, rdc, eta_eng=1.0, eta_wdr=1.0,
                               N_recovery=400):
    """Return late-window mean P(π₁) under strategy-specific decay."""
    rng = np.random.default_rng(seed)
    al1 = be1 = al2 = be2 = a0 / 2
    base = a0 / 2
    At = np.array([[p, 1-p], [1-p, p]])
    t_rev = Nh + T_ill
    N_total = t_rev + N_recovery
    late = []

    for t in range(N_total):
        if t < Nh:
            gt = gh
            dct = dch
        else:
            dt = t - Nh
            gt = max(gh - gr * dt, gf)
            dct = max(dch - dcr * dt, dcf)
        if t >= t_rev:
            gt = rg
            dct = rdc

        if eta_eng < 1.0 or eta_wdr < 1.0:
            ### DT ---> Engagement evidence and withdrawal evidence relax back
            ### DT ---> to baseline at different rates. This is the asymmetry
            ### DT ---> that symmetric decay (Step 4) cannot replicate.
            al1 = eta_eng * al1 + (1 - eta_eng) * base
            be1 = eta_eng * be1 + (1 - eta_eng) * base
            al2 = eta_wdr * al2 + (1 - eta_wdr) * base
            be2 = eta_wdr * be2 + (1 - eta_wdr) * base

        a = al1 / (al1 + be1)
        b = be2 / (al2 + be2)
        dG = compute_efe_difference(a, b, dct)
        Pp = policy_posterior_pi1(dG, gt)
        act = 0 if rng.random() < Pp else 1

        if act == 0:
            ob = 0 if rng.random() < At[0, 0] else 1
            if ob == 0:
                al1 += 1
            else:
                be1 += 1
        else:
            ob = 0 if rng.random() < At[0, 1] else 1
            if ob == 0:
                al2 += 1
            else:
                be2 += 1

        if t >= N_total - 100:
            late.append(Pp)

    return float(np.mean(late))


def run_condition(T_ill, eta_eng, eta_wdr, na=_NA_MAIN, seed_base=0,
                  dc_restored=None):
    p = REGIME["p"]
    a0 = REGIME["alpha_0"]
    gh = REGIME["gamma_healthy"]
    dch = REGIME["delta_c_healthy"]
    Nh = REGIME["N_healthy"]
    gr = REGIME["gamma_rate"]
    dcr = REGIME["delta_c_rate"]
    gf = REGIME["gamma_floor"]
    dcf = REGIME["delta_c_floor"]
    rg = gh
    rdc = dch if dc_restored is None else dc_restored

    vals = []
    for i in range(na):
        vals.append(run_agent_asymmetric_decay(
            p, a0, seed=seed_base + i, gh=gh, dch=dch, Nh=Nh,
            gr=gr, dcr=dcr, gf=gf, dcf=dcf, T_ill=T_ill,
            rg=rg, rdc=rdc, eta_eng=eta_eng, eta_wdr=eta_wdr
        ))
    return np.array(vals)


def find_dc_threshold(T_ill, eta_eng, eta_wdr, target=0.5):
    dc_grid = np.linspace(-0.02, 0.18, 17)
    meds = []
    for idx, dc in enumerate(dc_grid):
        vals = run_condition(
            T_ill, eta_eng, eta_wdr, na=_NA_THRESH,
            seed_base=40000 + T_ill * 50 + idx * 100,
            dc_restored=dc
        )
        meds.append(np.median(vals))
    meds = np.array(meds)
    above = np.where(meds >= target)[0]
    if len(above) == 0:
        return np.nan
    return float(dc_grid[above[0]])


def fig_asymmetric_hysteresis():
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.5))
    T_ills = [0, 50, 100, 200, 400, 600]
    ### DT ---> Labels include cumulative retention at T_ill=200 so the reader
    ### DT ---> sees that per-trial η close to 1 still implies large cumulative decay:
    ### DT ---> η=0.985 → 0.985^200 ≈ 0.05 (5% retained); η=0.997 → 0.55 (55% retained)
    conditions = [
        ('No decay', 1.0, 1.0, 'gray', '--'),
        ('Symmetric decay (η=0.992)', 0.992, 0.992, '#43A047', '-.'),
        (r'Asymmetric mild ($\eta_e$=0.990, $\eta_w$=0.997)', 0.990, 0.997, '#FB8C00', '-'),
        (r'Asymmetric strong ($\eta_e$=0.985, $\eta_w$=0.997)', 0.985, 0.997, '#8E24AA', '-'),
    ]

    ax = axes[0]
    summary_rows = []
    for idx, (label, eta_eng, eta_wdr, color, ls) in enumerate(conditions):
        meds, q1s, q3s = [], [], []
        for ti, T_ill in enumerate(T_ills):
            vals = run_condition(T_ill, eta_eng, eta_wdr, na=_NA_MAIN,
                                 seed_base=10000 + idx * 10000 + ti * 500)
            med = np.median(vals)
            q1, q3 = np.percentile(vals, [25, 75])
            meds.append(med)
            q1s.append(q1)
            q3s.append(q3)
            summary_rows.append([label, T_ill, eta_eng, eta_wdr, med, q1, q3])

        ax.plot(T_ills, meds, marker='o', linestyle=ls, color=color,
                ms=5, lw=2, label=label)
        ### DT ---> IQR band — helps distinguish noise from genuine T_ill trend
        ax.fill_between(T_ills, q1s, q3s, alpha=0.15, color=color)

    ax.axhline(y=0.5, color='gray', ls='--', alpha=0.5, label='Recovery threshold')
    ax.set_xlabel('Illness residence time ($T_\\mathrm{ill}$)')
    ax.set_ylabel(f'Median (IQR) late P(engaged) after full reversal\n(n = {_NA_MAIN} agents per cell)')
    ax.set_title('(a) Recovery declines with illness duration\n'
                 r'only when memory is asymmetric ($\eta_e < \eta_w$)')
    ax.legend(fontsize=7.5)
    ax.set_ylim(0, 1.05)

    ax = axes[1]
    eta_eng_grid = [1.000, 0.995, 0.990, 0.985, 0.980]
    T_grid = [0, 50, 100, 200, 400, 600]
    heat = np.zeros((len(eta_eng_grid), len(T_grid)))
    eta_wdr = 0.997
    for ei, eta_eng in enumerate(eta_eng_grid):
        for ti, T_ill in enumerate(T_grid):
            vals = run_condition(T_ill, eta_eng, eta_wdr, na=_NA_HEAT,
                                 seed_base=20000 + ei * 2000 + ti * 200)
            heat[ei, ti] = np.median(vals)
    im = ax.imshow(heat, aspect='auto', cmap='RdYlGn', vmin=0.2, vmax=1.0, origin='lower')
    ax.set_xticks(range(len(T_grid)))
    ax.set_xticklabels(T_grid)
    ax.set_yticks(range(len(eta_eng_grid)))
    ### DT ---> Show cumulative retention at T_ill=200 alongside η value so
    ### DT ---> readers see the per-trial and cumulative decay simultaneously.
    cum_labels = [
        rf'$\eta_e$={v:.3f}  ({v**200:.2f} retained at $T$=200)'
        for v in eta_eng_grid
    ]
    ax.set_yticklabels(cum_labels, fontsize=7.5)
    ax.set_xlabel('$T_\\mathrm{ill}$')
    ax.set_title(r'(b) Stronger $\eta_e$ asymmetry progressively lowers recovery'
                 '\n' + rf'($\eta_w = {eta_wdr:.3f}$, cumul. retained ≈0.55 at $T$=200; n = {_NA_HEAT})')
    for ei in range(len(eta_eng_grid)):
        for ti in range(len(T_grid)):
            val = heat[ei, ti]
            col = 'white' if val < 0.45 else 'black'
            ax.text(ti, ei, f'{val:.2f}', ha='center', va='center', fontsize=8, color=col)
    cb = plt.colorbar(im, ax=ax)
    cb.set_label('Median late P(engaged)')

    ax = axes[2]
    threshold_curves = [
        (r'Symmetric ($\eta_e$=0.992, $\eta_w$=0.992)', 0.992, 0.992, '#43A047'),
        (r'Asymmetric mild ($\eta_e$=0.990, $\eta_w$=0.997)', 0.990, 0.997, '#FB8C00'),
        (r'Asymmetric strong ($\eta_e$=0.985, $\eta_w$=0.997)', 0.985, 0.997, '#8E24AA'),
    ]
    for label, eta_eng, eta_wdr, color in threshold_curves:
        dc_needed = [find_dc_threshold(T_ill, eta_eng, eta_wdr) for T_ill in T_ills]
        ax.plot(T_ills, dc_needed, marker='o', color=color, lw=2, label=label)
    ax.axhline(y=REGIME["delta_c_healthy"], color='gray', ls=':', alpha=0.6,
               label=r'Healthy $\Delta c$')
    ax.set_xlabel('Illness residence time ($T_\\mathrm{ill}$)')
    ax.set_ylabel(r'Min. $\Delta c$ needed for recovery (median P ≥ 0.5)')
    ax.set_title(f'(c) Recovery threshold rises with illness duration\n'
                 f'only under asymmetric retention (n = {_NA_THRESH} per cell)')
    ax.legend(fontsize=8)
    ax.set_ylim(-0.02, 0.19)

    fig.suptitle(
        'Strategy-specific evidence retention creates illness-duration-dependent path asymmetry\n'
        r'(null result under symmetric decay, Step 4; mild asymmetry is weak — interpret with caution)',
        fontsize=12, y=1.04
    )
    plt.tight_layout()
    save_figure(fig, "fig_S5_asymmetric_memory")
    plt.close()
    print("Saved fig_S5_asymmetric_memory.[png|pdf]")

    with open(os.path.join(TABLE_OUT, "asymmetric_memory_summary.csv"), "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["condition", "T_ill", "eta_eng", "eta_wdr",
                          "median_late_p_engaged", "q1", "q3"])
        writer.writerows(summary_rows)
    print("Saved asymmetric_memory_summary.csv")


def run_all():
    """Entry point called by main.py."""
    print("Running Step 5: asymmetric memory extension...")
    fig_asymmetric_hysteresis()


if __name__ == "__main__":
    print("=" * 60)
    print("Step 5: Asymmetric Memory Extension")
    print("=" * 60)
    fig_asymmetric_hysteresis()
    print("=" * 60)
    print(f"Figure output: {OUT}")
    print(f"Table output: {TABLE_OUT}")
    print("=" * 60)
