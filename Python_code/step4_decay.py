"""
Symmetric Evidence Decay Does Not Create Regime-2 Hysteresis
=============================================================
Single-parameter extension: at each trial, all Dirichlet concentrations
decay toward the prior baseline:
    alpha_i -> eta * alpha_i + (1-eta) * alpha_0/2
 
This makes recent evidence count more than distant evidence.  The
original hypothesis was that an agent who has been withdrawn for a long
time would gradually FORGET engagement, creating an ambiguity asymmetry
that resists full same-variable reversal.
 
Result: symmetric decay does NOT produce this effect.  Because the
decay applies identically to all concentration parameters, it erases
illness-phase evidence during recovery just as effectively as it erases
health-phase evidence during illness.  Under full reversal with adequate
recovery time, faster decay can actually HELP recovery by making the
agent more responsive to the restored positive preference field.
 
This is a structural finding: the minimal two-state model with symmetric
Dirichlet learning and symmetric decay cannot produce asymmetric
hysteresis.  Genuine Regime-2 trapping would require strategy-dependent
decay rates, domain-specific forgetting, or other structural asymmetries
between the two strategies.
"""

import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, sys
from core_functions import compute_efe_difference, policy_posterior_pi1
from psychopathology_regimes import recovery_regime

OUT = os.path.join(os.path.dirname(__file__), "Figs_psychopathology")
os.makedirs(OUT, exist_ok=True)

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


def save_figure(fig, stem):
    fig.savefig(os.path.join(OUT, f"{stem}.png"))
    fig.savefig(os.path.join(OUT, f"{stem}.pdf"))


def run_agent_decay(p, a0, seed, gh, dch, Nh, gr, dcr, gf, dcf,
                    T_ill, rg, rdc, eta=1.0, N_recovery=400):
    """Agent with evidence decay.  Returns late-window mean P(pi1)."""
    rng = np.random.default_rng(seed)
    al1=be1=al2=be2=a0/2; base=a0/2
    At = np.array([[p,1-p],[1-p,p]])
    t_rev = Nh + T_ill
    N_total = t_rev + N_recovery
    Pp_late = []
 
    for t in range(N_total):
        # Schedule
        if t < Nh: gt=gh; dct=dch
        else:
            dt=t-Nh; gt=max(gh-gr*dt,gf); dct=max(dch-dcr*dt,dcf)
        if t >= t_rev: gt=rg; dct=rdc
 
        # Evidence decay (applied every trial)
        if eta < 1.0:
            al1 = eta*al1 + (1-eta)*base
            be1 = eta*be1 + (1-eta)*base
            al2 = eta*al2 + (1-eta)*base
            be2 = eta*be2 + (1-eta)*base
 
        # Policy
        a=al1/(al1+be1); b=be2/(al2+be2)
        dG = compute_efe_difference(a, b, dct)
        Pp = policy_posterior_pi1(dG, gt)
        act = 0 if rng.random() < Pp else 1
 
        # Environment + Dirichlet update
        if act==0:
            ob = 0 if rng.random() < At[0,0] else 1
            if ob==0: al1+=1
            else: be1+=1
        else:
            ob = 0 if rng.random() < At[0,1] else 1
            if ob==0: al2+=1
            else: be2+=1
 
        if t >= N_total - 100:
            Pp_late.append(Pp)
 
    return np.mean(Pp_late)
 
 
# ====================================================================
# FIGURE D1: Recovery under decay — does T_ill dependence emerge?
# ====================================================================
def fig_hysteresis_vs_till():
    """
    Panel (a): Full reversal recovery as f(T_ill) for several eta values.
    Panel (b): Same, with weaker Dc (0.05).
    Panel (c): Heatmap of recovery across (eta, T_ill).
 
    Result: symmetric decay does not produce T_ill-dependent trapping.
    Faster decay can facilitate recovery by erasing illness memories.
    """
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.5))
    p=REGIME["p"]; a0=REGIME["alpha_0"]; gh=REGIME["gamma_healthy"]
    Nh=REGIME["N_healthy"]; gr=REGIME["gamma_rate"]; gf=REGIME["gamma_floor"]; na=80
    T_ills = [0, 50, 100, 200, 300, 400, 600]
    etas = [1.0, 0.998, 0.995, 0.990, 0.980]
    eta_colours = ['gray', '#2196F3', '#4CAF50', '#FF9800', '#E91E63']
 
    for panel, (dc, ax) in enumerate(zip([0.10, 0.05], axes[:2])):
        dcr = dc/200; dcf = -dc
        for eta, col in zip(etas, eta_colours):
            hl = np.log(2)/np.log(1/eta) if eta<1 else float('inf')
            meds = []
            for Ti in T_ills:
                vals = []
                for i in range(na):
                    Pp = run_agent_decay(p, a0, seed=80000+panel*10000+i,
                        gh=gh, dch=dc, Nh=Nh, gr=gr, dcr=dcr, gf=gf, dcf=dcf,
                        T_ill=Ti, rg=gh, rdc=dc, eta=eta)
                    vals.append(Pp)
                meds.append(np.median(vals))
            label = f'η={eta}' if eta==1.0 else f'η={eta} (hl={hl:.0f})'
            ax.plot(T_ills, meds, 'o-', color=col, ms=5, lw=1.5, label=label)
        ax.axhline(y=0.5, color='gray', ls='--', alpha=0.5)
        ax.set_xlabel('Illness residence time (T_ill)')
        ax.set_ylabel('Median late P(engaged) after full reversal')
        ax.set_title(f'({"a" if panel==0 else "b"}) Full reversal, '
                     rf'$\Delta c_0 = {dc}$')
        ax.legend(fontsize=7); ax.set_ylim(0, 1.05)
 
    # Panel (c): Heatmap
    ax = axes[2]
    dc = 0.08; dcr=dc/200; dcf=-dc
    eta_grid = [1.0, 0.998, 0.995, 0.990, 0.985, 0.980]
    T_grid = [0, 50, 100, 200, 400, 800]
    rec = np.zeros((len(eta_grid), len(T_grid)))
    for ei, eta in enumerate(eta_grid):
        for ti, Ti in enumerate(T_grid):
            vals = []
            for i in range(na):
                Pp = run_agent_decay(p, a0, seed=90000+ei*1000+ti*100+i,
                    gh=gh, dch=dc, Nh=Nh, gr=gr, dcr=dcr, gf=gf, dcf=dcf,
                    T_ill=Ti, rg=gh, rdc=dc, eta=eta)
                vals.append(Pp)
            rec[ei, ti] = np.median(vals)
    im = ax.imshow(rec, aspect='auto', cmap='RdYlGn', vmin=0.2, vmax=1.0,
                   origin='lower')
    ax.set_xticks(range(len(T_grid)))
    ax.set_xticklabels(T_grid)
    ax.set_yticks(range(len(eta_grid)))
    ax.set_yticklabels([f'{e:.3f}' for e in eta_grid])
    ax.set_xlabel('T_ill'); ax.set_ylabel(r'Decay rate $\eta$')
    for ei in range(len(eta_grid)):
        for ti in range(len(T_grid)):
            v = rec[ei, ti]
            c = 'white' if v < 0.45 else 'black'
            ax.text(ti, ei, f'{v:.2f}', ha='center', va='center', fontsize=8, color=c)
    cb = plt.colorbar(im, ax=ax); cb.set_label('Median late P(engaged)')
    ax.set_title(rf'(c) Recovery heatmap ($\Delta c_0={dc}$)')
 
    fig.suptitle('Symmetric evidence decay produces no illness-duration-dependent hysteresis\n'
                 r'(null result: full same-variable reversal, $\gamma$ and $\Delta c$ both restored)',
                 fontsize=13, y=1.04)
    plt.tight_layout()
    save_figure(fig, "fig_D1_decay_hysteresis")
    plt.close()
    print("Saved fig_D1_decay_hysteresis.[png|pdf]")
 
 
# ====================================================================
# FIGURE D2: Contrast — with vs without decay
# ====================================================================
def fig_decay_contrast():
    """
    Direct comparison: the same reversal protocol with and without decay.
    Result: symmetric decay does not hinder recovery; faster decay can
    facilitate it by erasing illness-phase evidence during recovery.
    """
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    p=REGIME["p"]; a0=REGIME["alpha_0"]; gh=REGIME["gamma_healthy"]
    Nh=REGIME["N_healthy"]; gr=REGIME["gamma_rate"]; gf=REGIME["gamma_floor"]; na=80
    dc=0.08; dcr=dc/200; dcf=-dc
    T_ills = [0, 50, 100, 200, 400, 800]
 
    # Panel (a): bar comparison at T_ill=400
    ax = axes[0]
    Ti = 400
    conditions = [
        ('No decay\n(η=1.0)', 1.0, 'gray'),
        ('Slow decay\n(η=0.995)', 0.995, '#4CAF50'),
        ('Moderate decay\n(η=0.990)', 0.990, '#FF9800'),
        ('Fast decay\n(η=0.980)', 0.980, '#E91E63'),
    ]
    for ci, (name, eta, col) in enumerate(conditions):
        vals = []
        for i in range(na):
            Pp = run_agent_decay(p, a0, seed=95000+ci*1000+i,
                gh=gh, dch=dc, Nh=Nh, gr=gr, dcr=dcr, gf=gf, dcf=dcf,
                T_ill=Ti, rg=gh, rdc=dc, eta=eta)
            vals.append(Pp)
        vals = np.array(vals)
        med = np.median(vals); q1,q3 = np.percentile(vals,[25,75])
        ax.bar(ci, med, color=col, alpha=0.8, edgecolor='black')
        ax.errorbar(ci, med, yerr=[[med-q1],[q3-med]], fmt='none',
                    capsize=6, color='black', lw=1.5)
    ax.axhline(y=0.5, color='gray', ls='--', alpha=0.5)
    ax.set_xticks(range(4))
    ax.set_xticklabels([c[0] for c in conditions], fontsize=9)
    ax.set_ylabel('Median late P(engaged)')
    ax.set_title(f'(a) Full reversal at T_ill={Ti}\n'
                 rf'($\Delta c_0={dc}$, $\gamma$=16)')
    ax.set_ylim(0, 1.05)
 
    # Panel (b): T_ill curves, no-decay vs moderate-decay
    ax = axes[1]
    for eta, col, ls, label in [(1.0, 'gray', '--', 'No decay (η=1.0)'),
                                 (0.990, '#FF9800', '-', 'With decay (η=0.990)')]:
        meds = []
        for Ti in T_ills:
            vals = []
            for i in range(na):
                Pp = run_agent_decay(p, a0, seed=96000+int(eta*1000)+i,
                    gh=gh, dch=dc, Nh=Nh, gr=gr, dcr=dcr, gf=gf, dcf=dcf,
                    T_ill=Ti, rg=gh, rdc=dc, eta=eta)
                vals.append(Pp)
            meds.append(np.median(vals))
        ax.plot(T_ills, meds, marker='o', linestyle=ls, color=col, ms=6, lw=2, label=label)
    ax.axhline(y=0.5, color='gray', ls='--', alpha=0.5, label='Recovery threshold')
    ax.set_xlabel('Illness residence time (T_ill)')
    ax.set_ylabel('Median late P(engaged) after full reversal')
    ax.set_title(f'(b) T_ill dependence under decay:\n'
                 'symmetric decay does not create trapping')
    ax.legend(fontsize=9); ax.set_ylim(0, 1.05)
 
    fig.suptitle('Symmetric Decay Is Not the Missing Mechanism for Regime-2 Hysteresis\n'
                 '(decay erases illness memories during recovery as effectively as '
                 'health memories during illness)',
                 fontsize=13, y=1.05)
    plt.tight_layout()
    save_figure(fig, "fig_D2_decay_contrast")
    plt.close()
    print("Saved fig_D2_decay_contrast.[png|pdf]")
 
 
def run_all():
    """
    Entry point called by main.py.

    Produces only fig_D1_decay_hysteresis.png (the informative negative
    result that symmetric decay does not create T_ill-dependent trapping).
    fig_D2_decay_contrast is retained as a callable function for ad-hoc
    use but is not part of the main publication sweep.
    """
    print("Running Step 4: symmetric evidence decay (negative result)...")
    fig_hysteresis_vs_till()


if __name__ == "__main__":
    print("=" * 60)
    print("Step 4: Symmetric Evidence Decay (Informative Negative Result)")
    print("=" * 60)
    print("\n1. Recovery under decay across T_ill and eta...")
    fig_hysteresis_vs_till()
    print("\n2. Contrast: with vs without decay (supplementary, not in main sweep)...")
    fig_decay_contrast()
    print("\n" + "=" * 60)
