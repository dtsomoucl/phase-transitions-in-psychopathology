"""
Step 2: Recovery Boundary Under Partial and Full Control Restoration
====================================================================
Building on the Step 1 consolidation, we map the recovery landscape
in the (gamma_restored, Dc_restored) control space.

Key findings:
  - In the minimal two-state model with permanent Dirichlet evidence,
    evidence-driven hysteresis is WEAK (O(0.001) in ambiguity difference
    once both strategies have been explored for >50 trials).
    This is a substantive negative result.
  - The dominant barrier to recovery is the ONGOING ADVERSE PREFERENCE
    FIELD (low Dc), not accumulated maladaptive memory.
  - Restoring Dc is necessary for recovery; restoring gamma alone fails.

This Step establishes:
  (a) A heuristic approximation of the recovery boundary in (γ, Δc)
      space, using a reduced model of the post-illness Dirichlet state,
      validated against the full stochastic simulator.
  (b) The weak T_ill dependence — honestly documented.
  (c) The central result: recovery requires restoring the preference
      field, not just cognitive precision.
  (d) The setup for Step 3: what if Δc CANNOT be restored directly?

Regime disclosure
-----------------
All simulations in this file use RECOVERY_REGIME (α₀=2, γ_rate=0.05),
which has a lower prior precision and faster γ erosion than ONSET_REGIME
(α₀=40, γ_rate=0.005).  The field-dominance conclusion is robust to this
choice — see fig_S2_cross_regime.png (step2_robustness.py).
"""

import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, sys
from core_functions import (compute_efe_difference, policy_posterior_pi1,
    binary_entropy, coupling_function, find_gamma_c, find_fixed_points,
    self_consistency_rhs)
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
 
 
# ── Fast agent for parameter sweeps ──────────────────────────────────
def sweep_agent(p, a0, N, seed, gh, dch, Nh, gr, dcr, gf, dcf,
                t_rev, rg, rdc, N_recovery=300):
    """Run onset + reversal, return mean P(pi1) over last 100 recovery trials."""
    rng = np.random.default_rng(seed)
    al1=be1=al2=be2=a0/2; At=np.array([[p,1-p],[1-p,p]])
    N_total = t_rev + N_recovery
    late_sum = 0.0; late_count = 0
    for t in range(N_total):
        if t<Nh: gt=gh; dct=dch
        else:
            dt=t-Nh; gt=max(gh-gr*dt,gf); dct=max(dch-dcr*dt,dcf)
        if t >= t_rev:
            gt = rg; dct = rdc
        a=al1/(al1+be1); b=be2/(al2+be2)
        dG=compute_efe_difference(a,b,dct)
        Pp=policy_posterior_pi1(dG,gt)
        act=0 if rng.random()<Pp else 1
        if act==0:
            ob=0 if rng.random()<At[0,0] else 1
            if ob==0: al1+=1
            else: be1+=1
        else:
            ob=0 if rng.random()<At[0,1] else 1
            if ob==0: al2+=1
            else: be2+=1
        if t >= N_total - 100:
            late_sum += Pp; late_count += 1
    return late_sum / max(late_count, 1)
 
 
# ── Semi-analytical recovery boundary (reduced approximation) ────────
def analytical_recovery_Pp1(gamma, dc, a, b):
    """Given learned (a, b) and restored (gamma, dc), return P(pi1)."""
    dG = compute_efe_difference(a, b, dc)
    return policy_posterior_pi1(dG, gamma)
 
def post_illness_discriminabilities(p, a0, Nh, T_ill, frac_eng_healthy=0.85,
                                     frac_wd_ill=0.80):
    """Approximate learned a, b after healthy phase + illness.

    ### DT --> The healthy phase is allowed to contribute learning to both columns,
    ### DT --> so T_ill = 0 does not collapse the withdrawal discriminability to 0.5.
    """
    n1 = Nh * frac_eng_healthy
    n2 = Nh * max(1.0 - frac_eng_healthy, 0.0) + T_ill * frac_wd_ill
    a = (a0/2 + n1*p) / (a0 + n1)
    b = (a0/2 + n2*p) / (a0 + n2) if n2 > 0 else 0.5
    return a, b
 
 
# ====================================================================
# FIGURE S2a: Recovery boundary in (gamma, Dc) space
# ====================================================================
def fig_recovery_boundary():
    """
    Panel (a): Semi-analytical recovery boundary — contour of P(engaged)=0.5
               in (gamma_restored, Dc_restored) space, for several T_ill,
               using a reduced heuristic approximation of the post-illness state.
    Panel (b): Simulation verification of the boundary.
    Panel (c): Decomposition — why Dc matters more than T_ill.
    """
    fig, axes = plt.subplots(
        1, 3, figsize=(12.2, 4.6),
        gridspec_kw={"width_ratios": [1.35, 1.0, 1.0]}
    )
    p = REGIME["p"]; a0 = REGIME["alpha_0"]
    Nh = REGIME["N_healthy"]
    gh = REGIME["gamma_healthy"]; dch = REGIME["delta_c_healthy"]
    gr = REGIME["gamma_rate"]; dcr = REGIME["delta_c_rate"]
    gf = REGIME["gamma_floor"]; dcf = REGIME["delta_c_floor"]
 
    # ── Panel (a): Analytical boundary ──
    ax = axes[0]
    gamma_range = np.linspace(4, 20, 80)
    dc_range = np.linspace(-0.2, 0.5, 80)
    
    a_ref, b_ref = post_illness_discriminabilities(p, a0, 200, 200)
    ref_grid = np.zeros((len(dc_range), len(gamma_range)))
    for i, dc in enumerate(dc_range):
        for j, gam in enumerate(gamma_range):
            ref_grid[i, j] = analytical_recovery_Pp1(gam, dc, a_ref, b_ref)
    ax.contourf(
        gamma_range, dc_range, ref_grid,
        levels=[0.0, 0.5, 1.0],
        colors=["#f5e7e7", "#e6f4ea"],
        alpha=0.65
    )

    for T_ill, color, ls in [(0, '#7a7a7a', ':'), (100, '#4c78a8', '--'),
                              (200, '#2ca25f', '-'), (400, '#c44e52', '-.')]:
        a, b = post_illness_discriminabilities(p, a0, 200, T_ill)
        # Compute P(pi1) on grid
        Pp_grid = np.zeros((len(dc_range), len(gamma_range)))
        for i, dc in enumerate(dc_range):
            for j, gam in enumerate(gamma_range):
                Pp_grid[i, j] = analytical_recovery_Pp1(gam, dc, a, b)
        ax.contour(gamma_range, dc_range, Pp_grid, levels=[0.5],
                   colors=[color], linestyles=[ls], linewidths=2)
        # Label
        idx = np.argmin(np.abs(Pp_grid[:, -1] - 0.5))
        ax.annotate(f'T={T_ill}', xy=(gamma_range[-1], dc_range[idx]),
                    fontsize=8, color=color, ha='right')

    ax.set_xlabel(r'Restored $\gamma$')
    ax.set_ylabel(r'Restored $\Delta c$')
    ax.set_title('(a) Semi-analytical recovery boundary')
    ax.axhline(y=0, color='gray', ls='--', alpha=0.3)
    ax.text(4.8, 0.34, 'Recovery', fontsize=10, color='#1b7837',
            fontstyle='italic')
    ax.text(4.8, -0.13, 'Persistent withdrawal', fontsize=9, color='#a50f15',
            fontstyle='italic')
    ### DT ---> Add manual legend for contour lines (contour() does not
    ### DT ---> auto-label in this version of matplotlib)
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='gray',  linestyle=':',  lw=2, label='T_ill=0'),
        Line2D([0], [0], color='blue',  linestyle='--', lw=2, label='T_ill=100'),
        Line2D([0], [0], color='green', linestyle='-',  lw=2, label='T_ill=200'),
        Line2D([0], [0], color='red',   linestyle='-.', lw=2, label='T_ill=400'),
    ]
    ax.legend(handles=legend_elements, fontsize=8, loc='lower right')
 
    # ── Panel (b): Simulation verification ──
    ax = axes[1]
    na = 40
    dc_test = np.linspace(-0.1, 0.4, 11)
    gamma_test = [8.0, 12.0, 16.0]
    T_ill_sim = 200; t_rev = Nh + T_ill_sim

    ### DT ---> Seeds vary by both agent index i and dc value (di) to
    ### DT ---> avoid artificial correlation between Dc levels.
    for gi, (gam, marker, col) in enumerate(zip(gamma_test, ['o', 's', 'D'],
                                 ['steelblue', 'green', 'darkred'])):
        meds = []
        for di, dc in enumerate(dc_test):
            vals = []
            for i in range(na):
                Pp = sweep_agent(p, a0, None,
                    seed=10000 + gi * 10000 + di * 100 + i,
                    gh=gh, dch=dch, Nh=Nh, gr=gr, dcr=dcr,
                    gf=gf, dcf=dcf,
                    t_rev=t_rev, rg=gam, rdc=dc)
                vals.append(Pp)
            meds.append(np.median(vals))
        ax.plot(dc_test, meds, f'{marker}-', color=col, ms=5, lw=1.5,
                label=rf'$\gamma_{{res}}={gam}$')
 
    ax.axhline(y=0.5, color='gray', ls='--', alpha=0.5, label='Recovery threshold')
    ax.set_xlabel(r'Restored $\Delta c$')
    ax.set_ylabel('Median late P(engaged)')
    ax.set_title(f'(b) Simulation check (T = {T_ill_sim})')
    ax.legend(fontsize=8)
 
    # ── Panel (c): T_ill dependence (weak) ──
    ax = axes[2]
    T_ills = [50, 100, 200, 400, 600]
    # At Dc_restored = 0.08 (moderate), sweep T_ill
    dc_rest = 0.08
    for gam, col, ls in [(16.0, 'blue', '-'), (10.0, 'orange', '--'),
                          (7.0, 'red', ':')]:
        meds = []
        for Ti in T_ills:
            t_rev = Nh + Ti
            vals = []
            for i in range(na):
                Pp = sweep_agent(p, a0, None, seed=11000+i,
                    gh=gh, dch=dch, Nh=Nh, gr=gr, dcr=dcr,
                    gf=gf, dcf=dcf,
                    t_rev=t_rev, rg=gam, rdc=dc_rest)
                vals.append(Pp)
            meds.append(np.median(vals))
        ax.plot(T_ills, meds, 'o-', color=col, lw=1.5, ms=6,
                label=rf'$\gamma_{{res}}={gam},\ \Delta c_{{res}}={dc_rest}$')
 
    ax.axhline(y=0.5, color='gray', ls='--', alpha=0.5)
    ax.set_xlabel('Illness residence time (T_ill)')
    ax.set_ylabel('Median late P(engaged)')
    ax.set_title(r'(c) Weak illness-duration dependence')
    ax.legend(fontsize=7)
 
    fig.suptitle(
        'Recovery is governed primarily by the restored preference field',
        fontsize=12, y=1.03
    )
    plt.tight_layout()
    save_figure(fig, "fig_S2_recovery_boundary")
    plt.close()
    print("Saved fig_S2_recovery_boundary.[png|pdf]")
 
 
# ====================================================================
# FIGURE S2b: The clinical prediction — which variable matters?
# ====================================================================
def fig_clinical_prediction():
    """
    Clean comparison of four recovery strategies after fixed onset:
      1. Restore gamma only (Dc adverse)
      2. Restore Dc only (gamma stays low)
      3. Restore both (full same-variable reversal)
      4. No intervention (control)
    
    This directly answers: what is the critical variable for recovery?
    """
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.1))
    p = REGIME["p"]; a0 = REGIME["alpha_0"]; Nh = REGIME["N_healthy"]
    gh = REGIME["gamma_healthy"]; dch = REGIME["delta_c_healthy"]
    gr = REGIME["gamma_rate"]; dcr = REGIME["delta_c_rate"]
    gf = REGIME["gamma_floor"]; dcf = REGIME["delta_c_floor"]
    T_ill = 200; t_rev = Nh + T_ill
    na = 100; N_total = t_rev + 400
    
    conditions = {
        'No intervention': dict(rg=None, rdc=None),
        r'Restore $\gamma$ only': dict(rg=gh, rdc=None),
        r'Restore $\Delta c$ only': dict(rg=None, rdc=dch),
        r'Restore both': dict(rg=gh, rdc=dch),
    }
    colours = ['gray', 'steelblue', 'orange', 'darkgreen']
    
    # ── Panel (a): Recovery metric by condition ──
    ax = axes[0]
    for idx, (name, kw) in enumerate(conditions.items()):
        vals = []
        for i in range(na):
            h_rg = kw.get('rg'); h_rdc = kw.get('rdc')
            Pp = sweep_agent(p, a0, None, seed=12000+i,
                gh=gh, dch=dch, Nh=Nh, gr=gr, dcr=dcr,
                gf=gf, dcf=dcf,
                t_rev=t_rev if (h_rg or h_rdc) else 99999,
                rg=h_rg if h_rg else gf,
                rdc=h_rdc if h_rdc else dcf,
                N_recovery=400)
            vals.append(Pp)
        vals = np.array(vals)
        med = np.median(vals)
        q1, q3 = np.percentile(vals, [25, 75])
        ax.vlines(idx, q1, q3, color='black', lw=1.8, zorder=2)
        ax.scatter(idx, med, s=90, color=colours[idx], edgecolor='black', zorder=3)
 
    ax.set_xticks(range(4))
    ax.set_xticklabels(['Control', r'$\gamma$ only', r'$\Delta c$ only', 'Both'],
                       fontsize=10)
    ax.set_ylabel('Median late P(engaged)')
    ax.axhline(y=0.5, color='gray', ls='--', alpha=0.5)
    ax.set_title(f'(a) Recovery by intervention type (T = {T_ill})')
    ax.set_ylim(0, 1.05)
    
    # ── Panel (b): Trajectories for the four conditions ──
    ax = axes[1]
    from core_functions import compute_efe_difference, policy_posterior_pi1
    
    def full_trajectory(seed, rg, rdc, t_rev_actual):
        rng = np.random.default_rng(seed)
        al1=be1=al2=be2=a0/2; At=np.array([[p,1-p],[1-p,p]])
        Pp_trace = np.zeros(N_total)
        for t in range(N_total):
            if t<Nh: gt=gh; dct=dch
            else:
                dt=t-Nh; gt=max(gh-gr*dt,gf); dct=max(dch-dcr*dt,dcf)
            if t >= t_rev_actual:
                if rg is not None: gt=rg
                if rdc is not None: dct=rdc
            a=al1/(al1+be1); b=be2/(al2+be2)
            dG=compute_efe_difference(a,b,dct)
            Pp=policy_posterior_pi1(dG,gt)
            act=0 if rng.random()<Pp else 1
            if act==0:
                ob=0 if rng.random()<At[0,0] else 1
                if ob==0: al1+=1
                else: be1+=1
            else:
                ob=0 if rng.random()<At[0,1] else 1
                if ob==0: al2+=1
                else: be2+=1
            Pp_trace[t] = Pp
        return Pp_trace
    
    def causal_smooth(arr, w=20):
        out = np.empty_like(arr, dtype=float); cs = np.cumsum(arr)
        for t in range(len(arr)):
            ww = min(t+1, w); out[t] = (cs[t] - (cs[t-ww] if t>=ww else 0)) / ww
        return out
    seed = 12042  # fixed seed for reproducibility
    for name, kw, col in zip(
            ['Control', r'$\gamma$ only', r'$\Delta c$ only', 'Both'],
            [dict(rg=None,rdc=None), dict(rg=gh,rdc=None),
             dict(rg=None,rdc=dch), dict(rg=gh,rdc=dch)],
            colours):
        tr_rev = t_rev if (kw['rg'] or kw['rdc']) else 99999
        trace = full_trajectory(seed, kw['rg'], kw['rdc'], tr_rev)
        ax.plot(range(N_total), causal_smooth(trace, 15), color=col,
                lw=2 if 'Both' in name else 1.5,
                ls='-' if 'Both' in name else '--',
                alpha=0.9 if 'Both' in name else 0.7,
                label=name)
    
    ax.axvline(x=Nh, color='red', ls='--', alpha=0.4, label='Adversity onset')
    ax.axvline(x=t_rev, color='green', ls='--', alpha=0.4, label='Intervention')
    ax.set_xlabel('Trial')
    ax.set_ylabel('Latent P(engaged)')
    ax.set_title('(b) Example trajectories across conditions')
    ax.set_ylim(-0.05, 1.05)
    ax.legend(fontsize=7, loc='center left')
    
    fig.suptitle('Recovery is field-dominant', fontsize=12, y=1.03)
    plt.tight_layout()
    save_figure(fig, "fig_S2_clinical_prediction")
    plt.close()
    print("Saved fig_S2_clinical_prediction.[png|pdf]")
 
 
def run_all():
    """Entry point called by main.py. Produces both Step 2 figures."""
    print("Running Step 2: recovery boundary and field-dominance prediction...")
    fig_recovery_boundary()
    fig_clinical_prediction()


if __name__ == "__main__":
    print("=" * 60)
    print("Step 2: Hysteresis and the Recovery Boundary")
    print("=" * 60)
    print("\n1. Recovery boundary in (gamma, Dc) space...")
    fig_recovery_boundary()
    print("\n2. Clinical prediction: which variable matters?...")
    fig_clinical_prediction()
    print("\n" + "=" * 60)
    print(f"All figures in {OUT}")
    print("=" * 60)
 
