"""
Step 3 (exploratory): Orthogonal Interventions When Preferences
Cannot Be Fully Restored
==========================================================================
Clinical motivation: one cannot simply prescribe motivation.  If Dc can
only be PARTIALLY restored, do count tempering or environmental restructuring
(increasing p) close the remaining gap?
 
Findings after adding Step-2-style robustness checks:
  - Orthogonal interventions have LIMITED additional leverage because
    recovery is field-dominated in this minimal model.
  - Count tempering helps when Dc is already favourable; it can hurt
    when Dc is adverse.
  - p-increase amplifies whichever direction the field favours:
    helpful when Dc > 0, counterproductive when Dc < 0.
  - No orthogonal intervention substitutes for restoring Dc.
 
Status: internally consistent with the field-dominance result from
Step 2, now supplemented with threshold sensitivity and regime
dependence checks to show that any orthogonal benefit is secondary.
"""

import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, sys
from core_functions import compute_efe_difference, policy_posterior_pi1
from psychopathology_regimes import recovery_regime

OUT = os.path.join(os.path.dirname(__file__), "Figs_psychopathology")
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({'font.family':'serif','font.size':11,'axes.labelsize':13,
    'axes.titlesize':12,'figure.dpi':150,'savefig.dpi':150,'savefig.bbox':'tight'})

REGIME = recovery_regime()


def run_recovery_agent(p, a0, seed, gh, dch, Nh, gr, dcr, gf, dcf, T_ill,
                       rg, rdc, p_new=None, eta=None, N_recovery=400):
    """
    Run onset then recovery with optional orthogonal interventions.
    eta: count tempering factor (0<eta<1 rescales Dirichlet toward prior)
    p_new: new environmental discriminability during recovery
 
    Returns: mean P(pi1) over last 100 recovery trials (unified with Step 2).
    """
    rng = np.random.default_rng(seed)
    al1=be1=al2=be2=a0/2
    At = np.array([[p,1-p],[1-p,p]])
    t_rev = Nh + T_ill
    N_total = t_rev + N_recovery
    Pp_late = []
 
    for t in range(N_total):
        if t < Nh: gt=gh; dct=dch
        else:
            dt = t - Nh; gt=max(gh-gr*dt, gf); dct=max(dch-dcr*dt, dcf)
 
        if t == t_rev:
            # Apply interventions at reversal point
            if eta is not None and eta < 1.0:
                base = a0/2
                al1 = eta*al1 + (1-eta)*base
                be1 = eta*be1 + (1-eta)*base
                al2 = eta*al2 + (1-eta)*base
                be2 = eta*be2 + (1-eta)*base
            if p_new is not None:
                At = np.array([[p_new,1-p_new],[1-p_new,p_new]])
 
        if t >= t_rev:
            gt = rg; dct = rdc
 
        a = al1/(al1+be1); b = be2/(al2+be2)
        dG = compute_efe_difference(a, b, dct)
        Pp = policy_posterior_pi1(dG, gt)
        act = 0 if rng.random() < Pp else 1
 
        if act == 0:
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
 
 
def run_condition(na, seed_base, **kw):
    vals = []
    for i in range(na):
        vals.append(run_recovery_agent(seed=seed_base+i, **kw))
    return np.array(vals)
 
 
# ====================================================================
# FIGURE S3a: Count tempering x partial Dc restoration
# ====================================================================
def fig_tempering_x_dc():
    """
    Heatmap: median recovery as a function of (eta, Dc_restored).
    gamma is fully restored.  Shows where tempering adds value.
    """
    p=REGIME["p"]; a0=REGIME["alpha_0"]; Nh=REGIME["N_healthy"]; T_ill=200; na=60
    gh=REGIME["gamma_healthy"]; dch=REGIME["delta_c_healthy"]
    gr=REGIME["gamma_rate"]; dcr=REGIME["delta_c_rate"]
    gf=REGIME["gamma_floor"]; dcf=REGIME["delta_c_floor"]
 
    eta_vals = [1.0, 0.5, 0.3, 0.1, 0.01]
    dc_vals = np.linspace(-0.10, 0.15, 11)
 
    recovery_grid = np.zeros((len(eta_vals), len(dc_vals)))
 
    for ei, eta in enumerate(eta_vals):
        for di, dc in enumerate(dc_vals):
            vals = run_condition(na, seed_base=50000+ei*1000+di*100,
                p=p, a0=a0, gh=gh, dch=dch, Nh=Nh, gr=gr, dcr=dcr,
                gf=gf, dcf=dcf, T_ill=T_ill,
                rg=gh, rdc=dc, eta=eta if eta<1.0 else None)
            recovery_grid[ei, di] = np.median(vals)
 
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
 
    # Panel (a): heatmap
    ax = axes[0]
    im = ax.imshow(recovery_grid, aspect='auto', cmap='RdYlGn',
                   vmin=0, vmax=1, origin='lower')
    ax.set_xticks(range(len(dc_vals)))
    ax.set_xticklabels([f'{v:.2f}' for v in dc_vals], rotation=45, fontsize=8)
    ax.set_yticks(range(len(eta_vals)))
    ax.set_yticklabels([f'{v}' for v in eta_vals])
    ax.set_xlabel(r'Restored $\Delta c$')
    ax.set_ylabel(r'Count tempering $\eta$ (1=none, 0.01=near-reset)')
    for ei in range(len(eta_vals)):
        for di in range(len(dc_vals)):
            v = recovery_grid[ei, di]
            c = 'white' if v < 0.4 or v > 0.7 else 'black'
            ax.text(di, ei, f'{v:.2f}', ha='center', va='center', fontsize=7, color=c)
    cb = plt.colorbar(im, ax=ax); cb.set_label('Median late P(engaged)')
    ax.set_title(r'(a) Recovery: count tempering $\times$ restored $\Delta c$'
                 f'\n(T_ill={T_ill}, ' r'$\gamma$' f' restored to {gh})')
 
    # Panel (b): marginal effect of tempering at selected Dc values
    ax = axes[1]
    dc_select = [-0.05, 0.0, 0.05, 0.10]
    colours = ['red', 'orange', 'steelblue', 'darkgreen']
    for dc, col in zip(dc_select, colours):
        di = np.argmin(np.abs(dc_vals - dc))
        ax.plot(eta_vals, recovery_grid[:, di], 'o-', color=col, ms=6, lw=1.5,
                label=rf'$\Delta c_{{res}}={dc:.2f}$')
    ax.set_xlabel(r'Count tempering $\eta$')
    ax.set_ylabel('Median late P(engaged)')
    ax.axhline(y=0.5, color='gray', ls='--', alpha=0.3)
    ax.set_title('(b) Marginal effect of count tempering\nat different ' r'$\Delta c$' ' levels')
    ax.legend(fontsize=8)
    ax.invert_xaxis()  # eta=1 (no tempering) on right, eta=0.01 (strong) on left
 
    fig.suptitle('Step 3: Count Tempering as Orthogonal Intervention\n'
                 r'($\gamma$ fully restored; varying $\eta$ and $\Delta c_{res}$)',
                 fontsize=13, y=1.04)
    plt.tight_layout()
    plt.savefig(f"{OUT}/fig_S3_tempering.png")
    plt.close()
    print("Saved fig_S3_tempering.png")
 
 
# ====================================================================
# FIGURE S3b: Environmental restructuring (p-increase)
# ====================================================================
def fig_p_increase():
    """
    Effect of increasing environmental discriminability during recovery.
    """
    p=REGIME["p"]; a0=REGIME["alpha_0"]; Nh=REGIME["N_healthy"]; T_ill=200; na=60
    gh=REGIME["gamma_healthy"]; dch=REGIME["delta_c_healthy"]
    gr=REGIME["gamma_rate"]; dcr=REGIME["delta_c_rate"]
    gf=REGIME["gamma_floor"]; dcf=REGIME["delta_c_floor"]
 
    p_new_vals = [0.75, 0.85, 0.90, 0.95, 0.99]
    dc_vals = np.linspace(-0.10, 0.15, 11)
 
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for pn, col in zip(p_new_vals, plt.cm.viridis(np.linspace(0.2, 0.9, len(p_new_vals)))):
        meds = []
        for dc in dc_vals:
            vals = run_condition(na, seed_base=60000+int(pn*100)*100+int((dc+0.1)*100),
                p=p, a0=a0, gh=gh, dch=dch, Nh=Nh, gr=gr, dcr=dcr,
                gf=gf, dcf=dcf, T_ill=T_ill,
                rg=gh, rdc=dc, p_new=pn)
            meds.append(np.median(vals))
        ax.plot(dc_vals, meds, 'o-', color=col, ms=5, lw=1.5,
                label=rf'$p_{{new}}={pn}$')
 
    ax.axhline(y=0.5, color='gray', ls='--', alpha=0.3)
    ax.set_xlabel(r'Restored $\Delta c$')
    ax.set_ylabel('Median late P(engaged)')
    ax.set_title('Effect of Environmental Restructuring (increasing p)\n'
                 r'($\gamma$ fully restored, T_ill=200)')
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(f"{OUT}/fig_S3_p_increase.png")
    plt.close()
    print("Saved fig_S3_p_increase.png")
 
 
# ====================================================================
# FIGURE S3c: Combined intervention comparison
# ====================================================================
def fig_combined():
    """
    Bar chart: recovery under realistic clinical scenarios where Dc
    is only PARTIALLY restored (to 0.0 rather than 0.15).
    
    Conditions:
      1. Partial Dc only (Dc→0, gamma stays low)
      2. Partial Dc + gamma restored
      3. Partial Dc + gamma + count tempering (eta=0.1)
      4. Partial Dc + gamma + p-increase (p→0.95)
      5. Partial Dc + gamma + both orthogonal
      6. Full restoration (benchmark)
    """
    p=REGIME["p"]; a0=REGIME["alpha_0"]; Nh=REGIME["N_healthy"]; T_ill=200; na=100
    gh=REGIME["gamma_healthy"]; dch=REGIME["delta_c_healthy"]
    gr=REGIME["gamma_rate"]; dcr=REGIME["delta_c_rate"]
    gf=REGIME["gamma_floor"]; dcf=REGIME["delta_c_floor"]
    dc_partial = 0.0  # Dc restored to neutral, not positive
 
    conditions = [
        (r'Partial $\Delta c$ only',        dict(rg=gf, rdc=dc_partial, eta=None, p_new=None)),
        (r'+ restore $\gamma$',             dict(rg=gh, rdc=dc_partial, eta=None, p_new=None)),
        (r'+ $\gamma$ + temper ($\eta$=0.1)', dict(rg=gh, rdc=dc_partial, eta=0.1, p_new=None)),
        (r'+ $\gamma$ + p=0.95',            dict(rg=gh, rdc=dc_partial, eta=None, p_new=0.95)),
        (r'+ $\gamma$ + temper + p=0.95',   dict(rg=gh, rdc=dc_partial, eta=0.1, p_new=0.95)),
        ('Full restoration\n(benchmark)',    dict(rg=gh, rdc=dch, eta=None, p_new=None)),
    ]
    colours = ['#E8E8E8', '#90CAF9', '#64B5F6', '#FFA726', '#66BB6A', '#2E7D32']
 
    fig, ax = plt.subplots(figsize=(12, 6))
    for ci, (name, kw) in enumerate(conditions):
        vals = run_condition(na, seed_base=70000+ci*1000,
            p=p, a0=a0, gh=gh, dch=dch, Nh=Nh, gr=gr, dcr=dcr,
            gf=gf, dcf=dcf, T_ill=T_ill, **kw)
        med = np.median(vals)
        q1, q3 = np.percentile(vals, [25, 75])
        ax.bar(ci, med, color=colours[ci], alpha=0.85, edgecolor='black', lw=0.8)
        ax.errorbar(ci, med, yerr=[[med-q1], [q3-med]], fmt='none',
                    capsize=6, color='black', lw=1.5)
 
    ax.set_xticks(range(len(conditions)))
    ax.set_xticklabels([c[0] for c in conditions], fontsize=9, ha='center')
    ax.set_ylabel('Median late P(engaged)')
    ax.axhline(y=0.5, color='gray', ls='--', alpha=0.5, label='Recovery threshold')
    ax.set_ylim(0, 1.1)
    ax.set_title(
        rf'Orthogonal additions help only after $\Delta c$ is at least partially restored'
        '\n'
        rf'($\Delta c_{{res}}$={dc_partial} [neutral]; T_ill={T_ill}; γ and $\Delta c$ floors from RECOVERY_REGIME)',
        fontsize=11
    )
    ax.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(f"{OUT}/fig_S3_combined.png")
    plt.close()
    print("Saved fig_S3_combined.png")


# ====================================================================
# FIGURE S3d: Threshold sensitivity
# ====================================================================
def fig_threshold_sensitivity():
    p=REGIME["p"]; a0=REGIME["alpha_0"]; Nh=REGIME["N_healthy"]; T_ill=200; na=80
    gh=REGIME["gamma_healthy"]; dch=REGIME["delta_c_healthy"]
    gr=REGIME["gamma_rate"]; dcr=REGIME["delta_c_rate"]
    gf=REGIME["gamma_floor"]; dcf=REGIME["delta_c_floor"]
    dc_partial = 0.0

    conditions = [
        ('Partial\nDc only', dict(rg=gf, rdc=dc_partial, eta=None, p_new=None)),
        ('+ gamma', dict(rg=gh, rdc=dc_partial, eta=None, p_new=None)),
        ('+ temper', dict(rg=gh, rdc=dc_partial, eta=0.1, p_new=None)),
        ('+ p', dict(rg=gh, rdc=dc_partial, eta=None, p_new=0.95)),
        ('+ both', dict(rg=gh, rdc=dc_partial, eta=0.1, p_new=0.95)),
        ('Full', dict(rg=gh, rdc=dch, eta=None, p_new=None)),
    ]
    colours = ['#E8E8E8', '#90CAF9', '#64B5F6', '#FFA726', '#66BB6A', '#2E7D32']
    thresholds = [0.5, 0.6, 0.7]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    for ti, th in enumerate(thresholds):
        ax = axes[ti]
        for ci, (_, kw) in enumerate(conditions):
            vals = run_condition(
                na, seed_base=81000 + ti * 10000 + ci * 1000,
                p=p, a0=a0, gh=gh, dch=dch, Nh=Nh, gr=gr, dcr=dcr,
                gf=gf, dcf=dcf, T_ill=T_ill, **kw
            )
            frac_recovered = np.mean(vals > th)
            ax.bar(ci, frac_recovered, color=colours[ci], alpha=0.8, edgecolor='black')
            ax.text(ci, frac_recovered + 0.02, f'{frac_recovered:.0%}',
                    ha='center', fontsize=8)
        ax.set_xticks(range(len(conditions)))
        ax.set_xticklabels([name for name, _ in conditions], fontsize=8)
        ax.set_ylabel('Fraction recovered')
        ax.set_ylim(0, 1.12)
        ax.set_title(f'Threshold = {th}')

    fig.suptitle('Step 3 robustness: orthogonal-intervention ranking across thresholds\n'
                 '(orthogonal additions help only after partial field restoration)',
                 fontsize=13, y=1.03)
    plt.tight_layout()
    plt.savefig(f"{OUT}/fig_S3_threshold_sensitivity.png")
    plt.close()
    print("Saved fig_S3_threshold_sensitivity.png")


# ====================================================================
# FIGURE S3e: Regime dependence
# ====================================================================
def fig_regime_dependence():
    """
    Heatmap of the best orthogonal gain over the gamma+partial-Dc baseline
    across (p, alpha_0). Positive values indicate that orthogonal additions
    help, but the magnitude shows whether that help is substantive.
    """
    Nh = REGIME["N_healthy"]; T_ill = 200; na = 40
    gh = REGIME["gamma_healthy"]; dch = REGIME["delta_c_healthy"]
    gr = REGIME["gamma_rate"]; gf = REGIME["gamma_floor"]
    p_vals = np.array([0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95])
    a0_vals = np.array([1.0, 2.0, 4.0, 8.0])
    dc_partial = 0.0

    gain = np.zeros((len(a0_vals), len(p_vals)))
    for ai, a0 in enumerate(a0_vals):
        for pi, p in enumerate(p_vals):
            dcr = dch / 200
            dcf = -dch
            p_new = min(0.95, max(p, 0.90))

            base = np.median(run_condition(
                na, seed_base=90000 + ai * 1000 + pi * 100,
                p=p, a0=a0, gh=gh, dch=dch, Nh=Nh, gr=gr, dcr=dcr,
                gf=gf, dcf=dcf, T_ill=T_ill,
                rg=gh, rdc=dc_partial, eta=None, p_new=None
            ))
            temper = np.median(run_condition(
                na, seed_base=91000 + ai * 1000 + pi * 100,
                p=p, a0=a0, gh=gh, dch=dch, Nh=Nh, gr=gr, dcr=dcr,
                gf=gf, dcf=dcf, T_ill=T_ill,
                rg=gh, rdc=dc_partial, eta=0.1, p_new=None
            ))
            restructure = np.median(run_condition(
                na, seed_base=92000 + ai * 1000 + pi * 100,
                p=p, a0=a0, gh=gh, dch=dch, Nh=Nh, gr=gr, dcr=dcr,
                gf=gf, dcf=dcf, T_ill=T_ill,
                rg=gh, rdc=dc_partial, eta=None, p_new=p_new
            ))
            combo = np.median(run_condition(
                na, seed_base=93000 + ai * 1000 + pi * 100,
                p=p, a0=a0, gh=gh, dch=dch, Nh=Nh, gr=gr, dcr=dcr,
                gf=gf, dcf=dcf, T_ill=T_ill,
                rg=gh, rdc=dc_partial, eta=0.1, p_new=p_new
            ))
            gain[ai, pi] = max(temper, restructure, combo) - base

    fig, ax = plt.subplots(figsize=(8, 5))
    im = ax.imshow(gain, aspect='auto', cmap='YlGnBu', vmin=-0.05, vmax=0.35, origin='lower')
    ax.set_xticks(range(len(p_vals)))
    ax.set_xticklabels([f'{v:.2f}' for v in p_vals])
    ax.set_yticks(range(len(a0_vals)))
    ax.set_yticklabels([f'{v:.1f}' for v in a0_vals])
    ax.set_xlabel('Environmental discriminability p')
    ax.set_ylabel(r'Prior concentration $\alpha_0$')
    for ai in range(len(a0_vals)):
        for pi in range(len(p_vals)):
            val = gain[ai, pi]
            colour = 'white' if val > 0.18 else 'black'
            ax.text(pi, ai, f'{val:+.2f}', ha='center', va='center',
                    fontsize=9, color=colour, fontweight='bold')
    cb = plt.colorbar(im, ax=ax)
    cb.set_label('Best orthogonal lift over gamma + partial Dc baseline')
    ax.set_title('Step 3 robustness: orthogonal gains are secondary across regimes\n'
                 '(positive values are usually modest and do not replace field restoration)',
                 fontsize=12)
    plt.tight_layout()
    plt.savefig(f"{OUT}/fig_S3_regime_dependence.png")
    plt.close()
    print("Saved fig_S3_regime_dependence.png")
 
 
def run_all():
    """
    Entry point called by main.py.

    Produces only fig_S3_combined.png — the combined intervention bar chart
    that directly backs up the field-dominance claim (partial Dc restoration
    without gamma or orthogonal additions leaves agents withdrawn; restoring
    both Dc and gamma is the key lever).

    The supplementary figures (tempering heatmap, p-increase, threshold
    sensitivity, regime dependence) remain callable individually but are not
    part of the main publication sweep.
    """
    print("Running Step 3: orthogonal interventions (combined figure)...")
    fig_combined()


if __name__ == "__main__":
    print("=" * 60)
    print("Step 3: Orthogonal Interventions")
    print("=" * 60)
    print("\n(a) Count tempering x Dc restoration...")
    fig_tempering_x_dc()
    print("\n(b) Environmental restructuring (p-increase)...")
    fig_p_increase()
    print("\n(c) Combined intervention comparison...")
    fig_combined()
    print("\n(d) Threshold sensitivity (supplementary, not in main sweep)...")
    fig_threshold_sensitivity()
    print("\n(e) Regime dependence (supplementary, not in main sweep)...")
    fig_regime_dependence()
    print("\n" + "=" * 60)
 
