"""
Synthetic internal consistency check — bridging latent field/precision
dynamics to the empirical regression and classification designs.

This script does NOT alter the main simulation notebooks. It shows, in a
transparent synthetic-data setting, that:

  1. A field-dominant data-generating process (DGP) calibrated to the
     observed MCS coefficients (field β ≈ −0.30, precision β ≈ +0.02)
     recovers larger lagged field than precision coefficients in the same
     regression design used empirically — and that this holds across a
     range of reliability mismatches.

  2. A balance index can outperform distress-only and motivation-only
     models when the binary outcome is genuinely driven by a ratio-like
     phase coordinate.

Calibration note
----------------
The DGP parameters (field_strength, precision_strength) are set to match
the observed MCS primary-MI coefficients: field β ≈ −0.30, precision β
≈ +0.02 (n.s.).  The n.s. precision result is reproduced by setting
precision_strength close to zero (0.02).  The simulation therefore serves
as an internal consistency check — it confirms that the observed ratio of
coefficients is compatible with the field-dominant latent dynamics
postulated by the model — NOT as independent evidence for those dynamics.
"""

import csv
import os

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit


OUT = os.path.join(os.path.dirname(__file__), "Figs_psychopathology")
TABLE_OUT = os.path.join(os.path.dirname(__file__), "..", "outputs", "tables")
os.makedirs(OUT, exist_ok=True)
os.makedirs(TABLE_OUT, exist_ok=True)

plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 11,
    'axes.labelsize': 13,
    'axes.titlesize': 12,
    'figure.dpi': 150,
    'savefig.dpi': 150,
    'savefig.bbox': 'tight',
})

### DT ---> MCS observed coefficients — field β ≈ −0.305, precision β ≈ +0.016
### DT ---> These values are used to calibrate the DGP so the bridge is a
### DT ---> genuine consistency check rather than a self-confirming demonstration.
_MCS_FIELD_BETA = 0.305    # abs value; sign = protective
_MCS_PRECISION_BETA = 0.016  # abs value; n.s. in MCS


def zscore(x):
    x = np.asarray(x, dtype=float)
    sd = np.std(x)
    return (x - np.mean(x)) / (sd if sd > 1e-12 else 1.0)


def make_proxy(latent, reliability, rng):
    noise = rng.normal(size=len(latent))
    return np.sqrt(reliability) * zscore(latent) + np.sqrt(1 - reliability) * zscore(noise)


def ols_standardized(y, *xs):
    y = zscore(y)
    X = np.column_stack([np.ones(len(y))] + [zscore(x) for x in xs])
    beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    return beta[1:]


def auc_score(y_true, scores):
    y_true = np.asarray(y_true, dtype=int)
    scores = np.asarray(scores, dtype=float)
    order = np.argsort(scores)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(scores) + 1)
    pos = y_true == 1
    n_pos = np.sum(pos)
    n_neg = len(y_true) - n_pos
    if n_pos == 0 or n_neg == 0:
        return np.nan
    rank_sum = np.sum(ranks[pos])
    return (rank_sum - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)


def fit_logistic_1d(x, y):
    x = zscore(x)
    X = np.column_stack([np.ones(len(x)), x])
    y = np.asarray(y, dtype=float)

    def objective(beta):
        p = expit(X @ beta)
        p = np.clip(p, 1e-8, 1 - 1e-8)
        return -np.sum(y * np.log(p) + (1 - y) * np.log(1 - p))

    res = minimize(objective, np.zeros(X.shape[1]), method='BFGS')
    p = expit(X @ res.x)
    return {
        'coef': float(res.x[1]),
        'auc': float(auc_score(y, p)),
        'brier': float(np.mean((y - p) ** 2)),
    }


def simulate_lagged_dataset(n, field_reliability, precision_reliability, seed,
                            field_strength=_MCS_FIELD_BETA,
                            precision_strength=_MCS_PRECISION_BETA):
    """
    Simulate the lagged-regression design used empirically in MCS/ABCD.

    field_strength and precision_strength are calibrated to MCS observed
    coefficients by default.  The function remains parameterisable so the
    caller can confirm sensitivity to these choices.
    """
    rng = np.random.default_rng(seed)
    adversity = rng.normal(size=n)
    field_latent = 0.35 * adversity + rng.normal(scale=0.9, size=n)
    precision_latent = -0.10 * adversity + rng.normal(scale=0.9, size=n)

    baseline_distress = (
        0.75 * adversity
        - 0.42 * field_latent
        - 0.18 * precision_latent
        + rng.normal(scale=0.8, size=n)
    )
    followup_distress = (
        0.62 * baseline_distress
        + 0.28 * adversity
        - field_strength * field_latent
        - precision_strength * precision_latent
        + rng.normal(scale=0.8, size=n)
    )

    field_proxy = make_proxy(field_latent, field_reliability, rng)
    precision_proxy = make_proxy(precision_latent, precision_reliability, rng)
    betas = ols_standardized(followup_distress, baseline_distress, field_proxy, precision_proxy)
    return {
        'beta_baseline': float(betas[0]),
        'beta_field': float(betas[1]),
        'beta_precision': float(betas[2]),
        'beta_gap_abs': float(abs(betas[1]) - abs(betas[2])),
    }


def simulate_balance_dataset(n, seed, distress_rel=0.80, motivation_rel=0.80):
    rng = np.random.default_rng(seed)
    adversity = rng.normal(size=n)
    field_latent = 0.30 * adversity + rng.normal(scale=0.9, size=n)

    distress_true = np.exp(0.55 * adversity - 0.35 * field_latent + rng.normal(scale=0.45, size=n))
    motivation_true = np.exp(0.60 * field_latent - 0.20 * adversity + rng.normal(scale=0.45, size=n))

    distress_obs = np.exp(make_proxy(np.log(distress_true), distress_rel, rng) / 2)
    motivation_obs = np.exp(make_proxy(np.log(motivation_true), motivation_rel, rng) / 2)

    balance = np.log((distress_obs + 0.05) / (motivation_obs + 0.05))
    case_prob = expit(-1.15 + 1.85 * zscore(balance))
    case = rng.binomial(1, case_prob)

    return {
        'distress': fit_logistic_1d(np.log(distress_obs + 0.05), case),
        'motivation': fit_logistic_1d(np.log(motivation_obs + 0.05), case),
        'balance': fit_logistic_1d(balance, case),
    }


def run_lagged_bridge():
    ### DT ---> Scenarios bracket the MCS measurement properties.
    ### DT ---> "MCS-calibrated" uses matched reliabilities typical of multi-item scales.
    ### DT ---> The two mismatch scenarios show that field dominance widens with
    ### DT ---> poorer precision measurement but is not caused by it.
    scenarios = [
        ('Matched reliability\n(field=0.75, prec=0.75)', 0.75, 0.75),
        ('Precision attenuated\n(field=0.80, prec=0.45)', 0.80, 0.45),
        ('Precision noise high\n(field=0.80, prec=0.30)', 0.80, 0.30),
    ]
    rows = []
    for label, field_rel, precision_rel in scenarios:
        stats = []
        for rep in range(200):
            stats.append(simulate_lagged_dataset(
                n=5000,
                field_reliability=field_rel,
                precision_reliability=precision_rel,
                seed=1000 + rep + int(field_rel * 100) + int(precision_rel * 1000),
            ))
        beta_field = np.array([s['beta_field'] for s in stats])
        beta_precision = np.array([s['beta_precision'] for s in stats])
        beta_gap = np.array([s['beta_gap_abs'] for s in stats])
        rows.append({
            'scenario': label.replace('\n', ' '),
            'field_reliability': field_rel,
            'precision_reliability': precision_rel,
            'field_beta_mean': float(np.mean(beta_field)),
            'field_beta_sd': float(np.std(beta_field)),
            'precision_beta_mean': float(np.mean(beta_precision)),
            'precision_beta_sd': float(np.std(beta_precision)),
            'abs_gap_mean': float(np.mean(beta_gap)),
            'p_abs_field_gt_precision': float(np.mean(beta_gap > 0)),
        })
    return rows


def run_balance_bridge():
    reps = []
    for rep in range(150):
        reps.append(simulate_balance_dataset(n=5000, seed=9000 + rep))

    rows = []
    for model in ['distress', 'motivation', 'balance']:
        aucs = np.array([r[model]['auc'] for r in reps])
        briers = np.array([r[model]['brier'] for r in reps])
        rows.append({
            'model': model,
            'auc_mean': float(np.mean(aucs)),
            'auc_sd': float(np.std(aucs)),
            'brier_mean': float(np.mean(briers)),
            'brier_sd': float(np.std(briers)),
        })
    return rows


def save_rows(filename, rows):
    if not rows:
        return
    with open(os.path.join(TABLE_OUT, filename), 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def plot_bridge(lagged_rows, balance_rows):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6.5))

    # --- Panel (a): field vs precision coefficients across reliability scenarios ---
    ax = axes[0]
    x = np.arange(len(lagged_rows))
    field_means = [r['field_beta_mean'] for r in lagged_rows]
    field_sds = [r['field_beta_sd'] for r in lagged_rows]
    prec_means = [r['precision_beta_mean'] for r in lagged_rows]
    prec_sds = [r['precision_beta_sd'] for r in lagged_rows]
    width = 0.34

    ### DT ---> Determine y ceiling before drawing so annotations sit above bars
    y_tops = [max(abs(fm), abs(pm)) + sd + 0.05
              for fm, pm, sd in zip(field_means, prec_means, field_sds)]
    y_ceil = max(y_tops) + 0.14

    ax.bar(x - width / 2, field_means, width, color='#2E7D32', alpha=0.85,
           label='Field proxy (motivational)')
    ax.bar(x + width / 2, prec_means, width, color='#1565C0', alpha=0.85,
           label='Precision proxy')
    ax.errorbar(x - width / 2, field_means, yerr=field_sds,
                fmt='none', ecolor='black', capsize=5, lw=1.2)
    ax.errorbar(x + width / 2, prec_means, yerr=prec_sds,
                fmt='none', ecolor='black', capsize=5, lw=1.2)
    ax.axhline(0, color='gray', lw=0.8, ls='--')
    ax.set_xticks(x)
    ### DT ---> Rotate 20° so two-line labels don't overlap; anchor to right edge.
    ax.set_xticklabels([r['scenario'] for r in lagged_rows], fontsize=8.5,
                       rotation=20, ha='right', rotation_mode='anchor')
    ax.set_ylabel('Mean standardised lagged coefficient (±1 SD)')
    ax.set_ylim(None, y_ceil)
    ax.set_title('(a) Field dominance holds across reliability conditions\n'
                 '(DGP calibrated to MCS: field β≈−0.30, precision β≈+0.02)')
    ax.legend(fontsize=8.5, loc='lower right')
    for xi, row, ytop in zip(x, lagged_rows, y_tops):
        ax.text(xi, ytop,
                f"P(|field|>|prec|) = {row['p_abs_field_gt_precision']:.2f}",
                ha='center', va='bottom', fontsize=8)

    # --- Panel (b): balance index vs distress/motivation in classification ---
    ax = axes[1]
    auc_means = [r['auc_mean'] for r in balance_rows]
    auc_sds = [r['auc_sd'] for r in balance_rows]
    briers = [r['brier_mean'] for r in balance_rows]
    x2 = np.arange(len(balance_rows))
    bars = ax.bar(x2, auc_means, color=['#C62828', '#2E7D32', '#6A1B9A'], alpha=0.82)
    ax.errorbar(x2, auc_means, yerr=auc_sds, fmt='none', ecolor='black', capsize=5, lw=1.2)
    ax.set_xticks(x2)
    ax.set_xticklabels(['Distress only', 'Motivation only', 'Balance index'], fontsize=9)
    ax.set_ylabel('Mean AUC (±1 SD, 150 replicates)')
    ax.set_ylim(0.45, 0.95)
    ax.set_title('(b) Balance index predicts best\n'
                 'when caseness is driven by a ratio-like phase coordinate')
    for xi, auc, sd, brier in zip(x2, auc_means, auc_sds, briers):
        ax.text(xi, auc + sd + 0.025, f'Brier={brier:.3f}',
                ha='center', va='bottom', fontsize=8)

    fig.suptitle(
        'Synthetic internal consistency check: field-dominant DGP\n'
        'reproduces the empirical coefficient ordering and balance-index advantage\n'
        '(calibrated to observed MCS coefficients — not independent evidence)',
        fontsize=11.5, y=1.03
    )
    ### DT ---> pad=1.5 gives headroom for the three-line x-tick labels;
    ### DT ---> savefig.bbox=tight in rcParams expands the saved canvas automatically.
    plt.tight_layout(pad=1.5)
    plt.savefig(f"{OUT}/fig_S_bridge_empirical.png")
    plt.close()
    print("Saved fig_S_bridge_empirical.png")


def run_all():
    """Entry point called by main.py."""
    print("Running empirical bridge (synthetic consistency check)...")
    lagged_rows = run_lagged_bridge()
    balance_rows = run_balance_bridge()
    save_rows("simulation_bridge_lagged.csv", lagged_rows)
    save_rows("simulation_bridge_balance.csv", balance_rows)
    plot_bridge(lagged_rows, balance_rows)
    print("  Lagged bridge summary:")
    for row in lagged_rows:
        print(f"    {row['scenario'][:40]}: "
              f"field={row['field_beta_mean']:.3f}, "
              f"prec={row['precision_beta_mean']:.3f}, "
              f"P(|f|>|p|)={row['p_abs_field_gt_precision']:.2f}")
    print("  Balance bridge summary:")
    for row in balance_rows:
        print(f"    {row['model']}: AUC={row['auc_mean']:.3f} ± {row['auc_sd']:.3f}")


if __name__ == "__main__":
    print("=" * 60)
    print("Standalone empirical bridge (synthetic consistency check)")
    print("=" * 60)
    run_all()
    print("=" * 60)
