"""
Core mathematical functions for Phase Transitions in Active Inference II
========================================================================
Reconstructed from the verified sim_bifurcation.py of original paper/code:
see https://github.com/dtsomoucl/phase-transitions-in-active-inference
These functions are the foundation on which this extension is built.
"""

import numpy as np
from scipy.optimize import brentq
 
 
# ── Entropy helpers ──────────────────────────────────────────────────
 
def binary_entropy(x):
    """Binary entropy H(x) = -x ln(x) - (1-x) ln(1-x)."""
    x = np.clip(x, 1e-12, 1 - 1e-12)
    return -x * np.log(x) - (1 - x) * np.log(1 - x)
 
 
def binary_entropy_deriv(x):
    """H'(x) = ln((1-x)/x)."""
    x = np.clip(x, 1e-12, 1 - 1e-12)
    return np.log((1 - x) / x)
 
 
# ── Analytical functions ─────────────────────────────────────────────
 
def coupling_function(tau, p):
    """
    Coupling function G(tau; p) from Eq. (9) of the paper.
 
    Parameters
    ----------
    tau : float   – rescaled developmental time N/(2*alpha_0)
    p   : float   – true discriminability in (0.5, 1)
 
    Returns
    -------
    G   : float   – coupling strength at this point in development
    """
    if tau < 1e-12:
        return 0.0
    a_bar = 0.5 + (p - 0.5) * tau / (1 + tau)
    lam = (p - 0.5) * tau / (1 + tau) ** 2
    log_ratio = np.abs(np.log((1 - a_bar) / a_bar))
    return lam * log_ratio
 
 
def find_gamma_c(p, tau_range=np.linspace(0.01, 50, 5000)):
    """Find critical precision gamma_c = 1 / max_tau G(tau; p)."""
    G_vals = np.array([coupling_function(t, p) for t in tau_range])
    G_max = np.max(G_vals)
    tau_max = tau_range[np.argmax(G_vals)]
    if G_max < 1e-12:
        return np.inf, 0, 0
    return 1.0 / G_max, tau_max, G_max
 
 
# ── EFE and policy selection ────────────────────────────────────────
 
def compute_efe_difference(a, b, delta_c):
    """
    ΔG = G(π₁) − G(π₂) = (1−a−b)Δc + H(a) − H(b).
    Eq. (5) of the paper.
    """
    risk_diff = (1 - a - b) * delta_c
    ambiguity_diff = binary_entropy(a) - binary_entropy(b)
    return risk_diff + ambiguity_diff
 
 
def policy_posterior_pi1(delta_G, gamma):
    """P(π₁) = σ(−γ ΔG) = 1 / (1 + exp(γ ΔG)).  Eq. (4)."""
    x = gamma * delta_G
    if x > 500:
        return 0.0
    elif x < -500:
        return 1.0
    else:
        return 1.0 / (1.0 + np.exp(x))
 
 
def self_consistency_rhs(z, tau, p, gamma, delta_c):
    """
    RHS of the self-consistency equation (8):
        tanh(−γ ΔG(z) / 2)
    where a(z) and b(z) depend on z and tau via Eq. (7).
    """
    n1 = (1 + z) * tau
    n2 = (1 - z) * tau
    a = 0.5 + (p - 0.5) * n1 / (1 + n1) if n1 > 0 else 0.5
    b = 0.5 + (p - 0.5) * n2 / (1 + n2) if n2 > 0 else 0.5
    dG = compute_efe_difference(a, b, delta_c)
    return np.tanh(-gamma / 2 * dG)
 
 
# ── Fixed-point finder ──────────────────────────────────────────────
 
def find_fixed_points(tau, p, gamma, delta_c,
                      z_grid=np.linspace(-0.999, 0.999, 2000)):
    """
    Find fixed points of z = tanh(−γ ΔG(z)/2) by scanning + Brent.
    Returns list of (z_fp, stable) tuples.
    """
    f_vals = np.array([
        self_consistency_rhs(z, tau, p, gamma, delta_c) - z
        for z in z_grid
    ])
    fps = []
    for i in range(len(f_vals) - 1):
        if f_vals[i] * f_vals[i + 1] < 0:
            try:
                z_fp = brentq(
                    lambda z: self_consistency_rhs(z, tau, p, gamma, delta_c) - z,
                    z_grid[i], z_grid[i + 1],
                )
                eps = 1e-5
                rhs_plus = self_consistency_rhs(z_fp + eps, tau, p, gamma, delta_c)
                rhs_minus = self_consistency_rhs(z_fp - eps, tau, p, gamma, delta_c)
                deriv = (rhs_plus - rhs_minus) / (2 * eps)
                stable = abs(deriv) < 1
                fps.append((z_fp, stable))
            except Exception:
                pass
    # Deduplicate
    fps.sort(key=lambda item: item[0])
    unique = []
    for z_fp, stable in fps:
        if not unique or abs(z_fp - unique[-1][0]) > 1e-4:
            unique.append((z_fp, stable))
    return unique
 
 
# ── Single-agent simulation (original) ──────────────────────────────
 
def run_single_agent(p, alpha_0, gamma, delta_c, N_total, seed=None):
    """
    Run one agent with fixed parameters.  Faithful reconstruction of
    the original sim_bifurcation.run_single_agent().
    """
    rng = np.random.default_rng(seed)
 
    alpha1, beta1 = alpha_0 / 2, alpha_0 / 2
    alpha2, beta2 = alpha_0 / 2, alpha_0 / 2
    A_true = np.array([[p, 1 - p], [1 - p, p]])
 
    history = {
        'a': np.zeros(N_total), 'b': np.zeros(N_total),
        'phi': np.zeros(N_total), 'z': np.zeros(N_total),
        'delta_G': np.zeros(N_total), 'P_pi1': np.zeros(N_total),
        'action': np.zeros(N_total, dtype=int),
        'n1': np.zeros(N_total), 'n2': np.zeros(N_total),
    }
    n1_count = 0
    n2_count = 0
 
    for t in range(N_total):
        a = alpha1 / (alpha1 + beta1)
        b = beta2 / (alpha2 + beta2)
        dG = compute_efe_difference(a, b, delta_c)
        P_pi1 = policy_posterior_pi1(dG, gamma)
        action = 0 if rng.random() < P_pi1 else 1
 
        if action == 0:
            obs = 0 if rng.random() < A_true[0, 0] else 1
            if obs == 0:
                alpha1 += 1
            else:
                beta1 += 1
            n1_count += 1
        else:
            obs = 0 if rng.random() < A_true[0, 1] else 1
            if obs == 0:
                alpha2 += 1
            else:
                beta2 += 1
            n2_count += 1
 
        N_so_far = n1_count + n2_count
        history['a'][t] = a
        history['b'][t] = b
        history['phi'][t] = n1_count / max(N_so_far, 1)
        history['z'][t] = 2 * history['phi'][t] - 1
        history['delta_G'][t] = dG
        history['P_pi1'][t] = P_pi1
        history['action'][t] = action
        history['n1'][t] = n1_count
        history['n2'][t] = n2_count
 
    return history
