"""Core expected-free-energy, coupling, and fixed-point utilities."""

import numpy as np
from scipy.optimize import brentq
 
 
 
def binary_entropy(x):
    """Binary entropy H(x) = -x ln(x) - (1-x) ln(1-x)."""
    x = np.clip(x, 1e-12, 1 - 1e-12)
    return -x * np.log(x) - (1 - x) * np.log(1 - x)
 
 
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
                stable = deriv < 1
                fps.append((z_fp, stable))
            except Exception:
                pass
    fps.sort(key=lambda item: item[0])
    unique = []
    for z_fp, stable in fps:
        if not unique or abs(z_fp - unique[-1][0]) > 1e-4:
            unique.append((z_fp, stable))
    return unique
