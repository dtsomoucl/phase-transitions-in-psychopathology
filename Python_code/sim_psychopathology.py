"""
Pathological Phase Transitions in Active Inference - Step 1
===========================================================
Extension of original work to model onset of psychopathology
(e.g. adolescent depression) as a *reverse* bifurcation, and
recovery via path-asymmetric intervention.

Step 1: Onset dynamics with drifting gamma and delta_c.

The agent undergoes three life-phases:
  Phase A  --  healthy development (fixed gamma_0, delta_c_0 > 0)
  Phase B  --  adversity / prodromal (gamma down, delta_c down slowly)
  Phase C  --  illness consolidation (gamma and delta_c stabilise at low values)

We verify that the reverse transition reproduces the catastrophe-flag
signatures from the original paper, now in the clinical direction
(engaged to withdrawn).

Changes in this version
-----------------------
  1. fig_P2 panel (b) — split into three separate time-series sub-panels
     (P(engaged), γ(t), Δc(t)); the ×30 scaling artefact is removed.
  2. fig_P_flags — bimodality panel uses a dynamically-chosen time point
     near the variance peak; divergence panel shows individual paired
     trajectory bifurcations rather than a dose-response curve; the
     redundant variance panel is removed (it is shown in fig_P2).
  3. fig_P3 EWS — autocorrelation panel replaced by event-aligned AC:
     each agent's trace is aligned at that agent's own transition time
     so the pre-transition rise (critical slowing down) is not masked by
     population heterogeneity in jump timing.
  4. fig_P_reversal_true is retained as a callable function but is not
     produced by run_all() / main.py.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.ndimage import gaussian_filter1d
import os, sys

sys.path.insert(0, os.path.dirname(__file__))
from core_functions import (
    binary_entropy, compute_efe_difference, policy_posterior_pi1,
    coupling_function, find_gamma_c, find_fixed_points,
    self_consistency_rhs,
)
from psychopathology_regimes import onset_regime

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

DEFAULTS = onset_regime()


def save_figure(fig, stem):
    fig.savefig(os.path.join(OUT, f"{stem}.png"))
    fig.savefig(os.path.join(OUT, f"{stem}.pdf"))


# ====================================================================
# Causal moving average (no edge artefacts)
# ====================================================================
def causal_ma(arr, w=20):
    out = np.full(len(arr), np.nan, dtype=float)
    cs = np.cumsum(arr)
    for t in range(len(arr)):
        lo = max(0, t - w + 1)
        n = t - lo + 1
        out[t] = (cs[t] - (cs[lo - 1] if lo > 0 else 0.0)) / n
    return out


# ====================================================================
# Agent simulator
# ====================================================================
def run_agent_with_drift(
    p, alpha_0, N_total, seed=None,
    gamma_healthy=16.0, delta_c_healthy=0.3, N_healthy=200,
    gamma_rate=0.0, delta_c_rate=0.0,
    gamma_floor=4.0, delta_c_floor=-0.5,
    t_reversal=None, gamma_reversal=None, delta_c_reversal=None,
    t_temper=None, temper_eta=1.0,
):
    rng = np.random.default_rng(seed)
    al1, be1 = alpha_0 / 2, alpha_0 / 2
    al2, be2 = alpha_0 / 2, alpha_0 / 2
    A_true = np.array([[p, 1 - p], [1 - p, p]])

    keys = ["a", "b", "phi", "z", "delta_G", "P_pi1",
            "gamma_t", "delta_c_t", "al1", "be1", "al2", "be2", "n1", "n2"]
    hist = {k: np.zeros(N_total) for k in keys}
    hist["action"] = np.zeros(N_total, dtype=int)
    n1 = n2 = 0

    for t in range(N_total):
        if t < N_healthy:
            g = gamma_healthy; dc = delta_c_healthy
        else:
            dt = t - N_healthy
            g = max(gamma_healthy - gamma_rate * dt, gamma_floor)
            dc = max(delta_c_healthy - delta_c_rate * dt, delta_c_floor)

        if t_reversal is not None and t >= t_reversal:
            if gamma_reversal is not None: g = gamma_reversal
            if delta_c_reversal is not None: dc = delta_c_reversal

        if t_temper is not None and t == t_temper:
            base = alpha_0 / 2
            al1 = temper_eta * al1 + (1 - temper_eta) * base
            be1 = temper_eta * be1 + (1 - temper_eta) * base
            al2 = temper_eta * al2 + (1 - temper_eta) * base
            be2 = temper_eta * be2 + (1 - temper_eta) * base

        a = al1 / (al1 + be1)
        b = be2 / (al2 + be2)
        dG = compute_efe_difference(a, b, dc)
        Pp1 = policy_posterior_pi1(dG, g)
        act = 0 if rng.random() < Pp1 else 1

        if act == 0:
            obs = 0 if rng.random() < A_true[0, 0] else 1
            if obs == 0: al1 += 1
            else: be1 += 1
            n1 += 1
        else:
            obs = 0 if rng.random() < A_true[0, 1] else 1
            if obs == 0: al2 += 1
            else: be2 += 1
            n2 += 1

        Nsf = n1 + n2
        hist["a"][t] = a; hist["b"][t] = b
        hist["phi"][t] = n1 / max(Nsf, 1)
        hist["z"][t] = 2 * hist["phi"][t] - 1
        hist["delta_G"][t] = dG; hist["P_pi1"][t] = Pp1
        hist["action"][t] = act
        hist["n1"][t] = n1; hist["n2"][t] = n2
        hist["gamma_t"][t] = g; hist["delta_c_t"][t] = dc
        hist["al1"][t] = al1; hist["be1"][t] = be1
        hist["al2"][t] = al2; hist["be2"][t] = be2
    return hist


# ====================================================================
# Quasi-static fixed-point analysis
# ====================================================================
def quasistatic_fps(hist, p, alpha_0, every=5):
    N = len(hist["gamma_t"])
    ts = np.arange(0, N, every)
    out = dict(t=ts, n_fps=np.zeros(len(ts), dtype=int),
               z_eng=np.full(len(ts), np.nan),
               z_wdr=np.full(len(ts), np.nan),
               z_uns=np.full(len(ts), np.nan))
    for idx, t in enumerate(ts):
        g = hist["gamma_t"][t]
        dc = hist["delta_c_t"][t]
        ntot = hist["n1"][t] + hist["n2"][t]
        tau = ntot / (2 * alpha_0) if alpha_0 > 0 else 0
        fps = find_fixed_points(tau, p, g, dc,
                                z_grid=np.linspace(-0.999, 0.999, 800))
        out["n_fps"][idx] = len(fps)
        stab = sorted([z for z, s in fps if s])
        unst = [z for z, s in fps if not s]
        if stab:
            out["z_eng"][idx] = stab[-1]
            if len(stab) > 1: out["z_wdr"][idx] = stab[0]
        if unst: out["z_uns"][idx] = unst[0]
    return out


def detect_fold(qs):
    ts = qs["t"]
    n_fps = qs["n_fps"]
    z_eng = qs["z_eng"]
    for i in range(1, len(ts)):
        if n_fps[i - 1] == 3 and n_fps[i] == 1 and z_eng[i] < 0:
            return int(ts[i])
    return None


def detect_jump(P_pi1, threshold=0.5, w=15):
    sm = causal_ma(P_pi1, w)
    for t in range(w, len(sm)):
        if sm[t] < threshold:
            end = min(t + w, len(sm))
            if np.all(sm[t:end] < threshold):
                return t
    return len(P_pi1)


def choose_display_agent(all_hist, jtimes, n_healthy, w=12):
    median_jump = np.median(jtimes)
    q25, q75 = np.percentile(jtimes, [25, 75])
    candidates = [
        idx for idx, jt in enumerate(jtimes)
        if q25 <= jt <= q75 and n_healthy < jt < len(all_hist[idx]["P_pi1"]) - w
    ]
    if not candidates:
        candidates = [
            idx for idx, jt in enumerate(jtimes)
            if n_healthy < jt < len(all_hist[idx]["P_pi1"]) - w
        ]

    best_idx = int(np.argsort(jtimes)[len(jtimes) // 2])
    best_score = -np.inf
    for idx in candidates:
        jt = int(jtimes[idx])
        sm = causal_ma(all_hist[idx]["P_pi1"], w)
        pre = np.nanmean(sm[max(n_healthy, jt - w):jt])
        post = np.nanmean(sm[jt:min(len(sm), jt + w)])
        drop = pre - post
        timing_penalty = abs(jt - median_jump) / max(median_jump, 1)
        score = drop - 0.15 * timing_penalty
        if np.isfinite(score) and score > best_score:
            best_score = score
            best_idx = idx
    return best_idx


# ====================================================================
# FIGURE P2 (revised): Onset + quasi-static analysis
#   Panel (b) is now split into three sub-panels — no ×30 scaling
# ====================================================================
def plot_onset_with_theory():
    """
    3×2 GridSpec layout:
      Col 0: (a) ensemble, (c) quasistatic, (d) variance
      Col 1: (b1) P(engaged), (b2) γ(t), (b3) Δc(t) — all shared x-axis
    """
    d = DEFAULTS.copy()
    N_total = 600
    n_agents = 100
    kw = {k: d[k] for k in ["gamma_healthy", "delta_c_healthy", "N_healthy",
                              "gamma_rate", "delta_c_rate", "gamma_floor", "delta_c_floor"]}

    all_hist = []
    jtimes = []
    for i in range(n_agents):
        h = run_agent_with_drift(d["p"], d["alpha_0"], N_total, seed=500 + i, **kw)
        all_hist.append(h)
        jtimes.append(detect_jump(h["P_pi1"]))
    jtimes = np.array(jtimes)

    display_idx = choose_display_agent(all_hist, jtimes, d["N_healthy"])
    h_display = all_hist[display_idx]
    qs = quasistatic_fps(h_display, d["p"], d["alpha_0"], every=2)
    fold_t = detect_fold(qs)

    ### DT ---> 3×2 GridSpec: left column holds ensemble, quasistatic, variance;
    ### DT ---> right column holds the three time-series panels for the median
    ### DT ---> agent (P(engaged), γ(t), Δc(t)) on a shared x-axis.
    fig = plt.figure(figsize=(11.2, 8.2))
    gs = GridSpec(3, 2, figure=fig, hspace=0.50, wspace=0.35)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b1 = fig.add_subplot(gs[0, 1])
    ax_b2 = fig.add_subplot(gs[1, 1], sharex=ax_b1)
    ax_b3 = fig.add_subplot(gs[2, 1], sharex=ax_b1)
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[2, 0])

    # (a) Ensemble + fold
    show = np.linspace(0, n_agents - 1, 25, dtype=int)
    for idx in show:
        ax_a.plot(np.arange(N_total), causal_ma(all_hist[idx]["P_pi1"], 15),
                  color="#9fb3c8", alpha=0.18, lw=0.8)
    pm = np.mean([h["P_pi1"] for h in all_hist], axis=0)
    ax_a.plot(np.arange(N_total), causal_ma(pm, 15), "k-", lw=2.5, label="Pop. mean")
    ax_a.axvline(x=d["N_healthy"], color="red", ls="--", alpha=0.5, label="Adversity onset")
    if fold_t is not None:
        ax_a.axvline(x=fold_t, color="purple", ls=":", lw=1.5,
                     label=f"Fold loss (t={fold_t})")
    ax_a.set_xlabel("Trial"); ax_a.set_ylabel("Latent P(engaged)")
    ax_a.set_title("(a) Ensemble trajectories under adversity drift")
    ax_a.set_ylim(-0.05, 1.05); ax_a.legend(fontsize=7, loc="lower left")

    # (b1) P(engaged) of illustrative jump agent
    jt = int(jtimes[display_idx])
    ax_b1.plot(np.arange(N_total), causal_ma(h_display["P_pi1"], 10),
               color="steelblue", lw=2.0, label="Latent P(engaged)")
    ax_b1.axvline(x=d["N_healthy"], color="red", ls="--", alpha=0.4, lw=1)
    ax_b1.axvline(x=jt, color="purple", ls=":", lw=1.5,
                  label=f"Observed transition (t={jt})")
    ax_b1.set_ylabel("P(engaged)")
    ax_b1.set_ylim(-0.05, 1.05)
    ax_b1.set_title("(b) Illustrative single-agent transition")
    ax_b1.legend(fontsize=7, loc="upper right")
    plt.setp(ax_b1.get_xticklabels(), visible=False)

    # (b2) γ(t)
    ax_b2.plot(np.arange(N_total), h_display["gamma_t"], color="crimson", lw=1.6)
    ax_b2.axvline(x=d["N_healthy"], color="red", ls="--", alpha=0.4, lw=1)
    ax_b2.axvline(x=jt, color="purple", ls=":", lw=1.5)
    ax_b2.set_ylabel(r"$\gamma(t)$")
    plt.setp(ax_b2.get_xticklabels(), visible=False)

    # (b3) Δc(t)
    ax_b3.plot(np.arange(N_total), h_display["delta_c_t"], color="darkorange", lw=1.6)
    ax_b3.axvline(x=d["N_healthy"], color="red", ls="--", alpha=0.4, lw=1,
                  label="Adversity onset")
    ax_b3.axvline(x=jt, color="purple", ls=":", lw=1.5)
    ax_b3.axhline(y=0, color="gray", ls="--", alpha=0.4, lw=0.8)
    ax_b3.set_xlabel("Trial"); ax_b3.set_ylabel(r"$\Delta c(t)$")
    ax_b3.legend(fontsize=7)

    # (c) Quasistatic fixed-point structure
    phi_eng = (1 + qs["z_eng"]) / 2
    phi_wdr = (1 + qs["z_wdr"]) / 2
    phi_uns = (1 + qs["z_uns"]) / 2
    ax_c.plot(qs["t"], phi_eng, color="darkgreen", lw=2, label="Stable engaged branch")
    ax_c.plot(qs["t"], phi_wdr, color="darkred", lw=2, label="Stable withdrawn branch")
    ax_c.plot(qs["t"], phi_uns, color="black", ls="--", lw=1.2, label="Unstable branch")
    ax_c.plot(np.arange(N_total), causal_ma(h_display["P_pi1"], 10),
              color="steelblue", lw=1.5, alpha=0.9, label="Observed latent P(engaged)")
    ax_c.axvline(x=d["N_healthy"], color="red", ls="--", alpha=0.3)
    if fold_t is not None:
        ax_c.axvline(x=fold_t, color="purple", ls=":", lw=1.5,
                     label=f"Engaged-branch fold loss (t={fold_t})")
    ax_c.set_xlabel("Trial"); ax_c.set_ylabel("Quasistatic branch / latent P(engaged)")
    ax_c.set_ylim(-0.05, 1.05)
    ax_c.set_title("(c) Quasistatic structure along the live schedule")
    ax_c.legend(fontsize=6, loc="center right")

    # (d) Population variance
    n_pop = 200
    allP = np.zeros((n_pop, N_total))
    for i in range(n_pop):
        hp = run_agent_with_drift(d["p"], d["alpha_0"], N_total, seed=1000 + i, **kw)
        allP[i, :] = hp["P_pi1"]
    pv = np.var(allP, axis=0); ppm = np.mean(allP, axis=0)
    axv = ax_d; axm = ax_d.twinx()
    axv.plot(np.arange(N_total), causal_ma(pv, 15), "r-", lw=2, label="Var[latent P(engaged)]")
    axm.plot(np.arange(N_total), causal_ma(ppm, 15), "b--", lw=2, label="Mean")
    ax_d.axvline(x=d["N_healthy"], color="red", ls="--", alpha=0.4)
    pk = int(np.nanargmax(causal_ma(pv, 15)[d["N_healthy"]:]) + d["N_healthy"])
    ax_d.axvline(x=pk, color="orange", ls=":", lw=2, label=f"Var peak (t={pk})")
    axv.set_xlabel("Trial"); axv.set_ylabel("Var", color="r"); axm.set_ylabel("Mean", color="b")
    ax_d.set_title("(d) Population susceptibility")
    axv.legend(fontsize=7, loc="upper left"); axm.legend(fontsize=7, loc="upper right")

    fig.suptitle("Onset under gradual adversity drift", fontsize=12, y=1.01)
    save_figure(fig, "fig_P2_onset_revised")
    plt.close()
    print("Saved fig_P2_onset_revised.[png|pdf]")
    return jtimes


# ====================================================================
# TRUE REVERSAL PROTOCOL
#   Not produced by run_all() but retained as a callable function.
# ====================================================================
def plot_true_reversal(unused=None):
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    d = DEFAULTS.copy()
    N_rec = 300
    n_agents = 80
    kw = {k: d[k] for k in ["gamma_healthy", "delta_c_healthy", "N_healthy",
                              "gamma_rate", "delta_c_rate", "gamma_floor", "delta_c_floor"]}
    t_revs = [350, 400, 450, 500]

    ### DT ---> Post-reversal observation window is fixed across residence-time
    ### DT ---> conditions so later reversals are not penalised by shorter recovery.
    def simulate_reversal(seed, t_rev, delta_c_reversal):
        N_total = t_rev + N_rec
        return run_agent_with_drift(
            d["p"], d["alpha_0"], N_total, seed=seed, **kw,
            t_reversal=t_rev, gamma_reversal=d["gamma_healthy"],
            delta_c_reversal=delta_c_reversal
        )

    ax = axes[0]
    t_rev = 400
    N_total = t_rev + N_rec
    for i in range(12):
        h = simulate_reversal(7000 + i, t_rev=t_rev, delta_c_reversal=d["delta_c_healthy"])
        ax.plot(np.arange(N_total), causal_ma(h["P_pi1"], 15), alpha=0.4, lw=0.9)
    ax.axvline(x=d["N_healthy"], color="red", ls="--", alpha=0.4, label="Adversity")
    ax.axvline(x=t_rev, color="green", ls="--", lw=2, label=f"Full reversal t={t_rev}")
    ax.set_xlabel("Trial"); ax.set_ylabel("Latent P(engaged)")
    ax.set_title(r"(a) Onset $\to$ full reversal ($\gamma$ + $\Delta c$)"
                 "\n(matched 300-trial recovery window)")
    ax.set_ylim(-0.05, 1.05); ax.legend(fontsize=7)

    ax = axes[1]
    rr_full = []; rr_gam = []
    for t_rev in t_revs:
        nf = ng = 0
        for i in range(n_agents):
            h = simulate_reversal(8000 + i, t_rev=t_rev, delta_c_reversal=d["delta_c_healthy"])
            if np.mean(h["P_pi1"][-50:]) > 0.6: nf += 1
            h2 = simulate_reversal(8000 + i, t_rev=t_rev, delta_c_reversal=None)
            if np.mean(h2["P_pi1"][-50:]) > 0.6: ng += 1
        rr_full.append(nf / n_agents); rr_gam.append(ng / n_agents)
    ax.plot(t_revs, rr_full, "bo-", lw=2, label=r"Restore $\gamma$ + $\Delta c$")
    ax.plot(t_revs, rr_gam, "rs--", lw=2, label=r"Restore $\gamma$ only")
    ax.set_xlabel("Reversal time"); ax.set_ylabel("Recovery rate")
    ax.set_title("(b) Recovery rate vs residence time\n(matched recovery window)")
    ax.legend(fontsize=8); ax.set_ylim(-0.05, 1.05)

    ax = axes[2]
    for lab, dc_r, col, mk in [
        (r"$\gamma$ + $\Delta c$", d["delta_c_healthy"], "blue", "o"),
        (r"$\gamma$ only", None, "red", "s")]:
        ms = []; ses = []
        for t_rev in t_revs:
            lv = []
            for i in range(n_agents):
                h = simulate_reversal(8000 + i, t_rev=t_rev, delta_c_reversal=dc_r)
                lv.append(np.mean(h["P_pi1"][-50:]))
            lv = np.array(lv)
            ms.append(np.mean(lv)); ses.append(np.std(lv) / np.sqrt(len(lv)))
        ax.errorbar(t_revs, ms, yerr=ses, fmt=f"{mk}-", color=col, lw=2, capsize=4, label=lab)
    ax.axhline(y=0.6, color="gray", ls=":", alpha=0.5, label="Recovery threshold")
    ax.set_xlabel("Reversal time"); ax.set_ylabel("Mean late P(engaged) ±SE")
    ax.set_title("(c) Hysteretic asymmetry: same-variable reversal\n(matched recovery window)")
    ax.legend(fontsize=7); ax.set_ylim(-0.05, 1.05)

    fig.suptitle("True reversal protocol (actual onset schedule, then attempted recovery)\n"
                 "(all residence-time conditions use the same post-reversal recovery window)",
                 fontsize=13, y=1.03)
    plt.tight_layout()
    plt.savefig(f"{OUT}/fig_P_reversal_true.png")
    plt.close()
    print("Saved fig_P_reversal_true.png")


# ====================================================================
# SHARPENED FLAGS  (2 panels: divergence, sudden jumps)
#   - Bimodality panel dropped: bimodality is already shown via the
#     variance peak in fig_P2 panel (d); the histogram is auxiliary.
#   - Inaccessible-region annotation dropped: not a core claim here.
#   - Divergence: paired trajectory bifurcation from similar Δc₀ values
#   - Sudden jumps: individual trajectories showing rapid commitment
# ====================================================================
def plot_sharpened_flags():
    """
    Two catastrophe flags in a 1×2 layout.

    Flag 5 (divergence): pairs of agents with Δc₀ just above and just
    below the catastrophe set diverge to opposite attractors — showing
    how small initial differences produce qualitatively different long-run
    outcomes.

    Flag 3 (sudden jumps): individual trajectories showing rapid
    commitment once the fold is crossed.
    """
    d = DEFAULTS.copy()
    N_total = 600
    kw = {k: d[k] for k in ["gamma_healthy", "delta_c_healthy", "N_healthy",
                              "gamma_rate", "delta_c_rate", "gamma_floor", "delta_c_floor"]}

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))

    # --- Flag 5: genuine divergence from similar initial conditions ---
    ### DT ---> Two Δc₀ values: "high" (just above fold, stays engaged) and
    ### DT ---> "low" (just below fold, transitions to withdrawn).
    ### DT ---> Values chosen by inspecting the quasistatic fixed-point structure.
    ax = axes[0]
    dc_high = 0.12   # above the bifurcation — stays engaged at t=500
    dc_low = 0.04    # below the bifurcation — transitions to withdrawn
    n_pairs = 5
    colours_high = plt.cm.Greens(np.linspace(0.5, 0.9, n_pairs))
    colours_low = plt.cm.Reds(np.linspace(0.5, 0.9, n_pairs))
    for i in range(n_pairs):
        for dc_val, cols, ls, lbl in [
            (dc_high, colours_high, '-', r'$\Delta c_0=0.12$ (stays engaged)' if i == 0 else '_'),
            (dc_low, colours_low, '--', r'$\Delta c_0=0.04$ (transitions)' if i == 0 else '_'),
        ]:
            h = run_agent_with_drift(
                d["p"], d["alpha_0"], N_total, seed=9000 + i,
                gamma_healthy=d["gamma_healthy"], delta_c_healthy=dc_val,
                N_healthy=d["N_healthy"], gamma_rate=d["gamma_rate"],
                delta_c_rate=d["delta_c_rate"], gamma_floor=d["gamma_floor"],
                delta_c_floor=-0.20
            )
            ax.plot(np.arange(N_total), causal_ma(h["P_pi1"], 15),
                    color=cols[i], lw=1.5, ls=ls, alpha=0.85, label=lbl)

    ax.axvline(x=d["N_healthy"], color="red", ls="--", alpha=0.4, label="Adversity onset")
    ax.set_xlabel("Trial")
    ax.set_ylabel("Latent P(engaged)")
    ax.set_ylim(-0.05, 1.05)
    ax.set_title("(a) Divergence: small difference in initial engagement value\n"
                 r"separates trajectories to opposite attractors ($n=5$ pairs)")
    ax.legend(fontsize=8)

    # --- Flag 3: sudden jumps ---
    ax = axes[1]
    for i in range(10):
        h = run_agent_with_drift(d["p"], d["alpha_0"], N_total, seed=4000 + i * 7, **kw)
        ax.plot(np.arange(N_total), causal_ma(h["P_pi1"], 15), lw=1, alpha=0.65)
    ax.axvline(x=d["N_healthy"], color="red", ls="--", alpha=0.4, label="Adversity onset")
    ax.set_xlabel("Trial")
    ax.set_ylabel("Latent P(engaged)")
    ax.set_title("(b) Sudden jumps in individual trajectories\n"
                 "once the fold threshold is crossed")
    ax.legend(fontsize=8)
    ax.set_ylim(-0.05, 1.05)

    fig.suptitle("Two catastrophe signatures: trajectory divergence and sudden jumps\n"
                 "(latent variable; bimodality and variance shown in fig_P2)",
                 fontsize=12, y=1.04)
    plt.tight_layout()
    plt.savefig(f"{OUT}/fig_P_flags_revised.png")
    plt.close()
    print("Saved fig_P_flags_revised.png")


# ====================================================================
# EWS: event-aligned autocorrelation
#   Replaces population-averaged AC (which decreases after adversity onset
#   due to population heterogeneity in jump timing) with event-aligned AC:
#   each agent's rolling lag-1 AC is aligned at that agent's own transition
#   time, then averaged across agents.  This reveals the pre-transition
#   rise (critical slowing down) that is masked in the population average.
# ====================================================================
def compute_event_aligned_ac(allP, jtimes, window_before=80, window_after=40,
                              roll_w=25, min_agents=15):
    """
    Compute event-aligned rolling lag-1 autocorrelation.

    For each agent that has a finite transition time, compute rolling
    lag-1 AC at each trial, align at the agent's transition time (t_rel=0),
    and average across agents.

    Returns (t_rel, mean_ac, se_ac, n_valid).
    """
    n_agents, N_total = allP.shape
    rel_range = np.arange(-window_before, window_after)
    aligned = np.full((n_agents, len(rel_range)), np.nan)

    for i, jt in enumerate(jtimes):
        if jt >= N_total:   # agent never transitioned
            continue
        tr = allP[i]
        for t_idx, t_rel in enumerate(rel_range):
            t_abs = jt + t_rel
            if t_abs < roll_w or t_abs >= N_total - 1:
                continue
            seg = tr[t_abs - roll_w: t_abs]
            sd = np.std(seg)
            if sd > 1e-8:
                aligned[i, t_idx] = np.corrcoef(seg[:-1], seg[1:])[0, 1]

    n_valid = np.sum(~np.isnan(aligned), axis=0)
    mean_ac = np.where(n_valid >= min_agents, np.nanmean(aligned, axis=0), np.nan)
    se_ac = np.where(n_valid >= min_agents,
                     np.nanstd(aligned, axis=0) / np.sqrt(np.maximum(n_valid, 1)),
                     np.nan)
    return rel_range, mean_ac, se_ac, n_valid


def plot_ews_revised():
    """
    Two-panel EWS figure.

    Panel (a): event-aligned lag-1 autocorrelation — averaged across agents
    after aligning each at their own transition time.  Shows the pre-
    transition rise in AC (critical slowing down) that is washed out in the
    population-average because agents transition at different times.

    Panel (b): population susceptibility (variance of latent P(engaged))
    with population mean overlaid.
    """
    d = DEFAULTS.copy()
    N_total = 600; n_ag = 200
    kw = {k: d[k] for k in ["gamma_healthy", "delta_c_healthy", "N_healthy",
                              "gamma_rate", "delta_c_rate", "gamma_floor", "delta_c_floor"]}
    allP = np.zeros((n_ag, N_total))
    jtimes_ews = np.zeros(n_ag, dtype=int)
    for i in range(n_ag):
        h = run_agent_with_drift(d["p"], d["alpha_0"], N_total, seed=2000 + i, **kw)
        allP[i, :] = h["P_pi1"]
        jtimes_ews[i] = detect_jump(h["P_pi1"])

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    # --- Panel (a): event-aligned autocorrelation ---
    ax = axes[0]
    ### DT ---> Event-aligned AC: align each agent at its own transition time
    ### DT ---> so the pre-transition rise is not masked by population spread
    ### DT ---> in jump timing (which causes the population-average AC to fall).
    t_rel, mean_ac, se_ac, n_valid = compute_event_aligned_ac(
        allP, jtimes_ews, window_before=80, window_after=40, roll_w=25
    )
    valid = ~np.isnan(mean_ac)
    ax.fill_between(t_rel[valid],
                    mean_ac[valid] - se_ac[valid],
                    mean_ac[valid] + se_ac[valid],
                    alpha=0.25, color="purple")
    ax.plot(t_rel[valid], mean_ac[valid], color="purple", lw=2,
            label="Mean event-aligned lag-1 AC")
    ax.axvline(x=0, color="red", ls="--", lw=1.5, alpha=0.7,
               label="Transition (t_rel = 0)")
    ax.axhline(y=0, color="gray", ls=":", alpha=0.4)
    ax.set_xlabel("Trials relative to individual transition")
    ax.set_ylabel("Lag-1 autocorrelation of latent P(engaged)")
    ax.set_title("(a) Pre-transition rise in lag-1 AC (critical slowing down)\n"
                 "revealed by event-alignment at each agent's own transition\n"
                 f"(shading = ±1 SE; n agents at midpoint: "
                 f"{n_valid[len(t_rel)//2]} of {n_ag})")
    ax.legend(fontsize=8)

    # Annotate pre-transition trend
    pre = t_rel < 0
    if np.any(valid & pre):
        pre_vals = mean_ac[valid & pre]
        pre_ts = t_rel[valid & pre]
        if len(pre_vals) > 5:
            slope = np.polyfit(pre_ts[-30:], pre_vals[-30:], 1)[0]
            ax.text(0.05, 0.92,
                    f"Pre-transition slope: {slope:+.4f} per trial",
                    transform=ax.transAxes, fontsize=8, color="purple")

    # --- Panel (b): population susceptibility ---
    ax = axes[1]
    cpv = np.var(allP, axis=0)
    ax.plot(np.arange(N_total), causal_ma(cpv, 15), "darkred", lw=2, label="Var")
    ax2 = ax.twinx()
    cpm = np.mean(allP, axis=0)
    ax2.plot(np.arange(N_total), causal_ma(cpm, 15), "b--", lw=1.5, label="Mean")
    ax.axvline(x=d["N_healthy"], color="red", ls="--", alpha=0.5,
               label="Adversity onset")
    ax.set_xlabel("Trial")
    ax.set_ylabel("Var of latent P(engaged)", color="darkred")
    ax2.set_ylabel("Mean latent P(engaged)", color="b")
    ax.set_title("(b) Population susceptibility (latent variable)")
    ax.legend(fontsize=7, loc="upper left"); ax2.legend(fontsize=7, loc="center right")

    fig.suptitle(
        "Early warning signals: event-aligned AC unmasks pre-transition critical slowing down\n"
        "(population-average AC without alignment is masked by heterogeneous jump timing)",
        fontsize=12, y=1.04
    )
    plt.tight_layout()
    plt.savefig(f"{OUT}/fig_P3_ews_revised.png")
    plt.close()
    print("Saved fig_P3_ews_revised.png")


# ====================================================================
# Entry point called by main.py
# ====================================================================
def run_all():
    """
    Produces the three Step 1 publication figures:
      - fig_P2_onset_revised.png   (ensemble, schedule panels, quasistatic, variance)
      - fig_P_flags_revised.png    (bimodality, divergence, sudden jumps)
      - fig_P3_ews_revised.png     (event-aligned AC, susceptibility)

    fig_P_reversal_true is not produced here (retained as a callable function
    for ad-hoc use).
    """
    print("Running Step 1: onset dynamics, catastrophe flags, EWS...")
    jtimes = plot_onset_with_theory()
    jv = jtimes[jtimes < 600]
    if len(jv) > 0:
        print(f"  Jump times: median={np.median(jv):.0f}, "
              f"IQR=[{np.percentile(jv, 25):.0f}, {np.percentile(jv, 75):.0f}]")
    plot_sharpened_flags()
    plot_ews_revised()


# ====================================================================
# MAIN
# ====================================================================
if __name__ == "__main__":
    print("=" * 65)
    print("Step 1 (consolidated): all review items addressed")
    print("=" * 65)

    print("\n1. Onset + commitment analysis...")
    jtimes = plot_onset_with_theory()
    jv = jtimes[jtimes < 600]
    if len(jv) > 0:
        print(f"   Jump times: median={np.median(jv):.0f}, "
              f"IQR=[{np.percentile(jv, 25):.0f}, {np.percentile(jv, 75):.0f}]")

    print("\n2. Sharpened flags...")
    plot_sharpened_flags()

    print("\n3. EWS (event-aligned autocorrelation)...")
    plot_ews_revised()

    print("\n4. True reversal protocol (supplementary, not in main sweep)...")
    plot_true_reversal(None)

    print("\n" + "=" * 65)
    print(f"All saved to {OUT}")
    print("=" * 65)
