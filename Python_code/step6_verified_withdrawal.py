"""Verified-withdrawal sensitivity analyses for the recovery claims.

The main manuscript figures use a fixed adversity-exposure schedule.  This
module leaves those figures unchanged and answers two narrower questions:

1. What state have agents actually reached at the fixed intervention trial?
2. Does the field-versus-precision ordering, and the asymmetric-memory result,
   persist when intervention is anchored to a verified withdrawal state?

Withdrawal is verified only after both adverse controls have reached their
floors and latent P(engaged) has remained below .50 for 30 consecutive trials.
The duration sensitivity then counts additional adverse-schedule trials from
the end of that confirmation window; it does not require the stochastic latent
probability to remain below .50 at every subsequent trial. Intervention
conditions branch from identical learned evidence and identical random-number-
generator states, giving a paired comparison.
"""

import copy
import csv
import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from core_functions import compute_efe_difference, policy_posterior_pi1
from psychopathology_regimes import recovery_regime


OUT = os.path.join(os.path.dirname(__file__), "Figs_psychopathology")
TABLE_OUT = os.path.join(os.path.dirname(__file__), "..", "outputs", "tables")
os.makedirs(OUT, exist_ok=True)
os.makedirs(TABLE_OUT, exist_ok=True)

REGIME = recovery_regime()

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


def _initial_state(alpha_0):
    base = alpha_0 / 2.0
    return {
        "al1": base,
        "be1": base,
        "al2": base,
        "be2": base,
        "n_engaged": 0,
        "n_withdrawn": 0,
    }


def _controls_at(t, regime=REGIME):
    if t < regime["N_healthy"]:
        return regime["gamma_healthy"], regime["delta_c_healthy"]
    dt = t - regime["N_healthy"]
    gamma = max(
        regime["gamma_healthy"] - regime["gamma_rate"] * dt,
        regime["gamma_floor"],
    )
    delta_c = max(
        regime["delta_c_healthy"] - regime["delta_c_rate"] * dt,
        regime["delta_c_floor"],
    )
    return gamma, delta_c


def adverse_floor_trial(regime=REGIME):
    gamma_steps = (
        (regime["gamma_healthy"] - regime["gamma_floor"])
        / regime["gamma_rate"]
        if regime["gamma_rate"] > 0
        else 0
    )
    field_steps = (
        (regime["delta_c_healthy"] - regime["delta_c_floor"])
        / regime["delta_c_rate"]
        if regime["delta_c_rate"] > 0
        else 0
    )
    return regime["N_healthy"] + int(math.ceil(max(gamma_steps, field_steps)))


def _latent_probability(state, gamma, delta_c):
    a = state["al1"] / (state["al1"] + state["be1"])
    b = state["be2"] / (state["al2"] + state["be2"])
    delta_g = compute_efe_difference(a, b, delta_c)
    return float(policy_posterior_pi1(delta_g, gamma))


def _advance_one(
    state,
    rng,
    gamma,
    delta_c,
    eta_eng=1.0,
    eta_wdr=1.0,
    regime=REGIME,
):
    base = regime["alpha_0"] / 2.0
    if eta_eng < 1.0 or eta_wdr < 1.0:
        state["al1"] = eta_eng * state["al1"] + (1.0 - eta_eng) * base
        state["be1"] = eta_eng * state["be1"] + (1.0 - eta_eng) * base
        state["al2"] = eta_wdr * state["al2"] + (1.0 - eta_wdr) * base
        state["be2"] = eta_wdr * state["be2"] + (1.0 - eta_wdr) * base

    p_engaged = _latent_probability(state, gamma, delta_c)
    action = 0 if rng.random() < p_engaged else 1
    p = regime["p"]

    if action == 0:
        observation = 0 if rng.random() < p else 1
        if observation == 0:
            state["al1"] += 1.0
        else:
            state["be1"] += 1.0
        state["n_engaged"] += 1
    else:
        observation = 0 if rng.random() < (1.0 - p) else 1
        if observation == 0:
            state["al2"] += 1.0
        else:
            state["be2"] += 1.0
        state["n_withdrawn"] += 1

    return p_engaged


def fixed_schedule_state(seed, t_intervention=None, regime=REGIME):
    """State diagnostic for the fixed schedule used in the main recovery figure."""
    if t_intervention is None:
        t_intervention = regime["N_healthy"] + 200
    rng = np.random.default_rng(seed)
    state = _initial_state(regime["alpha_0"])
    history = []

    for t in range(t_intervention):
        gamma, delta_c = _controls_at(t, regime)
        history.append(_advance_one(state, rng, gamma, delta_c, regime=regime))

    gamma, delta_c = _controls_at(t_intervention, regime)
    p_intervention = _latent_probability(state, gamma, delta_c)
    total_actions = state["n_engaged"] + state["n_withdrawn"]
    z = (
        (state["n_engaged"] - state["n_withdrawn"]) / total_actions
        if total_actions
        else 0.0
    )
    return {
        "seed": seed,
        "intervention_trial": t_intervention,
        "adversity_exposure_trials": t_intervention - regime["N_healthy"],
        "gamma_at_intervention": gamma,
        "delta_c_at_intervention": delta_c,
        "p_engaged_at_intervention": p_intervention,
        "mean_p_engaged_prior_15": float(np.mean(history[-15:])),
        "mean_p_engaged_prior_100": float(np.mean(history[-100:])),
        "order_parameter_z": z,
    }


def verified_withdrawal_state(
    seed,
    residence_trials=200,
    confirm_trials=30,
    eta_eng=1.0,
    eta_wdr=1.0,
    max_trials=3000,
    regime=REGIME,
):
    """Return a paired branch state after confirmed withdrawal and a delay."""
    rng = np.random.default_rng(seed)
    state = _initial_state(regime["alpha_0"])
    floor_trial = adverse_floor_trial(regime)
    below_run = 0
    verification_trial = None
    intervention_trial = None
    post_verification_below = []
    recent_probabilities = []

    for t in range(max_trials):
        gamma, delta_c = _controls_at(t, regime)
        p_engaged = _advance_one(
            state,
            rng,
            gamma,
            delta_c,
            eta_eng=eta_eng,
            eta_wdr=eta_wdr,
            regime=regime,
        )
        recent_probabilities.append(p_engaged)
        if len(recent_probabilities) > confirm_trials:
            recent_probabilities.pop(0)

        if verification_trial is None:
            if t >= floor_trial and p_engaged < 0.5:
                below_run += 1
            else:
                below_run = 0
            if below_run >= confirm_trials:
                verification_trial = t + 1
                intervention_trial = verification_trial + residence_trials
        else:
            post_verification_below.append(float(p_engaged < 0.5))

        if verification_trial is not None and t + 1 >= intervention_trial:
            p_at_intervention = _latent_probability(
                state, regime["gamma_floor"], regime["delta_c_floor"]
            )
            recent_mean = float(np.mean(recent_probabilities))
            if p_at_intervention >= 0.5 or recent_mean >= 0.5:
                continue
            total_actions = state["n_engaged"] + state["n_withdrawn"]
            z = (
                (state["n_engaged"] - state["n_withdrawn"]) / total_actions
                if total_actions
                else 0.0
            )
            return {
                "seed": seed,
                "state": copy.deepcopy(state),
                "rng_state": copy.deepcopy(rng.bit_generator.state),
                "floor_trial": floor_trial,
                "verification_trial": verification_trial,
                "intervention_trial": t + 1,
                "requested_post_verification_adverse_trials": residence_trials,
                "actual_post_verification_trials": t + 1 - verification_trial,
                "confirm_trials": confirm_trials,
                "p_engaged_at_intervention": p_at_intervention,
                "mean_p_engaged_prior_confirmation_window": recent_mean,
                "post_verification_proportion_p_below_0_5": (
                    float(np.mean(post_verification_below))
                    if post_verification_below
                    else float("nan")
                ),
                "order_parameter_z": z,
                "eta_eng": eta_eng,
                "eta_wdr": eta_wdr,
            }

    raise RuntimeError(
        f"Seed {seed} did not meet the verified-withdrawal criterion by "
        f"trial {max_trials}."
    )


def continue_from_verified_state(
    branch,
    gamma,
    delta_c,
    eta_eng=1.0,
    eta_wdr=1.0,
    recovery_trials=400,
    late_window=100,
    regime=REGIME,
):
    state = copy.deepcopy(branch["state"])
    rng = np.random.default_rng()
    rng.bit_generator.state = copy.deepcopy(branch["rng_state"])
    probabilities = []

    for _ in range(recovery_trials):
        probabilities.append(
            _advance_one(
                state,
                rng,
                gamma,
                delta_c,
                eta_eng=eta_eng,
                eta_wdr=eta_wdr,
                regime=regime,
            )
        )

    return float(np.mean(probabilities[-late_window:]))


def _write_rows(path, rows):
    if not rows:
        raise ValueError(f"No rows supplied for {path}")
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _summarise(values):
    values = np.asarray(values, dtype=float)
    return {
        "median": float(np.median(values)),
        "q1": float(np.percentile(values, 25)),
        "q3": float(np.percentile(values, 75)),
        "mean": float(np.mean(values)),
        "sd": float(np.std(values, ddof=1)),
        "proportion_above_0_5": float(np.mean(values > 0.5)),
    }


def run_fixed_schedule_diagnostic(n_agents=100):
    rows = [fixed_schedule_state(12000 + i) for i in range(n_agents)]
    _write_rows(
        os.path.join(TABLE_OUT, "fixed_schedule_intervention_state.csv"), rows
    )
    return rows


def run_verified_field_dominance(n_agents=200, residence_trials=200):
    conditions = [
        ("Control", REGIME["gamma_floor"], REGIME["delta_c_floor"]),
        ("Restore gamma only", REGIME["gamma_healthy"], REGIME["delta_c_floor"]),
        ("Restore delta_c only", REGIME["gamma_floor"], REGIME["delta_c_healthy"]),
        ("Restore both", REGIME["gamma_healthy"], REGIME["delta_c_healthy"]),
    ]
    long_rows = []
    diagnostic_rows = []

    for i in range(n_agents):
        branch = verified_withdrawal_state(
            130000 + i, residence_trials=residence_trials
        )
        diagnostic_rows.append({
            key: branch[key]
            for key in [
                "seed",
                "floor_trial",
                "verification_trial",
                "intervention_trial",
                "requested_post_verification_adverse_trials",
                "actual_post_verification_trials",
                "confirm_trials",
                "p_engaged_at_intervention",
                "mean_p_engaged_prior_confirmation_window",
                "post_verification_proportion_p_below_0_5",
                "order_parameter_z",
            ]
        })
        for condition, gamma, delta_c in conditions:
            late_p = continue_from_verified_state(
                branch, gamma=gamma, delta_c=delta_c
            )
            long_rows.append({
                "seed": branch["seed"],
                "condition": condition,
                "gamma_restored": gamma,
                "delta_c_restored": delta_c,
                "post_verification_adverse_trials": residence_trials,
                "late_p_engaged": late_p,
            })

    summary_rows = []
    for condition, _, _ in conditions:
        vals = [r["late_p_engaged"] for r in long_rows if r["condition"] == condition]
        summary_rows.append({
            "condition": condition,
            "n_agents": n_agents,
            **_summarise(vals),
            "proportion_positive": "",
        })

    gamma_by_seed = {
        r["seed"]: r["late_p_engaged"]
        for r in long_rows
        if r["condition"] == "Restore gamma only"
    }
    field_by_seed = {
        r["seed"]: r["late_p_engaged"]
        for r in long_rows
        if r["condition"] == "Restore delta_c only"
    }
    paired_advantage = np.array([
        field_by_seed[seed] - gamma_by_seed[seed] for seed in sorted(gamma_by_seed)
    ])
    advantage_summary = _summarise(paired_advantage)
    advantage_summary["proportion_above_0_5"] = ""
    summary_rows.append({
        "condition": "Paired delta_c advantage over gamma",
        "n_agents": n_agents,
        **advantage_summary,
        "proportion_positive": float(np.mean(paired_advantage > 0.0)),
    })

    _write_rows(
        os.path.join(TABLE_OUT, "verified_withdrawal_field_long.csv"), long_rows
    )
    _write_rows(
        os.path.join(TABLE_OUT, "verified_withdrawal_field_summary.csv"),
        summary_rows,
    )
    _write_rows(
        os.path.join(TABLE_OUT, "verified_withdrawal_diagnostics.csv"),
        diagnostic_rows,
    )
    return long_rows, summary_rows, diagnostic_rows


def run_verified_memory_sensitivity(
    n_agents=200,
    n_threshold=120,
    residence_values=(0, 50, 100, 200, 400),
):
    conditions = [
        ("Permanent evidence", 1.0, 1.0),
        ("Symmetric decay", 0.992, 0.992),
        ("Asymmetric mild", 0.990, 0.997),
        ("Asymmetric strong", 0.985, 0.997),
    ]
    recovery_rows = []
    threshold_rows = []
    dc_grid = np.arange(-0.02, 0.1201, 0.005)

    for ci, (condition, eta_eng, eta_wdr) in enumerate(conditions):
        for residence in residence_values:
            branches = [
                verified_withdrawal_state(
                    200000 + ci * 10000 + i,
                    residence_trials=residence,
                    eta_eng=eta_eng,
                    eta_wdr=eta_wdr,
                )
                for i in range(n_agents)
            ]
            full_values = [
                continue_from_verified_state(
                    branch,
                    gamma=REGIME["gamma_healthy"],
                    delta_c=REGIME["delta_c_healthy"],
                    eta_eng=eta_eng,
                    eta_wdr=eta_wdr,
                )
                for branch in branches
            ]
            recovery_rows.append({
                "condition": condition,
                "eta_eng": eta_eng,
                "eta_wdr": eta_wdr,
                "post_verification_adverse_trials": residence,
                "n_agents": n_agents,
                **_summarise(full_values),
            })

            if condition == "Permanent evidence":
                continue

            threshold_branches = branches[:n_threshold]
            grid_medians = []
            for delta_c in dc_grid:
                values = [
                    continue_from_verified_state(
                        branch,
                        gamma=REGIME["gamma_healthy"],
                        delta_c=float(delta_c),
                        eta_eng=eta_eng,
                        eta_wdr=eta_wdr,
                    )
                    for branch in threshold_branches
                ]
                grid_medians.append(float(np.median(values)))
            above = np.flatnonzero(np.asarray(grid_medians) >= 0.5)
            min_delta_c = (
                float(np.round(dc_grid[above[0]], 3))
                if len(above) else float("nan")
            )
            threshold_rows.append({
                "condition": condition,
                "eta_eng": eta_eng,
                "eta_wdr": eta_wdr,
                "post_verification_adverse_trials": residence,
                "n_agents_per_grid_cell": n_threshold,
                "min_restored_delta_c": min_delta_c,
            })

    _write_rows(
        os.path.join(TABLE_OUT, "verified_withdrawal_memory_summary.csv"),
        recovery_rows,
    )
    _write_rows(
        os.path.join(TABLE_OUT, "verified_withdrawal_threshold_summary.csv"),
        threshold_rows,
    )
    return recovery_rows, threshold_rows


def make_figure(fixed_rows, field_rows, threshold_rows):
    fig, axes = plt.subplots(1, 3, figsize=(12.3, 4.25))

    ax = axes[0]
    diagnostic_specs = [
        ("Current\nP < .50", "p_engaged_at_intervention", 0.5),
        ("Prior-15 mean\nP < .50", "mean_p_engaged_prior_15", 0.5),
        ("Order parameter\nz < 0", "order_parameter_z", 0.0),
    ]
    diagnostic_props = [
        np.mean([row[key] < threshold for row in fixed_rows])
        for _, key, threshold in diagnostic_specs
    ]
    bars = ax.bar(
        range(3), diagnostic_props, color=["#607D8B", "#78909C", "#90A4AE"],
        edgecolor="black", linewidth=0.7
    )
    for bar, value in zip(bars, diagnostic_props):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value + 0.035,
            f"{value:.0%}",
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold",
        )
    ax.set_xticks(range(3))
    ax.set_xticklabels([item[0] for item in diagnostic_specs], fontsize=8)
    ax.set_ylabel("Proportion of agents")
    ax.set_ylim(0, 1.05)
    median_p = np.median([row["p_engaged_at_intervention"] for row in fixed_rows])
    ax.set_title("(a) Fixed-schedule intervention\nis not state-anchored")
    ax.text(
        0.03,
        0.96,
        f"Median current P(engaged) = {median_p:.3f}",
        transform=ax.transAxes,
        va="top",
        fontsize=8,
    )

    ax = axes[1]
    condition_order = [
        "Control",
        "Restore gamma only",
        "Restore delta_c only",
        "Restore both",
    ]
    labels = ["Control", r"$\gamma$ only", r"$\Delta c$ only", "Both"]
    colours = ["#8A8A8A", "#4C78A8", "#F28E2B", "#2E7D32"]
    for idx, (condition, colour) in enumerate(zip(condition_order, colours)):
        vals = np.array([
            row["late_p_engaged"] for row in field_rows if row["condition"] == condition
        ])
        med = np.median(vals)
        q1, q3 = np.percentile(vals, [25, 75])
        ax.vlines(idx, q1, q3, color="black", lw=1.8, zorder=2)
        ax.scatter(
            idx, med, s=75, color=colour, edgecolor="black", linewidth=0.7, zorder=3
        )
    ax.axhline(0.5, color="#777777", ls="--", lw=0.9)
    ax.set_xticks(range(4))
    ax.set_xticklabels(labels)
    ax.set_ylabel("Late P(engaged), median [IQR]")
    ax.set_ylim(0, 1.05)
    ax.set_title("(b) Recovery after verified withdrawal\nand 200 further adverse trials")

    ax = axes[2]
    line_specs = [
        ("Symmetric decay", "#2CA25F", "--"),
        ("Asymmetric mild", "#F28E2B", "-"),
        ("Asymmetric strong", "#7B3294", "-"),
    ]
    for condition, colour, linestyle in line_specs:
        rows = sorted(
            [row for row in threshold_rows if row["condition"] == condition],
            key=lambda row: row["post_verification_adverse_trials"],
        )
        ax.plot(
            [row["post_verification_adverse_trials"] for row in rows],
            [row["min_restored_delta_c"] for row in rows],
            marker="o",
            ms=4.5,
            lw=1.8,
            linestyle=linestyle,
            color=colour,
            label=condition,
        )
    ax.axhline(
        REGIME["delta_c_healthy"], color="#777777", ls=":", lw=1.0,
        label=r"Healthy $\Delta c$"
    )
    ax.set_xlabel("Post-verification adverse-schedule trials")
    ax.set_ylabel(r"Min. restored $\Delta c$ for median P $\geq$ .50")
    ax.set_ylim(-0.025, 0.19)
    ax.set_title("(c) Event-anchored memory sensitivity")
    ax.legend(frameon=False, fontsize=7.2, loc="upper left")

    fig.suptitle(
        "State-anchored sensitivity analyses for the recovery simulations",
        fontsize=12,
        y=1.02,
    )
    plt.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_S6_verified_withdrawal.png"))
    fig.savefig(os.path.join(OUT, "fig_S6_verified_withdrawal.pdf"))
    plt.close(fig)


def run_all():
    print("Running fixed-schedule state diagnostic...")
    fixed_rows = run_fixed_schedule_diagnostic()
    print("Running paired field-versus-precision analysis after verified withdrawal...")
    field_rows, field_summary, diagnostic_rows = run_verified_field_dominance()
    print("Running event-anchored memory sensitivity...")
    memory_rows, threshold_rows = run_verified_memory_sensitivity()
    make_figure(fixed_rows, field_rows, threshold_rows)
    print("Saved fig_S6_verified_withdrawal.[png|pdf] and five CSV tables.")
    return {
        "fixed": fixed_rows,
        "field": field_rows,
        "field_summary": field_summary,
        "diagnostics": diagnostic_rows,
        "memory": memory_rows,
        "thresholds": threshold_rows,
    }


if __name__ == "__main__":
    run_all()
