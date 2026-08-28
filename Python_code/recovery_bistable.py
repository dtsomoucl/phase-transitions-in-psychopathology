"""Bistable recovery, hysteresis, and retention analyses.

The reference protocol is selected by structural criteria rather than by a
desired recovery outcome. At intervention, the adverse control condition has
two stable mean-field branches separated by an unstable branch. Restoring the
preference field removes the withdrawn branch. Stochastic intervention
comparisons are paired by branching from the same verified learned state and
random-number-generator state.
"""

from __future__ import annotations

import copy
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from core_functions import (
    compute_efe_difference,
    coupling_function,
    find_fixed_points,
    policy_posterior_pi1,
)


ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = ROOT / "outputs" / "tables"
MANUSCRIPT_FIGURE_DIR = ROOT / "outputs" / "manuscript_figures"

for directory in (TABLE_DIR, MANUSCRIPT_FIGURE_DIR):
    directory.mkdir(parents=True, exist_ok=True)

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


@dataclass(frozen=True)
class RecoveryRegime:
    p: float = 0.85
    alpha_0: float = 40.0
    gamma_healthy: float = 20.0
    gamma_floor: float = 18.0
    delta_c_healthy: float = 0.15
    delta_c_floor: float = -0.05
    healthy_trials: int = 100
    drift_trials: int = 100
    intervention_trial: int = 300
    confirm_trials: int = 30
    recovery_trials: int = 250
    late_window: int = 100


REGIME = RecoveryRegime()


@dataclass
class Branch:
    seed: int
    state: np.ndarray
    rng_state: dict
    probabilities: np.ndarray
    actions: np.ndarray


def controls_at(trial: int, regime: RecoveryRegime = REGIME) -> tuple[float, float]:
    """Return the scheduled precision and preference field at one trial."""
    if trial < regime.healthy_trials:
        return regime.gamma_healthy, regime.delta_c_healthy
    progress = min(
        max((trial - regime.healthy_trials) / regime.drift_trials, 0.0),
        1.0,
    )
    gamma = regime.gamma_healthy + progress * (
        regime.gamma_floor - regime.gamma_healthy
    )
    delta_c = regime.delta_c_healthy + progress * (
        regime.delta_c_floor - regime.delta_c_healthy
    )
    return float(gamma), float(delta_c)


def initial_state(regime: RecoveryRegime = REGIME) -> np.ndarray:
    """Return the four symmetric Dirichlet concentrations."""
    return np.repeat(regime.alpha_0 / 2.0, 4).astype(float)


def latent_probability(
    state: np.ndarray,
    gamma: float,
    delta_c: float,
) -> float:
    """Evaluate the engagement posterior from the learned state."""
    alpha_1, beta_1, alpha_2, beta_2 = state
    a = alpha_1 / (alpha_1 + beta_1)
    b = beta_2 / (alpha_2 + beta_2)
    delta_g = compute_efe_difference(a, b, delta_c)
    return float(policy_posterior_pi1(delta_g, gamma))


def advance_one(
    state: np.ndarray,
    rng: np.random.Generator,
    gamma: float,
    delta_c: float,
    eta_engaged: float = 1.0,
    eta_withdrawn: float = 1.0,
    regime: RecoveryRegime = REGIME,
) -> tuple[float, int]:
    """Advance one stochastic learning-action step in place."""
    prior = regime.alpha_0 / 2.0
    state[:2] = eta_engaged * state[:2] + (1.0 - eta_engaged) * prior
    state[2:] = eta_withdrawn * state[2:] + (1.0 - eta_withdrawn) * prior
    probability = latent_probability(state, gamma, delta_c)
    action = 0 if rng.random() < probability else 1
    if action == 0:
        state[0 if rng.random() < regime.p else 1] += 1.0
    else:
        state[2 if rng.random() < (1.0 - regime.p) else 3] += 1.0
    return probability, action


def simulate_to_intervention(
    seed: int,
    regime: RecoveryRegime = REGIME,
) -> Branch:
    """Run the common healthy, drift, and adverse-hold schedule."""
    rng = np.random.default_rng(seed)
    state = initial_state(regime)
    probabilities = np.zeros(regime.intervention_trial)
    actions = np.zeros(regime.intervention_trial, dtype=int)
    for trial in range(regime.intervention_trial):
        gamma, delta_c = controls_at(trial, regime)
        probabilities[trial], actions[trial] = advance_one(
            state,
            rng,
            gamma,
            delta_c,
            regime=regime,
        )
    return Branch(
        seed=seed,
        state=state.copy(),
        rng_state=copy.deepcopy(rng.bit_generator.state),
        probabilities=probabilities,
        actions=actions,
    )


def branch_is_verified(
    branch: Branch,
    withdrawn: bool,
    regime: RecoveryRegime = REGIME,
) -> bool:
    """Apply the prespecified 30-trial state-verification criterion."""
    recent = branch.probabilities[-regime.confirm_trials :]
    if withdrawn:
        return bool(branch.probabilities[-1] < 0.5 and np.mean(recent) < 0.5)
    return bool(branch.probabilities[-1] > 0.5 and np.mean(recent) > 0.5)


def select_branches(
    n_agents: int,
    withdrawn: bool,
    seed_start: int,
    regime: RecoveryRegime = REGIME,
) -> tuple[list[Branch], int]:
    """Select the first deterministic set meeting the state criterion."""
    selected: list[Branch] = []
    seed = seed_start
    while len(selected) < n_agents:
        branch = simulate_to_intervention(seed, regime)
        if branch_is_verified(branch, withdrawn, regime):
            selected.append(branch)
        seed += 1
        if seed - seed_start > 100000:
            raise RuntimeError("State-verification search exceeded 100000 candidates")
    return selected, seed - seed_start


def continue_branch(
    branch: Branch,
    gamma: float,
    delta_c: float,
    trials: int | None = None,
    late_window: int | None = None,
    eta_engaged: float = 1.0,
    eta_withdrawn: float = 1.0,
    regime: RecoveryRegime = REGIME,
    return_trace: bool = False,
):
    """Continue one paired branch under fixed controls."""
    n_trials = regime.recovery_trials if trials is None else trials
    window = regime.late_window if late_window is None else late_window
    state = branch.state.copy()
    rng = np.random.default_rng()
    rng.bit_generator.state = copy.deepcopy(branch.rng_state)
    probabilities = np.zeros(n_trials)
    actions = np.zeros(n_trials, dtype=int)
    for trial in range(n_trials):
        probabilities[trial], actions[trial] = advance_one(
            state,
            rng,
            gamma,
            delta_c,
            eta_engaged=eta_engaged,
            eta_withdrawn=eta_withdrawn,
            regime=regime,
        )
    result = float(np.mean(probabilities[-min(window, n_trials) :]))
    if return_trace:
        return result, probabilities, actions, state, copy.deepcopy(rng.bit_generator.state)
    return result


def adverse_residence_branch(
    branch: Branch,
    residence_trials: int,
    eta_engaged: float,
    eta_withdrawn: float,
    regime: RecoveryRegime = REGIME,
) -> Branch:
    """Hold a verified branch at adverse controls for a specified duration."""
    if residence_trials == 0:
        return copy.deepcopy(branch)
    _, probabilities, actions, state, rng_state = continue_branch(
        branch,
        regime.gamma_floor,
        regime.delta_c_floor,
        trials=residence_trials,
        late_window=min(residence_trials, regime.confirm_trials),
        eta_engaged=eta_engaged,
        eta_withdrawn=eta_withdrawn,
        regime=regime,
        return_trace=True,
    )
    return Branch(
        seed=branch.seed,
        state=state,
        rng_state=rng_state,
        probabilities=np.concatenate([branch.probabilities, probabilities]),
        actions=np.concatenate([branch.actions, actions]),
    )


def condition_spec(regime: RecoveryRegime = REGIME):
    """Return the four paired intervention conditions."""
    return [
        ("Control", regime.gamma_floor, regime.delta_c_floor),
        ("Restore gamma only", regime.gamma_healthy, regime.delta_c_floor),
        ("Restore delta_c only", regime.gamma_floor, regime.delta_c_healthy),
        ("Restore both", regime.gamma_healthy, regime.delta_c_healthy),
    ]


def write_rows(path: Path, rows: list[dict]) -> None:
    """Write a list of homogeneous records."""
    if not rows:
        raise ValueError(f"No rows supplied for {path}")
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def summarise(values: np.ndarray) -> dict[str, float]:
    """Return the preregistered distribution summaries."""
    values = np.asarray(values, dtype=float)
    return {
        "median": float(np.median(values)),
        "q1": float(np.percentile(values, 25)),
        "q3": float(np.percentile(values, 75)),
        "mean": float(np.mean(values)),
        "sd": float(np.std(values, ddof=1)),
        "proportion_above_0_5": float(np.mean(values > 0.5)),
    }


def save_figure(fig: plt.Figure, stem: str, manuscript_number: int) -> None:
    """Save the canonical vector and raster manuscript figures."""
    name = f"Fig_{manuscript_number}"
    fig.savefig(MANUSCRIPT_FIGURE_DIR / f"{name}.png", dpi=300)
    fig.savefig(MANUSCRIPT_FIGURE_DIR / f"{name}.pdf")


def intervention_fixed_points(regime: RecoveryRegime = REGIME) -> list[dict]:
    """Enumerate mean-field branches for each intervention condition."""
    tau = regime.intervention_trial / (2.0 * regime.alpha_0)
    rows: list[dict] = []
    for condition, gamma, delta_c in condition_spec(regime):
        for branch_index, (z_value, stable) in enumerate(
            find_fixed_points(tau, regime.p, gamma, delta_c),
            start=1,
        ):
            rows.append({
                "condition": condition,
                "tau": tau,
                "gamma": gamma,
                "delta_c": delta_c,
                "branch_index": branch_index,
                "z_fixed_point": float(z_value),
                "stable_flow": bool(stable),
            })
    return rows


def regime_audit(regime: RecoveryRegime = REGIME) -> list[dict]:
    """Audit coupling strength and fixed-point count along the live schedule."""
    rows: list[dict] = []
    horizon = regime.intervention_trial + regime.recovery_trials
    for accumulated in range(1, horizon + 1):
        gamma, delta_c = controls_at(accumulated - 1, regime)
        tau = accumulated / (2.0 * regime.alpha_0)
        fixed_points = find_fixed_points(tau, regime.p, gamma, delta_c)
        rows.append({
            "accumulated_actions": accumulated,
            "tau": tau,
            "gamma": gamma,
            "delta_c": delta_c,
            "gamma_times_coupling": gamma * coupling_function(tau, regime.p),
            "fixed_point_count": len(fixed_points),
            "stable_fixed_point_count": sum(stable for _, stable in fixed_points),
        })
    return rows


def structural_sensitivity(regime: RecoveryRegime = REGIME) -> list[dict]:
    """Check a small neighbouring grid without selecting on outcomes."""
    rows: list[dict] = []
    for alpha_0 in (35.0, 40.0, 45.0):
        for intervention_trial in (280, 300, 320):
            tau = intervention_trial / (2.0 * alpha_0)
            for gamma in (17.0, 18.0, 19.0):
                for delta_c in (-0.06, -0.05, -0.04):
                    fixed_points = find_fixed_points(
                        tau,
                        regime.p,
                        gamma,
                        delta_c,
                    )
                    rows.append({
                        "alpha_0": alpha_0,
                        "intervention_trial": intervention_trial,
                        "tau": tau,
                        "gamma": gamma,
                        "delta_c": delta_c,
                        "fixed_point_count": len(fixed_points),
                        "stable_fixed_point_count": sum(
                            stable for _, stable in fixed_points
                        ),
                        "structurally_bistable": int(
                            len(fixed_points) == 3
                            and sum(stable for _, stable in fixed_points) == 2
                        ),
                    })
    return rows


def run_primary_recovery(
    withdrawn_branches: list[Branch],
    candidates_examined: int,
    regime: RecoveryRegime = REGIME,
) -> tuple[list[dict], dict[str, np.ndarray]]:
    """Run the paired four-condition recovery comparison."""
    long_rows: list[dict] = []
    values_by_condition: dict[str, np.ndarray] = {}
    for condition, gamma, delta_c in condition_spec(regime):
        values = np.array([
            continue_branch(branch, gamma, delta_c, regime=regime)
            for branch in withdrawn_branches
        ])
        values_by_condition[condition] = values
        for branch, value in zip(withdrawn_branches, values):
            long_rows.append({
                "seed": branch.seed,
                "condition": condition,
                "gamma": gamma,
                "delta_c": delta_c,
                "late_p_engaged": value,
            })
    summary_rows: list[dict] = []
    for condition, gamma, delta_c in condition_spec(regime):
        summary_rows.append({
            "condition": condition,
            "n_verified_agents": len(withdrawn_branches),
            "candidate_schedules_examined": candidates_examined,
            "gamma": gamma,
            "delta_c": delta_c,
            **summarise(values_by_condition[condition]),
        })
    write_rows(TABLE_DIR / "bistable_recovery_long.csv", long_rows)
    write_rows(TABLE_DIR / "bistable_recovery_summary.csv", summary_rows)
    return summary_rows, values_by_condition


def make_figure_3(
    withdrawn_branches: list[Branch],
    summary_rows: list[dict],
    values_by_condition: dict[str, np.ndarray],
    regime: RecoveryRegime = REGIME,
) -> None:
    """Create the state-anchored intervention comparison."""
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.8))
    order = [item[0] for item in condition_spec(regime)]
    labels = ["Control", r"$\gamma$ only", r"$\Delta c$ only", "Both"]
    colours = ["#8A9499", "#4C78A8", "#F28E2B", "#2E7D32"]
    ax = axes[0]
    for index, (row, colour) in enumerate(zip(summary_rows, colours)):
        ax.vlines(index, row["q1"], row["q3"], color="#27343A", lw=2.0)
        ax.scatter(
            index,
            row["median"],
            s=72,
            color=colour,
            edgecolor="#27343A",
            linewidth=0.8,
            zorder=3,
        )
    ax.axhline(0.5, color="#9AA3A7", ls="--", lw=1.0)
    ax.set_xticks(range(4), labels)
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("Late-window P(engaged)")
    ax.set_title("(a)  Recovery from verified withdrawal")

    control_values = values_by_condition["Control"]
    example_index = int(np.argmin(np.abs(control_values - np.median(control_values))))
    branch = withdrawn_branches[example_index]
    ax = axes[1]
    pre_trials = np.arange(regime.intervention_trial)
    ax.plot(
        pre_trials,
        branch.probabilities,
        color="#AAB2B6",
        lw=1.4,
        label="Shared pre-intervention path",
    )
    post_trials = np.arange(
        regime.intervention_trial,
        regime.intervention_trial + regime.recovery_trials,
    )
    for (condition, gamma, delta_c), colour in zip(
        condition_spec(regime),
        colours,
    ):
        _, trace, _, _, _ = continue_branch(
            branch,
            gamma,
            delta_c,
            regime=regime,
            return_trace=True,
        )
        display_label = {
            "Restore gamma only": r"Restore $\gamma$ only",
            "Restore delta_c only": r"Restore $\Delta c$ only",
        }.get(condition, condition)
        ax.plot(post_trials, trace, color=colour, lw=1.5, label=display_label)
    ax.axvline(regime.intervention_trial, color="#176B87", ls="--", lw=1.2)
    ax.axhline(0.5, color="#9AA3A7", ls=":", lw=0.9)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Trial")
    ax.set_ylabel("Latent P(engaged)")
    ax.set_title("(b)  Paired trajectories from one learned state")
    ax.legend(frameon=False, fontsize=6.8, loc="center left", bbox_to_anchor=(1.0, 0.5))
    plt.tight_layout()
    save_figure(fig, "fig_S2_bistable_recovery", 3)
    plt.close(fig)


def mean_field_grid(regime: RecoveryRegime = REGIME):
    """Evaluate fixed-point multiplicity across precision-field space."""
    gamma_values = np.linspace(12.0, 24.0, 97)
    delta_values = np.linspace(-0.16, 0.16, 129)
    tau = regime.intervention_trial / (2.0 * regime.alpha_0)
    count = np.zeros((len(delta_values), len(gamma_values)))
    for row, delta_c in enumerate(delta_values):
        for column, gamma in enumerate(gamma_values):
            count[row, column] = len(
                find_fixed_points(tau, regime.p, gamma, delta_c)
            )
    return gamma_values, delta_values, count


def finite_recovery_grid(
    withdrawn_branches: list[Branch],
    regime: RecoveryRegime = REGIME,
) -> tuple[list[dict], np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate late recovery from identical withdrawn histories."""
    gamma_values = np.arange(14.0, 22.01, 1.0)
    delta_values = np.arange(-0.10, 0.201, 0.025)
    subset = withdrawn_branches[:60]
    grid = np.zeros((len(delta_values), len(gamma_values)))
    rows: list[dict] = []
    for row, delta_c in enumerate(delta_values):
        for column, gamma in enumerate(gamma_values):
            values = np.array([
                continue_branch(branch, gamma, float(delta_c), regime=regime)
                for branch in subset
            ])
            grid[row, column] = np.median(values)
            rows.append({
                "gamma": gamma,
                "delta_c": float(delta_c),
                "n_verified_agents": len(subset),
                **summarise(values),
            })
    write_rows(TABLE_DIR / "bistable_recovery_grid.csv", rows)
    return rows, gamma_values, delta_values, grid


def make_figure_4(
    withdrawn_branches: list[Branch],
    audit_rows: list[dict],
    regime: RecoveryRegime = REGIME,
) -> None:
    """Create the fold map, finite-agent boundary, and regime audit."""
    gamma_values, delta_values, fixed_count = mean_field_grid(regime)
    _, finite_gamma, finite_delta, finite_grid = finite_recovery_grid(
        withdrawn_branches,
        regime,
    )
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.75))
    ax = axes[0]
    ax.contourf(
        gamma_values,
        delta_values,
        fixed_count == 3,
        levels=[-0.5, 0.5, 1.5],
        colors=["#F3F5F6", "#DCECEF"],
    )
    ax.contour(
        gamma_values,
        delta_values,
        fixed_count == 3,
        levels=[0.5],
        colors=["#176B87"],
        linewidths=1.6,
    )
    markers = ["o", "s", "D", "^"]
    colours = ["#8A9499", "#4C78A8", "#F28E2B", "#2E7D32"]
    for (condition, gamma, delta_c), marker, colour in zip(
        condition_spec(regime),
        markers,
        colours,
    ):
        ax.scatter(gamma, delta_c, marker=marker, s=45, color=colour, edgecolor="#27343A", label=condition)
    ax.set_xlabel(r"Policy precision, $\gamma$")
    ax.set_ylabel(r"Preference field, $\Delta c$")
    ax.set_title("(a)  Mean-field fold region at intervention")
    ax.legend(frameon=False, fontsize=6.4, loc="lower right")

    ax = axes[1]
    image = ax.imshow(
        finite_grid,
        origin="lower",
        aspect="auto",
        extent=[
            finite_gamma.min() - 0.5,
            finite_gamma.max() + 0.5,
            finite_delta.min() - 0.0125,
            finite_delta.max() + 0.0125,
        ],
        cmap="viridis",
        vmin=0,
        vmax=1,
    )
    ax.contour(
        finite_gamma,
        finite_delta,
        finite_grid,
        levels=[0.5],
        colors=["white"],
        linewidths=1.5,
    )
    ax.set_xlabel(r"Restored $\gamma$")
    ax.set_ylabel(r"Restored $\Delta c$")
    ax.set_title("(b)  State-anchored finite-agent recovery")
    colourbar = fig.colorbar(image, ax=ax, fraction=0.048, pad=0.04)
    colourbar.set_label("Median late P(engaged)", fontsize=8)
    colourbar.ax.tick_params(labelsize=7)

    ax = axes[2]
    accumulated = np.array([row["accumulated_actions"] for row in audit_rows])
    ratio = np.array([row["gamma_times_coupling"] for row in audit_rows])
    multiplicity = np.array([row["fixed_point_count"] for row in audit_rows])
    ax.plot(accumulated, ratio, color="#176B87", lw=1.8, label=r"$\gamma\mathcal{G}(\tau;p)$")
    ax.axhline(1.0, color="#27343A", ls="--", lw=1.0, label="Symmetric critical value")
    ax.axvline(regime.intervention_trial, color="#7B3294", ls=":", lw=1.4, label="Intervention")
    mask = multiplicity == 3
    ax.fill_between(
        accumulated,
        0,
        np.maximum(ratio.max() * 1.05, 1.1),
        where=mask,
        color="#DCECEF",
        alpha=0.75,
        label="Three fixed points on live schedule",
    )
    ax.set_xlim(0, accumulated.max())
    ax.set_ylim(0, max(ratio.max() * 1.05, 1.1))
    ax.set_xlabel("Accumulated actions")
    ax.set_ylabel(r"$\gamma\mathcal{G}(\tau;p)$")
    ax.set_title("(c)  Regime audit across the live schedule")
    ax.legend(
        frameon=False,
        fontsize=6.2,
        loc="center right",
        bbox_to_anchor=(1.01, 0.35),
    )
    plt.tight_layout()
    save_figure(fig, "fig_S2_bistable_boundary", 4)
    plt.close(fig)


def history_response_curves(
    withdrawn_branches: list[Branch],
    engaged_branches: list[Branch],
    regime: RecoveryRegime = REGIME,
) -> list[dict]:
    """Compare identical controls after engaged versus withdrawn histories."""
    rows: list[dict] = []
    delta_grid = np.arange(-0.10, 0.151, 0.02)
    cohorts = [
        ("Verified engaged history", engaged_branches[:100]),
        ("Verified withdrawn history", withdrawn_branches[:100]),
    ]
    for history, branches in cohorts:
        for delta_c in delta_grid:
            values = np.array([
                continue_branch(
                    branch,
                    regime.gamma_floor,
                    float(delta_c),
                    regime=regime,
                )
                for branch in branches
            ])
            rows.append({
                "history": history,
                "gamma": regime.gamma_floor,
                "delta_c": float(delta_c),
                "n_agents": len(branches),
                **summarise(values),
            })
    write_rows(TABLE_DIR / "bistable_history_response.csv", rows)
    return rows


def threshold_from_matrix(matrix: np.ndarray, delta_grid: np.ndarray) -> float:
    """Return the first field where the sample median reaches .50."""
    medians = np.median(matrix, axis=0)
    indices = np.flatnonzero(medians >= 0.5)
    return float(delta_grid[indices[0]]) if len(indices) else float("nan")


def retention_thresholds(
    withdrawn_branches: list[Branch],
    regime: RecoveryRegime = REGIME,
) -> list[dict]:
    """Estimate state-anchored recovery thresholds with agent bootstrap CIs."""
    retention_conditions = [
        ("Permanent evidence", 1.0, 1.0),
        ("Symmetric decay", 0.995, 0.995),
        ("Asymmetric mild", 0.990, 0.997),
        ("Asymmetric strong", 0.985, 0.997),
    ]
    residence_values = (0, 50, 100, 150)
    delta_grid = np.arange(-0.02, 0.3001, 0.01)
    branches = withdrawn_branches[:100]
    rows: list[dict] = []
    for condition_index, (condition, eta_engaged, eta_withdrawn) in enumerate(
        retention_conditions
    ):
        for residence in residence_values:
            resident = [
                adverse_residence_branch(
                    branch,
                    residence,
                    eta_engaged,
                    eta_withdrawn,
                    regime,
                )
                for branch in branches
            ]
            matrix = np.zeros((len(resident), len(delta_grid)))
            for column, delta_c in enumerate(delta_grid):
                matrix[:, column] = [
                    continue_branch(
                        branch,
                        regime.gamma_healthy,
                        float(delta_c),
                        eta_engaged=eta_engaged,
                        eta_withdrawn=eta_withdrawn,
                        regime=regime,
                    )
                    for branch in resident
                ]
            estimate = threshold_from_matrix(matrix, delta_grid)
            rng = np.random.default_rng(710000 + condition_index * 1000 + residence)
            bootstrap = np.zeros(1000)
            for iteration in range(len(bootstrap)):
                sampled = rng.integers(0, len(resident), len(resident))
                bootstrap[iteration] = threshold_from_matrix(
                    matrix[sampled, :],
                    delta_grid,
                )
            finite = bootstrap[np.isfinite(bootstrap)]
            ci_low = (
                float(np.percentile(finite, 2.5))
                if len(finite)
                else float("nan")
            )
            ci_high = (
                float(np.percentile(finite, 97.5))
                if len(finite)
                else float("nan")
            )
            rows.append({
                "condition": condition,
                "eta_engaged": eta_engaged,
                "eta_withdrawn": eta_withdrawn,
                "post_verification_adverse_trials": residence,
                "n_agents": len(resident),
                "delta_grid_step": 0.01,
                "min_restored_delta_c": estimate,
                "bootstrap_ci_low": ci_low,
                "bootstrap_ci_high": ci_high,
            })
    write_rows(TABLE_DIR / "bistable_retention_thresholds.csv", rows)
    return rows


def make_figure_5(
    history_rows: list[dict],
    threshold_rows: list[dict],
) -> None:
    """Create the path-dependence and retention-sensitivity figure."""
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.75))
    ax = axes[0]
    for history, colour in [
        ("Verified engaged history", "#2E7D32"),
        ("Verified withdrawn history", "#7B3294"),
    ]:
        selected = [row for row in history_rows if row["history"] == history]
        x = np.array([row["delta_c"] for row in selected])
        median = np.array([row["median"] for row in selected])
        q1 = np.array([row["q1"] for row in selected])
        q3 = np.array([row["q3"] for row in selected])
        ax.plot(x, median, marker="o", ms=3.5, lw=1.7, color=colour, label=history)
        ax.fill_between(x, q1, q3, color=colour, alpha=0.16)
    ax.axhline(0.5, color="#9AA3A7", ls="--", lw=1.0)
    ax.axvline(0.0, color="#9AA3A7", ls=":", lw=0.8)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel(r"Post-intervention $\Delta c$ at $\gamma=18$")
    ax.set_ylabel("Late-window P(engaged)")
    ax.set_title("(a)  The same controls yield history-dependent states")
    ax.legend(frameon=False, fontsize=7)

    ax = axes[1]
    style = [
        ("Permanent evidence", "#8A9499", "o"),
        ("Symmetric decay", "#2CA25F", "s"),
        ("Asymmetric mild", "#F28E2B", "D"),
        ("Asymmetric strong", "#7B3294", "^"),
    ]
    for condition, colour, marker in style:
        selected = [row for row in threshold_rows if row["condition"] == condition]
        x = np.array([row["post_verification_adverse_trials"] for row in selected])
        y = np.array([row["min_restored_delta_c"] for row in selected])
        finite = np.isfinite(y)
        lower = y[finite] - np.array(
            [row["bootstrap_ci_low"] for row in selected]
        )[finite]
        upper = np.array(
            [row["bootstrap_ci_high"] for row in selected]
        )[finite] - y[finite]
        ax.errorbar(
            x[finite],
            y[finite],
            yerr=[lower, upper],
            marker=marker,
            ms=4.5,
            lw=1.5,
            capsize=3,
            color=colour,
            label=condition,
        )
        if np.any(~finite):
            ax.scatter(
                x[~finite],
                np.full(np.sum(~finite), 0.30),
                marker=marker,
                s=32,
                facecolors="white",
                edgecolors=colour,
                linewidths=1.2,
                zorder=4,
            )
            for missing_x in x[~finite]:
                ax.annotate(
                    r"$>.30$",
                    (missing_x, 0.30),
                    xytext=(0, 7),
                    textcoords="offset points",
                    ha="center",
                    fontsize=6.3,
                    color=colour,
                )
    ax.set_xlabel("Additional adverse trials after verified withdrawal")
    ax.set_ylabel(r"Minimum restored $\Delta c$")
    ax.set_ylim(-0.005, 0.335)
    ax.set_title("(b)  Recovery thresholds with bootstrap intervals")
    ax.legend(frameon=False, fontsize=6.7)
    plt.tight_layout()
    save_figure(fig, "fig_S5_bistable_hysteresis", 5)
    plt.close(fig)


def write_verified_cohort(
    withdrawn_branches: list[Branch],
    candidates_examined: int,
    regime: RecoveryRegime = REGIME,
) -> None:
    """Record the exact state-selection diagnostics."""
    rows: list[dict] = []
    for branch in withdrawn_branches:
        recent = branch.probabilities[-regime.confirm_trials :]
        cumulative_z = 1.0 - 2.0 * float(np.mean(branch.actions))
        recent_z = 1.0 - 2.0 * float(np.mean(branch.actions[-regime.confirm_trials :]))
        rows.append({
            "seed": branch.seed,
            "intervention_trial": regime.intervention_trial,
            "confirm_trials": regime.confirm_trials,
            "p_engaged_at_intervention": float(branch.probabilities[-1]),
            "mean_p_engaged_prior_confirmation_window": float(np.mean(recent)),
            "cumulative_action_order_parameter_z": cumulative_z,
            "recent_action_order_parameter_z": recent_z,
            "candidates_examined_for_200_verified_agents": candidates_examined,
        })
    write_rows(TABLE_DIR / "bistable_verified_withdrawal_cohort.csv", rows)


def run_all() -> dict:
    """Run the complete bistable recovery analysis and save canonical outputs."""
    withdrawn, candidates_examined = select_branches(
        200,
        withdrawn=True,
        seed_start=10000,
    )
    engaged, _ = select_branches(
        200,
        withdrawn=False,
        seed_start=50000,
    )
    write_verified_cohort(withdrawn, candidates_examined)
    fixed_point_rows = intervention_fixed_points()
    audit_rows = regime_audit()
    sensitivity_rows = structural_sensitivity()
    write_rows(TABLE_DIR / "bistable_intervention_fixed_points.csv", fixed_point_rows)
    write_rows(TABLE_DIR / "bistable_regime_audit.csv", audit_rows)
    write_rows(TABLE_DIR / "bistable_structural_sensitivity.csv", sensitivity_rows)
    summary_rows, values_by_condition = run_primary_recovery(
        withdrawn,
        candidates_examined,
    )
    make_figure_3(withdrawn, summary_rows, values_by_condition)
    make_figure_4(withdrawn, audit_rows)
    history_rows = history_response_curves(withdrawn, engaged)
    threshold_rows = retention_thresholds(withdrawn)
    make_figure_5(history_rows, threshold_rows)
    print("Bistable recovery outputs written to", TABLE_DIR)
    return {
        "summary": summary_rows,
        "fixed_points": fixed_point_rows,
        "regime_audit": audit_rows,
        "history_response": history_rows,
        "retention_thresholds": threshold_rows,
    }


if __name__ == "__main__":
    run_all()
