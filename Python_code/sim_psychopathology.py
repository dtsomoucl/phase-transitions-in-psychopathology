"""Finite-agent onset simulation and quasistatic regime audit."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np

from core_functions import compute_efe_difference, find_fixed_points, policy_posterior_pi1
from psychopathology_regimes import onset_regime


ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = ROOT / "outputs" / "manuscript_figures"
TABLE_DIR = ROOT / "outputs" / "tables"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)

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


def causal_mean(values: np.ndarray, width: int = 15) -> np.ndarray:
    """Return a trailing mean without padding across the left boundary."""
    values = np.asarray(values, dtype=float)
    cumulative = np.cumsum(values)
    output = np.empty(len(values), dtype=float)
    for index in range(len(values)):
        start = max(0, index - width + 1)
        total = cumulative[index] - (cumulative[start - 1] if start else 0.0)
        output[index] = total / (index - start + 1)
    return output


def run_agent(
    seed: int,
    n_total: int = 600,
    regime: dict | None = None,
) -> dict[str, np.ndarray]:
    """Simulate one agent under the prespecified gradual-adversity schedule."""
    parameters = onset_regime() if regime is None else regime.copy()
    rng = np.random.default_rng(seed)
    alpha_0 = parameters["alpha_0"]
    alpha_1 = beta_1 = alpha_0 / 2.0
    alpha_2 = beta_2 = alpha_0 / 2.0
    n_1 = n_2 = 0
    histories = {
        key: np.zeros(n_total, dtype=float)
        for key in ("probability", "gamma", "delta_c", "z", "n_1", "n_2")
    }

    for trial in range(n_total):
        if trial < parameters["N_healthy"]:
            gamma = parameters["gamma_healthy"]
            delta_c = parameters["delta_c_healthy"]
        else:
            elapsed = trial - parameters["N_healthy"]
            gamma = max(
                parameters["gamma_healthy"] - parameters["gamma_rate"] * elapsed,
                parameters["gamma_floor"],
            )
            delta_c = max(
                parameters["delta_c_healthy"] - parameters["delta_c_rate"] * elapsed,
                parameters["delta_c_floor"],
            )

        a = alpha_1 / (alpha_1 + beta_1)
        b = beta_2 / (alpha_2 + beta_2)
        delta_g = compute_efe_difference(a, b, delta_c)
        probability = policy_posterior_pi1(delta_g, gamma)
        action = 0 if rng.random() < probability else 1

        if action == 0:
            engaged_outcome = rng.random() < parameters["p"]
            alpha_1 += float(engaged_outcome)
            beta_1 += float(not engaged_outcome)
            n_1 += 1
        else:
            engaged_outcome = rng.random() < (1.0 - parameters["p"])
            alpha_2 += float(engaged_outcome)
            beta_2 += float(not engaged_outcome)
            n_2 += 1

        histories["probability"][trial] = probability
        histories["gamma"][trial] = gamma
        histories["delta_c"][trial] = delta_c
        histories["z"][trial] = (n_1 - n_2) / (n_1 + n_2)
        histories["n_1"][trial] = n_1
        histories["n_2"][trial] = n_2
    return histories


def quasistatic_structure(
    history: dict[str, np.ndarray],
    every: int = 2,
    regime: dict | None = None,
) -> dict[str, np.ndarray]:
    """Enumerate mean-field fixed points along the live finite-agent schedule."""
    parameters = onset_regime() if regime is None else regime.copy()
    trials = np.arange(0, len(history["probability"]), every)
    result = {
        "trial": trials,
        "count": np.zeros(len(trials), dtype=int),
        "engaged": np.full(len(trials), np.nan),
        "withdrawn": np.full(len(trials), np.nan),
        "unstable": np.full(len(trials), np.nan),
    }
    for index, trial in enumerate(trials):
        accumulated = history["n_1"][trial] + history["n_2"][trial]
        tau = accumulated / (2.0 * parameters["alpha_0"])
        points = find_fixed_points(
            tau,
            parameters["p"],
            history["gamma"][trial],
            history["delta_c"][trial],
            z_grid=np.linspace(-0.999, 0.999, 800),
        )
        stable = sorted(point for point, is_stable in points if is_stable)
        unstable = [point for point, is_stable in points if not is_stable]
        result["count"][index] = len(points)
        if stable:
            result["engaged"][index] = stable[-1]
        if len(stable) > 1:
            result["withdrawn"][index] = stable[0]
        if unstable:
            result["unstable"][index] = unstable[0]
    return result


def fold_loss_trial(structure: dict[str, np.ndarray]) -> int | None:
    """Return the first sampled trial after loss of the engaged stable branch."""
    for index in range(1, len(structure["trial"])):
        if (
            structure["count"][index - 1] == 3
            and structure["count"][index] == 1
            and structure["engaged"][index] < 0
        ):
            return int(structure["trial"][index])
    return None


def transition_trial(probability: np.ndarray, width: int = 15) -> int | None:
    """Locate sustained withdrawal using the prespecified trailing-mean rule."""
    smoothed = causal_mean(probability, width)
    for trial in range(width, len(smoothed) - width + 1):
        if smoothed[trial] < 0.5 and np.all(smoothed[trial:trial + width] < 0.5):
            return trial
    return None


def high_to_low_duration(probability: np.ndarray, width: int = 10) -> float:
    """Measure the last .80-to-first-.20 transition duration when observable."""
    smoothed = causal_mean(probability, width)
    low = np.flatnonzero(smoothed <= 0.20)
    if not len(low):
        return float("nan")
    first_low = int(low[0])
    high = np.flatnonzero((np.arange(len(smoothed)) < first_low) & (smoothed >= 0.80))
    return float(first_low - high[-1]) if len(high) else float("nan")


def onset_diagnostics(
    histories: list[dict[str, np.ndarray]],
    fold_trial: int,
) -> tuple[list[dict], dict]:
    """Write operational timing and abruptness diagnostics for every agent."""
    rows = []
    for index, history in enumerate(histories):
        transition = transition_trial(history["probability"])
        probability = history["probability"]
        rows.append({
            "seed": 500 + index,
            "fold_loss_trial": fold_trial,
            "transition_trial": transition if transition is not None else "",
            "fold_to_transition_lag": transition - fold_trial if transition is not None else "",
            "maximum_one_trial_decline": float(np.max(probability[:-1] - probability[1:])),
            "maximum_ten_trial_decline": float(np.max(probability[:-10] - probability[10:])),
            "smoothed_80_to_20_duration": high_to_low_duration(probability),
        })
    with (TABLE_DIR / "onset_agent_diagnostics.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    transitions = np.array([
        row["transition_trial"] for row in rows if row["transition_trial"] != ""
    ], dtype=float)
    durations = np.array([
        row["smoothed_80_to_20_duration"] for row in rows
        if np.isfinite(row["smoothed_80_to_20_duration"])
    ])
    ten_trial = np.array([row["maximum_ten_trial_decline"] for row in rows])
    summary = {
        "n_agents": len(rows),
        "n_transitioned_by_trial_600": len(transitions),
        "fold_loss_trial": fold_trial,
        "transition_trial_median": float(np.median(transitions)),
        "transition_trial_q1": float(np.percentile(transitions, 25)),
        "transition_trial_q3": float(np.percentile(transitions, 75)),
        "fold_to_transition_lag_median": float(np.median(transitions - fold_trial)),
        "fold_to_transition_lag_q1": float(np.percentile(transitions - fold_trial, 25)),
        "fold_to_transition_lag_q3": float(np.percentile(transitions - fold_trial, 75)),
        "maximum_ten_trial_decline_median": float(np.median(ten_trial)),
        "n_with_observed_80_to_20_transition": len(durations),
        "smoothed_80_to_20_duration_median": float(np.median(durations)),
        "smoothed_80_to_20_duration_q1": float(np.percentile(durations, 25)),
        "smoothed_80_to_20_duration_q3": float(np.percentile(durations, 75)),
        "duration_as_fraction_of_600_trial_window": float(np.median(durations) / 600.0),
    }
    with (TABLE_DIR / "onset_diagnostics_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary.keys())
        writer.writeheader()
        writer.writerow(summary)
    return rows, summary


def make_figure() -> dict:
    """Generate the six-panel onset figure and its numerical audit tables."""
    parameters = onset_regime()
    histories = [run_agent(500 + index, regime=parameters) for index in range(200)]
    analysis_histories = histories[:100]
    representative = analysis_histories[0]
    structure = quasistatic_structure(representative, regime=parameters)
    fold_trial = fold_loss_trial(structure)
    if fold_trial is None:
        raise RuntimeError("The prespecified onset schedule did not lose its engaged branch.")
    rows, summary = onset_diagnostics(analysis_histories, fold_trial)
    target = summary["transition_trial_median"]
    display_index = min(
        range(len(analysis_histories)),
        key=lambda index: abs(
            (rows[index]["transition_trial"] if rows[index]["transition_trial"] != "" else 600)
            - target
        ),
    )
    display = analysis_histories[display_index]
    display_transition = int(rows[display_index]["transition_trial"])

    figure = plt.figure(figsize=(11.2, 8.1))
    grid = GridSpec(3, 2, figure=figure, hspace=0.50, wspace=0.35)
    axes = [figure.add_subplot(grid[row, column]) for row in range(3) for column in range(2)]
    ax_a, ax_b, ax_c, ax_d, ax_e, ax_f = axes
    trials = np.arange(600)

    for index in np.linspace(0, 99, 25, dtype=int):
        ax_a.plot(trials, causal_mean(analysis_histories[index]["probability"]), color="#9FB3C8", alpha=0.20, lw=0.8)
    population_mean = np.mean([history["probability"] for history in analysis_histories], axis=0)
    ax_a.plot(trials, causal_mean(population_mean), color="#27343A", lw=2.2, label="Population mean")
    ax_a.axvline(parameters["N_healthy"], color="#B23A48", ls="--", lw=1.1, label="Adversity onset")
    ax_a.axvline(fold_trial, color="#7B3294", ls=":", lw=1.4, label=f"Fold loss ({fold_trial})")
    ax_a.set(xlabel="Trial", ylabel="Latent P(engaged)", ylim=(-0.02, 1.02))
    ax_a.set_title("(a)  Ensemble trajectories")
    ax_a.legend(frameon=False, fontsize=6.8, loc="lower left")

    ax_b.plot(trials, causal_mean(display["probability"], 10), color="#4C78A8", lw=1.8)
    ax_b.axvline(fold_trial, color="#7B3294", ls=":", lw=1.4, label=f"Fold loss ({fold_trial})")
    ax_b.axvline(display_transition, color="#F28E2B", ls="--", lw=1.3, label=f"Sustained crossing ({display_transition})")
    ax_b.axhline(0.5, color="#9AA3A7", ls=":", lw=0.9)
    ax_b.set(xlabel="Trial", ylabel="Latent P(engaged)", ylim=(-0.02, 1.02))
    ax_b.set_title("(b)  Illustrative finite-agent trajectory")
    ax_b.legend(frameon=False, fontsize=6.8, loc="lower left")

    engaged = (1.0 + structure["engaged"]) / 2.0
    withdrawn = (1.0 + structure["withdrawn"]) / 2.0
    unstable = (1.0 + structure["unstable"]) / 2.0
    ax_c.plot(structure["trial"], engaged, color="#2E7D32", lw=1.8, label="Stable engaged branch")
    ax_c.plot(structure["trial"], withdrawn, color="#B23A48", lw=1.8, label="Stable withdrawn branch")
    ax_c.plot(structure["trial"], unstable, color="#27343A", ls="--", lw=1.0, label="Unstable branch")
    ax_c.plot(trials, causal_mean(display["probability"], 10), color="#4C78A8", lw=1.2, label="Finite-agent trajectory")
    ax_c.axvline(fold_trial, color="#7B3294", ls=":", lw=1.4)
    ax_c.set(xlabel="Trial", ylabel="Branch / latent P(engaged)", ylim=(-0.02, 1.02))
    ax_c.set_title("(c)  Quasistatic structure and kinetic lag")
    ax_c.legend(frameon=False, fontsize=6.1, loc="lower right")

    ax_d.plot(trials, display["gamma"], color="#B23A48", lw=1.8)
    ax_d.axvline(parameters["N_healthy"], color="#9AA3A7", ls="--", lw=1.0)
    ax_d.axvline(fold_trial, color="#7B3294", ls=":", lw=1.4)
    ax_d.set(xlabel="Trial", ylabel=r"Policy precision, $\gamma$")
    ax_d.set_title("(d)  Precision schedule")

    ax_e.plot(trials, display["delta_c"], color="#F28E2B", lw=1.8)
    ax_e.axhline(0, color="#9AA3A7", ls="--", lw=0.9)
    ax_e.axvline(parameters["N_healthy"], color="#9AA3A7", ls="--", lw=1.0)
    ax_e.axvline(fold_trial, color="#7B3294", ls=":", lw=1.4)
    ax_e.set(xlabel="Trial", ylabel=r"Preference field, $\Delta c$")
    ax_e.set_title("(e)  Preference-field schedule")

    all_probabilities = np.array([history["probability"] for history in histories])
    variance = np.var(all_probabilities, axis=0)
    mean = np.mean(all_probabilities, axis=0)
    variance_axis = ax_f
    mean_axis = ax_f.twinx()
    variance_axis.plot(trials, causal_mean(variance), color="#B23A48", lw=1.8, label="Variance")
    mean_axis.plot(trials, causal_mean(mean), color="#4C78A8", ls="--", lw=1.6, label="Mean")
    peak_trial = int(parameters["N_healthy"] + np.argmax(causal_mean(variance)[parameters["N_healthy"]:]))
    variance_axis.axvline(fold_trial, color="#7B3294", ls=":", lw=1.4, label="Fold loss")
    variance_axis.axvline(peak_trial, color="#F28E2B", ls="--", lw=1.3, label=f"Variance peak ({peak_trial})")
    variance_axis.set(xlabel="Trial", ylabel="Between-agent variance")
    mean_axis.set_ylabel("Population mean")
    variance_axis.set_title("(f)  Ensemble variance and mean")
    variance_axis.legend(frameon=False, fontsize=6.4, loc="upper left")
    mean_axis.legend(frameon=False, fontsize=6.4, loc="upper right")

    figure.subplots_adjust(left=0.08, right=0.92, bottom=0.08, top=0.95, hspace=0.50, wspace=0.38)
    figure.savefig(FIGURE_DIR / "Fig_2.png")
    figure.savefig(FIGURE_DIR / "Fig_2.pdf")
    plt.close(figure)
    summary["variance_peak_trial"] = peak_trial
    return summary


def run_all() -> dict:
    """Run the onset analysis used by the manuscript and SOM."""
    summary = make_figure()
    print(
        "Onset analysis complete: "
        f"fold={summary['fold_loss_trial']}, "
        f"transition median={summary['transition_trial_median']:.0f}."
    )
    return summary


if __name__ == "__main__":
    run_all()
