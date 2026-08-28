"""Clinical trajectory stress test for the active-inference transition model.

This module analyses the nine 209-week psychosis, depression, and mania
trajectories released with Jin et al. (2023).  It deliberately treats the
application as a small-sample stress test rather than validation.  The main
comparison uses leave-one-patient-out one-step forecasts from an intercept
baseline, a symptom-wise first-order Markov model, a linear burden model, a
two-state Bernoulli hidden Markov model, and a nonlinear mean-field map
derived from the active-inference/Ising reduction.

The supplied Lorenz DCM estimates are audited descriptively.  They are not
entered into the held-out comparison because the released SPM inversion is
an in-sample variational fit and does not provide held-out forecasts under a
common observation model.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat
from scipy.optimize import minimize
from scipy.special import expit


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_DIR = ROOT / "VALIDATION" / "data" / "jin_2023"
DEFAULT_OUTPUT_DIR = ROOT / "outputs" / "clinical_trajectories"
SYMPTOM_NAMES = ("Psychosis", "Depression", "Mania")
LORENZ_CHAOS_THRESHOLD = 24.7368
EPS = 1e-9


@dataclass(frozen=True)
class Prediction:
    model: str
    subject: int
    week: int
    symptom: str
    actual: float
    probability: float
    changed: bool


def _as_subject_list(raw: object) -> list[np.ndarray]:
    """Return the MATLAB cell/array representation as one array per patient."""
    values = np.asarray(raw)
    if values.dtype != object and values.ndim == 3:
        return [np.asarray(values[index], dtype=float) for index in range(values.shape[0])]
    return [np.asarray(item, dtype=float) for item in np.atleast_1d(raw)]


def load_trajectories(path: Path) -> tuple[list[np.ndarray], dict[str, str]]:
    """Load the released symptom trajectories and their provenance log."""
    raw = loadmat(path, simplify_cells=True)
    data = raw["D"]
    metadata = {key: str(value) for key, value in raw.get("log", {}).items()}
    return _as_subject_list(data["d"]), metadata


def _finite_pairs(sequences: Sequence[np.ndarray]) -> Iterable[tuple[int, int, np.ndarray, np.ndarray]]:
    """Yield adjacent observations with at least one finite current and future value."""
    for subject, values in enumerate(sequences, start=1):
        for week in range(values.shape[0] - 1):
            current = values[week]
            following = values[week + 1]
            if np.any(np.isfinite(current)) and np.any(np.isfinite(following)):
                yield subject, week + 1, current, following


def describe_sequences(sequences: Sequence[np.ndarray]) -> list[dict[str, object]]:
    """Create participant-level prevalence, missingness, and transition summaries."""
    rows: list[dict[str, object]] = []
    for subject, values in enumerate(sequences, start=1):
        row: dict[str, object] = {"subject": subject, "weeks": values.shape[0]}
        total_pairs = 0
        total_changes = 0
        for column, name in enumerate(SYMPTOM_NAMES):
            series = values[:, column]
            previous = series[:-1]
            following = series[1:]
            mask = np.isfinite(previous) & np.isfinite(following)
            changes = int(np.sum(previous[mask] != following[mask]))
            row[f"{name.lower()}_prevalence"] = float(np.nanmean(series))
            row[f"{name.lower()}_missing"] = int(np.isnan(series).sum())
            row[f"{name.lower()}_transitions"] = changes
            total_pairs += int(mask.sum())
            total_changes += changes
        row["adjacent_pairs"] = total_pairs
        row["transitions"] = total_changes
        row["persistence_accuracy"] = 1.0 - total_changes / total_pairs if total_pairs else math.nan
        rows.append(row)
    return rows


def _burden(values: np.ndarray) -> float:
    """Map observed binary symptoms to the order-parameter scale [-1, 1]."""
    finite = values[np.isfinite(values)]
    return float(2.0 * np.mean(finite) - 1.0) if finite.size else 0.0


def _training_pairs(sequences: Sequence[np.ndarray]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return current burden, next symptom, and symptom index for finite targets."""
    burdens: list[float] = []
    outcomes: list[float] = []
    symptoms: list[int] = []
    for _, _, current, following in _finite_pairs(sequences):
        z_value = _burden(current)
        for symptom, target in enumerate(following):
            if np.isfinite(target):
                burdens.append(z_value)
                outcomes.append(float(target))
                symptoms.append(symptom)
    return np.asarray(burdens), np.asarray(outcomes), np.asarray(symptoms, dtype=int)


def _bernoulli_loss(actual: np.ndarray, probability: np.ndarray) -> float:
    probabilities = np.clip(probability, EPS, 1.0 - EPS)
    return float(-np.sum(actual * np.log(probabilities) + (1.0 - actual) * np.log(1.0 - probabilities)))


def fit_iid(sequences: Sequence[np.ndarray], prior: float = 0.5) -> np.ndarray:
    """Fit pooled symptom prevalences with a Jeffreys Beta prior."""
    counts = np.zeros((3, 2), dtype=float)
    for values in sequences:
        for symptom in range(3):
            observed = values[:, symptom]
            observed = observed[np.isfinite(observed)].astype(int)
            counts[symptom, 0] += np.sum(observed == 0)
            counts[symptom, 1] += np.sum(observed == 1)
    return (counts[:, 1] + prior) / (counts.sum(axis=1) + 2.0 * prior)


def fit_markov(sequences: Sequence[np.ndarray], prior: float = 0.5) -> tuple[np.ndarray, np.ndarray]:
    """Fit symptom-wise first-order transition probabilities."""
    iid = fit_iid(sequences, prior=prior)
    counts = np.zeros((3, 2, 2), dtype=float)
    for values in sequences:
        for symptom in range(3):
            previous = values[:-1, symptom]
            following = values[1:, symptom]
            mask = np.isfinite(previous) & np.isfinite(following)
            for old, new in zip(previous[mask].astype(int), following[mask].astype(int)):
                counts[symptom, old, new] += 1.0
    probabilities = (counts[:, :, 1] + prior) / (counts.sum(axis=2) + 2.0 * prior)
    return probabilities, iid


def _fit_burden_model(sequences: Sequence[np.ndarray], nonlinear: bool, seed: int) -> np.ndarray:
    """Fit the linear or active-inference-derived nonlinear transition surface."""
    burden, actual, symptom = _training_pairs(sequences)
    prevalence = fit_iid(sequences)
    initial_intercepts = np.log(np.clip(prevalence, EPS, 1.0 - EPS) / np.clip(1.0 - prevalence, EPS, 1.0 - EPS))

    if nonlinear:
        bounds = [(-10.0, 10.0)] * 3 + [(0.0, 10.0)] * 3 + [(0.0, 10.0), (-10.0, 10.0)]

        def objective(parameters: np.ndarray) -> float:
            intercept = parameters[:3]
            loading = parameters[3:6]
            coupling, field = parameters[6:8]
            mapped = np.tanh(coupling * burden + field)
            probability = expit(intercept[symptom] + loading[symptom] * mapped)
            penalty = 0.01 * float(np.sum(parameters[3:] ** 2))
            return _bernoulli_loss(actual, probability) + penalty

        starts = [np.r_[initial_intercepts, np.ones(3), 1.1, 0.0]]
    else:
        bounds = [(-10.0, 10.0)] * 3 + [(0.0, 10.0)] * 3

        def objective(parameters: np.ndarray) -> float:
            intercept = parameters[:3]
            loading = parameters[3:6]
            probability = expit(intercept[symptom] + loading[symptom] * burden)
            penalty = 0.01 * float(np.sum(parameters[3:] ** 2))
            return _bernoulli_loss(actual, probability) + penalty

        starts = [np.r_[initial_intercepts, np.ones(3)]]

    rng = np.random.default_rng(seed)
    for _ in range(7):
        starts.append(np.asarray([rng.uniform(low, high) for low, high in bounds]))

    best = None
    for start in starts:
        result = minimize(objective, start, method="L-BFGS-B", bounds=bounds, options={"maxiter": 2000})
        if best is None or result.fun < best.fun:
            best = result
    if best is None or not np.isfinite(best.fun):
        raise RuntimeError("Burden-model optimisation failed")
    return np.asarray(best.x)


def bistable_wedge_diagnostic(coupling: float, field: float) -> tuple[float, bool]:
    """Return the fold-field magnitude and whether (J, h) lies inside the wedge."""
    if coupling <= 1.0:
        return 0.0, False
    fold_state = math.sqrt(1.0 - 1.0 / coupling)
    critical_field = coupling * fold_state - np.arctanh(fold_state)
    return float(critical_field), bool(abs(field) < critical_field)


def _hmm_emission_likelihood(sequence: np.ndarray, emission: np.ndarray) -> np.ndarray:
    """Evaluate Bernoulli emission likelihoods while ignoring missing symptoms."""
    likelihood = np.ones((sequence.shape[0], 2), dtype=float)
    for time, row in enumerate(sequence):
        for symptom, value in enumerate(row):
            if np.isfinite(value):
                probability = np.clip(emission[:, symptom], EPS, 1.0 - EPS)
                likelihood[time] *= probability if value == 1 else 1.0 - probability
    return np.clip(likelihood, EPS, None)


def _forward_backward(sequence: np.ndarray, initial: np.ndarray, transition: np.ndarray, emission: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Scaled forward-backward recursion for a two-state Bernoulli HMM."""
    likelihood = _hmm_emission_likelihood(sequence, emission)
    length = sequence.shape[0]
    alpha = np.zeros((length, 2), dtype=float)
    scale = np.zeros(length, dtype=float)
    alpha[0] = initial * likelihood[0]
    scale[0] = max(alpha[0].sum(), EPS)
    alpha[0] /= scale[0]
    for time in range(1, length):
        alpha[time] = (alpha[time - 1] @ transition) * likelihood[time]
        scale[time] = max(alpha[time].sum(), EPS)
        alpha[time] /= scale[time]
    beta = np.ones((length, 2), dtype=float)
    for time in range(length - 2, -1, -1):
        beta[time] = transition @ (likelihood[time + 1] * beta[time + 1])
        beta[time] /= max(scale[time + 1], EPS)
    gamma = alpha * beta
    gamma /= np.clip(gamma.sum(axis=1, keepdims=True), EPS, None)
    return alpha, gamma, float(np.log(scale).sum())


def fit_hmm(sequences: Sequence[np.ndarray], seed: int, n_starts: int = 8) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fit a pooled two-state Bernoulli HMM by multi-start EM."""
    rng = np.random.default_rng(seed)
    best: tuple[float, np.ndarray, np.ndarray, np.ndarray] | None = None
    for start in range(n_starts):
        if start == 0:
            initial = np.array([0.5, 0.5])
            transition = np.array([[0.98, 0.02], [0.02, 0.98]])
            emission = np.array([[0.03, 0.05, 0.02], [0.85, 0.75, 0.35]])
        else:
            initial = rng.dirichlet(np.ones(2))
            transition = np.vstack([rng.dirichlet(np.array([8.0, 1.0])), rng.dirichlet(np.array([1.0, 8.0]))])
            emission = rng.uniform(0.02, 0.98, size=(2, 3))

        previous_loglik = -np.inf
        for _ in range(500):
            initial_count = np.full(2, 0.5)
            transition_count = np.full((2, 2), 0.5)
            emission_success = np.full((2, 3), 0.5)
            emission_total = np.ones((2, 3))
            total_loglik = 0.0
            for sequence in sequences:
                alpha, gamma, loglik = _forward_backward(sequence, initial, transition, emission)
                likelihood = _hmm_emission_likelihood(sequence, emission)
                initial_count += gamma[0]
                total_loglik += loglik
                for time in range(sequence.shape[0] - 1):
                    xi = alpha[time][:, None] * transition * (likelihood[time + 1] * (gamma[time + 1] / np.clip(alpha[time + 1], EPS, None)))[None, :]
                    xi /= max(xi.sum(), EPS)
                    transition_count += xi
                for symptom in range(3):
                    mask = np.isfinite(sequence[:, symptom])
                    if np.any(mask):
                        weights = gamma[mask]
                        emission_total[:, symptom] += weights.sum(axis=0)
                        emission_success[:, symptom] += (weights * sequence[mask, symptom, None]).sum(axis=0)
            initial = initial_count / initial_count.sum()
            transition = transition_count / transition_count.sum(axis=1, keepdims=True)
            emission = emission_success / emission_total
            if abs(total_loglik - previous_loglik) < 1e-7:
                break
            previous_loglik = total_loglik
        if best is None or total_loglik > best[0]:
            best = (total_loglik, initial.copy(), transition.copy(), emission.copy())
    if best is None:
        raise RuntimeError("HMM optimisation failed")
    return best[1], best[2], best[3]


def _append_predictions(rows: list[Prediction], model: str, subject: int, sequence: np.ndarray, probabilities: np.ndarray) -> None:
    """Append finite one-step predictions to the long-form score table."""
    for week in range(sequence.shape[0] - 1):
        for symptom, name in enumerate(SYMPTOM_NAMES):
            target = sequence[week + 1, symptom]
            if np.isfinite(target):
                previous = sequence[week, symptom]
                changed = bool(np.isfinite(previous) and previous != target)
                rows.append(Prediction(model, subject, week + 2, name, float(target), float(probabilities[week, symptom]), changed))


def leave_one_subject_out_predictions(sequences: Sequence[np.ndarray], seed: int = 20260825) -> tuple[list[Prediction], list[dict[str, float]]]:
    """Run leave-one-patient-out one-step forecasts for all comparison models."""
    predictions: list[Prediction] = []
    parameter_rows: list[dict[str, float]] = []
    for held_out in range(len(sequences)):
        training = [values for index, values in enumerate(sequences) if index != held_out]
        test = sequences[held_out]
        iid = fit_iid(training)
        markov, markov_iid = fit_markov(training)
        linear = _fit_burden_model(training, nonlinear=False, seed=seed + held_out)
        cusp = _fit_burden_model(training, nonlinear=True, seed=seed + 100 + held_out)
        hmm_initial, hmm_transition, hmm_emission = fit_hmm(training, seed=seed + 200 + held_out)

        number_pairs = test.shape[0] - 1
        iid_probability = np.tile(iid, (number_pairs, 1))
        markov_probability = np.zeros((number_pairs, 3), dtype=float)
        linear_probability = np.zeros((number_pairs, 3), dtype=float)
        cusp_probability = np.zeros((number_pairs, 3), dtype=float)
        hmm_probability = np.zeros((number_pairs, 3), dtype=float)

        filtered = hmm_initial.copy()
        for week in range(number_pairs):
            current = test[week]
            burden = _burden(current)
            linear_probability[week] = expit(linear[:3] + linear[3:6] * burden)
            mapped = math.tanh(cusp[6] * burden + cusp[7])
            cusp_probability[week] = expit(cusp[:3] + cusp[3:6] * mapped)
            for symptom in range(3):
                previous = current[symptom]
                markov_probability[week, symptom] = markov[symptom, int(previous)] if np.isfinite(previous) else markov_iid[symptom]
            likelihood = _hmm_emission_likelihood(current[None, :], hmm_emission)[0]
            filtered = filtered * likelihood
            filtered /= max(filtered.sum(), EPS)
            predictive_state = filtered @ hmm_transition
            hmm_probability[week] = predictive_state @ hmm_emission
            filtered = predictive_state

        subject = held_out + 1
        _append_predictions(predictions, "Intercept", subject, test, iid_probability)
        _append_predictions(predictions, "Markov", subject, test, markov_probability)
        _append_predictions(predictions, "Linear burden", subject, test, linear_probability)
        _append_predictions(predictions, "Two-state HMM", subject, test, hmm_probability)
        _append_predictions(predictions, "AI-derived nonlinear", subject, test, cusp_probability)
        critical_field, inside_wedge = bistable_wedge_diagnostic(float(cusp[6]), float(cusp[7]))
        parameter_rows.append({
            "held_out_subject": float(subject),
            "nonlinear_coupling_J": float(cusp[6]),
            "nonlinear_field_h": float(cusp[7]),
            "critical_field_abs": critical_field,
            "inside_bistable_wedge": int(inside_wedge),
            "hmm_persist_low": float(hmm_transition[0, 0]),
            "hmm_persist_high": float(hmm_transition[1, 1]),
        })
    return predictions, parameter_rows


def summarise_predictions(predictions: Sequence[Prediction]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Summarise forecast performance overall and by held-out patient."""
    models = list(dict.fromkeys(row.model for row in predictions))
    overall: list[dict[str, object]] = []
    subject_rows: list[dict[str, object]] = []
    for model in models:
        selected = [row for row in predictions if row.model == model]
        actual = np.asarray([row.actual for row in selected])
        probability = np.asarray([row.probability for row in selected])
        changed = np.asarray([row.changed for row in selected])
        log_loss = _bernoulli_loss(actual, probability) / len(selected)
        brier = float(np.mean((actual - probability) ** 2))
        transition_loss = _bernoulli_loss(actual[changed], probability[changed]) / int(changed.sum()) if np.any(changed) else math.nan
        overall.append({
            "model": model,
            "predictions": len(selected),
            "transitions": int(changed.sum()),
            "log_loss": log_loss,
            "brier_score": brier,
            "transition_log_loss": transition_loss,
        })
        for subject in sorted({row.subject for row in selected}):
            subset = [row for row in selected if row.subject == subject]
            y_value = np.asarray([row.actual for row in subset])
            p_value = np.asarray([row.probability for row in subset])
            subject_rows.append({
                "model": model,
                "subject": subject,
                "predictions": len(subset),
                "log_loss": _bernoulli_loss(y_value, p_value) / len(subset),
                "brier_score": float(np.mean((y_value - p_value) ** 2)),
            })
    overall.sort(key=lambda row: float(row["log_loss"]))
    return overall, subject_rows


def audit_lorenz_estimates(path: Path) -> list[dict[str, object]]:
    """Audit the supplied in-sample Lorenz fits without treating them as forecasts."""
    raw = loadmat(path, simplify_cells=True)
    fits = list(np.atleast_1d(raw["GCM"]))
    rows: list[dict[str, object]] = []
    for subject, fit in enumerate(fits, start=1):
        actual = np.asarray(fit["Y"], dtype=float)
        predicted = np.asarray(fit["y"], dtype=float)
        mask = np.isfinite(actual) & np.isfinite(predicted)
        parameters = np.asarray(fit["Ep"]["A"], dtype=float).reshape(-1)
        inputs = np.asarray(fit["U"], dtype=float).reshape(-1)
        effective_rho = (1.0 - inputs * parameters[2]) * parameters[1]
        latent = np.asarray(fit["x"], dtype=float)
        signs = np.sign(latent[:, 0])
        crossings = int(np.sum(signs[1:] * signs[:-1] < 0))
        rows.append({
            "subject": subject,
            "free_energy": float(fit["F"]),
            "in_sample_rmse": float(np.sqrt(np.mean((actual[mask] - predicted[mask]) ** 2))),
            "baseline_rho": float(parameters[1]),
            "minimum_effective_rho": float(np.min(effective_rho)),
            "maximum_effective_rho": float(np.max(effective_rho)),
            "weeks_above_chaos_threshold": int(np.sum(effective_rho > LORENZ_CHAOS_THRESHOLD)),
            "latent_lobe_crossings": crossings,
        })
    return rows


def _write_csv(path: Path, rows: Sequence[dict[str, object]]) -> None:
    """Write a non-empty list of dictionaries to CSV."""
    if not rows:
        raise ValueError(f"No rows available for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def make_figure(sequences: Sequence[np.ndarray], descriptives: Sequence[dict[str, object]], overall: Sequence[dict[str, object]], subject_metrics: Sequence[dict[str, object]], lorenz_rows: Sequence[dict[str, object]], path: Path) -> None:
    """Create the four-panel manuscript figure for the clinical stress test."""
    colors = {"Psychosis": "#355C7D", "Depression": "#2A9D8F", "Mania": "#E76F51"}
    model_order = ["Intercept", "Markov", "Linear burden", "Two-state HMM", "AI-derived nonlinear"]
    model_colors = ["#B7B7B7", "#355C7D", "#7A5195", "#2A9D8F", "#E76F51"]
    figure, axes = plt.subplots(2, 2, figsize=(12.0, 8.7), constrained_layout=True)

    raster = np.full((len(sequences) * 3, max(values.shape[0] for values in sequences)), np.nan)
    for subject, values in enumerate(sequences):
        raster[subject * 3 : subject * 3 + 3, : values.shape[0]] = values.T
    axes[0, 0].imshow(raster, aspect="auto", interpolation="nearest", cmap="Greys", vmin=0, vmax=1)
    axes[0, 0].set_yticks([subject * 3 + 1 for subject in range(len(sequences))])
    axes[0, 0].set_yticklabels([f"P{subject + 1}" for subject in range(len(sequences))])
    axes[0, 0].set_xlabel("Week")
    axes[0, 0].set_title("(a)  Released symptom histories")

    transition_matrix = np.asarray([[row[f"{name.lower()}_transitions"] for name in SYMPTOM_NAMES] for row in descriptives])
    bottom = np.zeros(len(sequences))
    for symptom, name in enumerate(SYMPTOM_NAMES):
        axes[0, 1].bar(np.arange(1, len(sequences) + 1), transition_matrix[:, symptom], bottom=bottom, color=colors[name], label=name)
        bottom += transition_matrix[:, symptom]
    total_pairs = int(sum(int(row["adjacent_pairs"]) for row in descriptives))
    total_transitions = int(transition_matrix.sum())
    persistence = 1.0 - total_transitions / total_pairs
    axes[0, 1].text(0.98, 0.95, f"{total_transitions}/{total_pairs} changes\nPersistence = {persistence:.2%}", transform=axes[0, 1].transAxes, ha="right", va="top")
    axes[0, 1].set_xlabel("Patient")
    axes[0, 1].set_ylabel("Observed transitions")
    axes[0, 1].set_xticks(np.arange(1, len(sequences) + 1))
    axes[0, 1].legend(frameon=False, fontsize=8)
    axes[0, 1].set_title("(b)  Transition information")

    overall_map = {str(row["model"]): float(row["log_loss"]) for row in overall}
    positions = np.arange(len(model_order))
    axes[1, 0].bar(positions, [overall_map[name] for name in model_order], color=model_colors, alpha=0.85)
    for position, model in enumerate(model_order):
        values = [float(row["log_loss"]) for row in subject_metrics if row["model"] == model]
        jitter = np.linspace(-0.10, 0.10, len(values))
        axes[1, 0].scatter(np.full(len(values), position) + jitter, values, s=18, facecolor="white", edgecolor="black", linewidth=0.6, zorder=3)
    axes[1, 0].set_xticks(positions)
    axes[1, 0].set_xticklabels(["Intercept", "Markov", "Linear", "HMM", "AI nonlinear"], rotation=20, ha="right")
    axes[1, 0].set_ylabel("Leave-one-patient-out log loss")
    axes[1, 0].set_title("(c)  Held-out one-step prediction")

    subjects = np.arange(1, len(lorenz_rows) + 1)
    minimum = np.asarray([float(row["minimum_effective_rho"]) for row in lorenz_rows])
    maximum = np.asarray([float(row["maximum_effective_rho"]) for row in lorenz_rows])
    baseline = np.asarray([float(row["baseline_rho"]) for row in lorenz_rows])
    axes[1, 1].vlines(subjects, minimum, maximum, color="#355C7D", linewidth=3)
    axes[1, 1].scatter(subjects, baseline, color="white", edgecolor="#355C7D", zorder=3)
    axes[1, 1].axhline(LORENZ_CHAOS_THRESHOLD, color="#E76F51", linestyle="--", linewidth=1.5, label="Classical stability threshold")
    axes[1, 1].set_xlabel("Patient")
    axes[1, 1].set_ylabel("Effective Rayleigh parameter")
    axes[1, 1].set_xticks(subjects)
    axes[1, 1].legend(frameon=False, fontsize=8)
    axes[1, 1].set_title("(d)  Supplied Lorenz-fit audit")

    for axis in axes.ravel():
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def run_analysis(data_dir: Path = DEFAULT_DATA_DIR, output_dir: Path = DEFAULT_OUTPUT_DIR, seed: int = 20260825) -> dict[str, object]:
    """Run the full stress test and write manuscript-ready aggregate outputs."""
    trajectories, metadata = load_trajectories(data_dir / "Data_for_Model.mat")
    descriptives = describe_sequences(trajectories)
    predictions, parameters = leave_one_subject_out_predictions(trajectories, seed=seed)
    overall, subject_metrics = summarise_predictions(predictions)
    lorenz_rows = audit_lorenz_estimates(data_dir / "GCM_nosology.mat")

    output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(output_dir / "trajectory_descriptives.csv", descriptives)
    _write_csv(output_dir / "model_comparison.csv", overall)
    _write_csv(output_dir / "subject_model_metrics.csv", subject_metrics)
    _write_csv(output_dir / "cross_validation_parameters.csv", parameters)
    _write_csv(output_dir / "lorenz_fit_audit.csv", lorenz_rows)
    make_figure(trajectories, descriptives, overall, subject_metrics, lorenz_rows, output_dir / "clinical_trajectory_stress_test.png")

    summary = {
        "metadata": metadata,
        "patients": len(trajectories),
        "weeks_per_patient": [int(values.shape[0]) for values in trajectories],
        "adjacent_pairs": int(sum(int(row["adjacent_pairs"]) for row in descriptives)),
        "transitions": int(sum(int(row["transitions"]) for row in descriptives)),
        "models": overall,
    }
    return summary


def run_all() -> dict[str, object]:
    """Repository-orchestrator entry point."""
    return run_analysis()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--seed", type=int, default=20260825)
    arguments = parser.parse_args()
    summary = run_analysis(arguments.data_dir, arguments.output_dir, arguments.seed)
    print(f"Patients: {summary['patients']}")
    print(f"Finite adjacent symptom pairs: {summary['adjacent_pairs']}")
    print(f"Observed changes: {summary['transitions']}")
    print("Held-out model ranking (lower log loss is better):")
    for row in summary["models"]:
        print(f"  {row['model']}: log loss={row['log_loss']:.6f}, Brier={row['brier_score']:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
