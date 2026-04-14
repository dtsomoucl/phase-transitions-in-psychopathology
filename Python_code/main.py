"""
main.py — Central orchestrator for the psychopathology phase-transition project
================================================================================
Runs all simulation scripts in sequence and produces exactly the figures that
back up the claims in the analysis notebooks.

Figures produced (10 total)
---------------------------
Step 1 — Onset dynamics (sim_psychopathology.py)
  fig_P2_onset_revised.png        Ensemble, schedule panels, quasistatic FPs,
                                   population variance.  Panel (b) is a three-
                                   panel split (P(eng) | γ(t) | Δc(t)).
  fig_P_flags_revised.png         Catastrophe flags: bimodality at variance-peak
                                   trial, trajectory divergence, sudden jumps.
  fig_P3_ews_revised.png          EWS: event-aligned lag-1 AC (critical slowing
                                   down) and population susceptibility.

Step 2 — Recovery boundary and field dominance (step2_hysteresis.py)
  fig_S2_recovery_boundary.png    Heuristic recovery boundary in (γ, Δc) space
                                   and simulation verification.
  fig_S2_clinical_prediction.png  Four-condition recovery comparison and
                                   individual trajectories.

Step 2 robustness (step2_robustness.py)
  fig_S2_cross_regime.png         Cross-regime check: field dominance holds under
                                   both RECOVERY_REGIME (α₀=2) and
                                   ONSET_AS_RECOVERY_REGIME (α₀=40).  [NEW]

Step 3 — Orthogonal interventions (step3_orthogonal.py)
  fig_S3_combined.png             Bar chart of combined-intervention recovery
                                   when Δc is only partially restored.

Step 4 — Symmetric evidence decay (step4_decay.py)
  fig_D1_decay_hysteresis.png     Null result: symmetric decay does not produce
                                   T_ill-dependent trapping.

Step 5 — Asymmetric memory extension (step5_asymmetric_memory.py)
  fig_S5_asymmetric_memory.png    Strategy-specific retention creates genuine path
                                   asymmetry (with IQR bands; n=150 per cell).

Empirical bridge (empirical_bridge.py)
  fig_S_bridge_empirical.png      Synthetic internal consistency check:
                                   field-dominant DGP calibrated to MCS
                                   coefficients; balance vs component classifiers.

Figures retained as callable functions but NOT produced here
------------------------------------------------------------
  fig_P_reversal_true             True reversal protocol (Step 1 supplementary)
  fig_S2_threshold_sensitivity    Step 2 threshold robustness
  fig_S2_regime_dependence        Step 2 regime heatmap
  fig_S3_tempering                Step 3 count-tempering heatmap
  fig_S3_p_increase               Step 3 environmental restructuring
  fig_S3_threshold_sensitivity    Step 3 threshold sensitivity
  fig_S3_regime_dependence        Step 3 regime heatmap
  fig_D2_decay_contrast           Step 4 decay contrast (supplementary)

Usage
-----
    cd Python_code
    python main.py

All figures are written to Python_code/Figs_psychopathology/.
Tables are written to outputs/tables/.
"""

import os
import sys
import time

### DT ---> Ensure Python_code/ is on the path so imports resolve correctly
### DT ---> regardless of where the script is called from.
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

FIG_DIR = os.path.join(_HERE, "Figs_psychopathology")
TABLE_DIR = os.path.join(_HERE, "..", "outputs", "tables")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)


def separator(title):
    width = 70
    print("\n" + "=" * width)
    print(f"  {title}")
    print("=" * width)


def run_step(name, module_run_all):
    t0 = time.time()
    try:
        module_run_all()
        elapsed = time.time() - t0
        print(f"  [{name}] done in {elapsed:.1f}s")
        return True
    except Exception as exc:
        print(f"\n  [ERROR in {name}]: {exc}")
        import traceback
        traceback.print_exc()
        return False


def main():
    t_total = time.time()
    errors = []

    print("=" * 70)
    print("  Psychopathology Phase-Transition Project — full simulation sweep")
    print("=" * 70)
    print(f"  Output directory: {FIG_DIR}")
    print(f"  Table directory:  {TABLE_DIR}")

    # ── Step 1: Onset dynamics ─────────────────────────────────────────
    separator("Step 1: Onset dynamics (sim_psychopathology.py)")
    import sim_psychopathology
    if not run_step("Step 1", sim_psychopathology.run_all):
        errors.append("Step 1: sim_psychopathology")

    # ── Step 2a: Recovery boundary and field dominance ─────────────────
    separator("Step 2a: Recovery boundary (step2_hysteresis.py)")
    import step2_hysteresis
    if not run_step("Step 2a", step2_hysteresis.run_all):
        errors.append("Step 2a: step2_hysteresis")

    # ── Step 2b: Cross-regime robustness ───────────────────────────────
    separator("Step 2b: Robustness — cross-regime field dominance (step2_robustness.py)")
    import step2_robustness
    if not run_step("Step 2b", step2_robustness.run_all):
        errors.append("Step 2b: step2_robustness")

    # ── Step 3: Orthogonal interventions ──────────────────────────────
    separator("Step 3: Orthogonal interventions (step3_orthogonal.py)")
    import step3_orthogonal
    if not run_step("Step 3", step3_orthogonal.run_all):
        errors.append("Step 3: step3_orthogonal")

    # ── Step 4: Symmetric evidence decay ──────────────────────────────
    separator("Step 4: Symmetric evidence decay — null result (step4_decay.py)")
    import step4_decay
    if not run_step("Step 4", step4_decay.run_all):
        errors.append("Step 4: step4_decay")

    # ── Step 5: Asymmetric memory extension ───────────────────────────
    separator("Step 5: Asymmetric memory extension (step5_asymmetric_memory.py)")
    import step5_asymmetric_memory
    if not run_step("Step 5", step5_asymmetric_memory.run_all):
        errors.append("Step 5: step5_asymmetric_memory")

    # ── Empirical bridge ───────────────────────────────────────────────
    separator("Empirical bridge: synthetic consistency check (empirical_bridge.py)")
    import empirical_bridge
    if not run_step("Empirical bridge", empirical_bridge.run_all):
        errors.append("Empirical bridge")

    # ── Summary ────────────────────────────────────────────────────────
    separator("Summary")
    elapsed_total = time.time() - t_total

    ### DT ---> Check which of the 10 expected figures were actually created
    expected_figures = [
        "fig_P2_onset_revised.png",
        "fig_P_flags_revised.png",
        "fig_P3_ews_revised.png",
        "fig_S2_recovery_boundary.png",
        "fig_S2_clinical_prediction.png",
        "fig_S2_cross_regime.png",
        "fig_S3_combined.png",
        "fig_D1_decay_hysteresis.png",
        "fig_S5_asymmetric_memory.png",
        "fig_S_bridge_empirical.png",
    ]

    print(f"\n  Total elapsed: {elapsed_total / 60:.1f} min\n")
    print("  Figure check:")
    n_found = 0
    for fname in expected_figures:
        fpath = os.path.join(FIG_DIR, fname)
        exists = os.path.isfile(fpath)
        status = "OK" if exists else "MISSING"
        print(f"    [{status}]  {fname}")
        if exists:
            n_found += 1

    print(f"\n  {n_found}/{len(expected_figures)} expected figures found.")

    if errors:
        print("\n  Steps with errors:")
        for e in errors:
            print(f"    - {e}")
        print("\n  Re-run the affected step script individually for the full traceback.")
    else:
        print("\n  All steps completed without errors.")

    print("\n" + "=" * 70)
    return 0 if not errors else 1


if __name__ == "__main__":
    sys.exit(main())
