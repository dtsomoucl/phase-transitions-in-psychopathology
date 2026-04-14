"""
Shared parameter regimes for the psychopathology simulations.

Three regimes are defined:

ONSET_REGIME
    High prior precision (α₀ = 40), slow erosion of γ and Δc during
    adversity.  Used in Step 1 (sim_psychopathology.py) and in the
    cross-regime robustness check.  Chosen so the live adversity schedule
    enters the quasistatic multistable window while policy precision remains
    above the critical level — this supports a genuine catastrophe-style
    onset claim under the actual simulated schedule.

RECOVERY_REGIME
    Low prior precision (α₀ = 2), fast erosion of γ during illness with
    slow erosion of Δc.  Used in Steps 2–5 (hysteresis, robustness,
    orthogonal interventions, decay, asymmetric memory).  The qualitative
    difference from ONSET_REGIME (α₀ factor of 20, γ_rate factor of 10)
    reflects a less-informed prior at illness onset; the field-dominance
    and hysteresis conclusions are robust across both regimes (see
    fig_S2_cross_regime.png).

ONSET_AS_RECOVERY_REGIME
    ONSET_REGIME parameters re-run through the recovery protocol to confirm
    that field dominance of recovery is not an artefact of the lower-prior
    RECOVERY_REGIME.  Produced by fig_cross_regime_field_dominance() in
    step2_robustness.py.
"""

### DT ---> The onset regime is chosen so the live adversity schedule enters
### DT ---> the quasistatic multistable window while policy precision remains
### DT ---> above the critical level. This supports a genuine catastrophe-style
### DT ---> onset claim under the actual simulated schedule.
ONSET_REGIME = dict(
    p=0.85,
    alpha_0=40.0,
    gamma_healthy=16.0,
    delta_c_healthy=0.3,
    N_healthy=200,
    gamma_rate=0.005,
    delta_c_rate=0.003,
    gamma_floor=14.0,
    delta_c_floor=-0.15,
)

### DT ---> The recovery regime keeps the original lower-prior, stronger-erosion
### DT ---> schedule used in the field-dominance analyses so the Step 2-4
### DT ---> recovery results remain directly comparable.
RECOVERY_REGIME = dict(
    p=0.85,
    alpha_0=2.0,
    gamma_healthy=16.0,
    delta_c_healthy=0.15,
    N_healthy=200,
    gamma_rate=0.05,
    delta_c_rate=0.00075,
    gamma_floor=5.0,
    delta_c_floor=-0.15,
)

### DT ---> This regime runs the onset-regime parameters through the recovery
### DT ---> protocol. Purpose: confirm that field dominance does not depend on
### DT ---> the low-prior assumption of RECOVERY_REGIME. See W1 in review notes.
ONSET_AS_RECOVERY_REGIME = dict(
    p=0.85,
    alpha_0=40.0,           # same as ONSET_REGIME
    gamma_healthy=16.0,
    delta_c_healthy=0.15,   # recovery-protocol starting Δc
    N_healthy=200,
    gamma_rate=0.05,        # same erosion rate as RECOVERY_REGIME
    delta_c_rate=0.00075,
    gamma_floor=5.0,
    delta_c_floor=-0.15,
)


def onset_regime():
    return ONSET_REGIME.copy()


def recovery_regime():
    return RECOVERY_REGIME.copy()


def onset_as_recovery_regime():
    return ONSET_AS_RECOVERY_REGIME.copy()
