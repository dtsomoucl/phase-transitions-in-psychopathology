"""Prespecified parameter schedule for the onset analysis."""


ONSET_REGIME = {
    "p": 0.85,
    "alpha_0": 40.0,
    "gamma_healthy": 16.0,
    "delta_c_healthy": 0.30,
    "N_healthy": 200,
    "gamma_rate": 0.005,
    "delta_c_rate": 0.003,
    "gamma_floor": 14.0,
    "delta_c_floor": -0.15,
}


def onset_regime() -> dict:
    """Return an independent copy of the onset parameters."""
    return ONSET_REGIME.copy()
