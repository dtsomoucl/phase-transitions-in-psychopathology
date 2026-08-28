# Supplementary Online Material (SOM)

## Evidence map

| Level | Evidence | Supported conclusion | Principal limit | Appendix route |
|---|---|---|---|---|
| Formal | Mean-field derivation and finite-agent simulations | Action-dependent learning can produce an intermediate multistable window, fold loss, history-dependent recovery, and conditional entrenchment | The simulations are theoretical probes, not fitted clinical models | Notebooks 01–02 and `Python_code/` |
| Primary empirical | Wichers N-of-1 intensive study | A blind boundary search selects the same weekly occasion as the published analysis; a discontinuity is supported against a constant mean | Same series as the original report; the step is only modestly favoured over delayed discontinuation after search penalty | Appendix A |
| Secondary intensive | Fisher and Jang | Accumulated action history predicts subsequent behaviour | Simpler history models perform best; mechanism proxies and temporal order are incomplete | Appendices B–C |
| Secondary boundary | Jin sparse symptom histories | Persistence is a stringent benchmark and model performance is poor at rare changes | Symptoms do not measure actions, outcomes, preferences, or policy precision | Appendix D |

The paper therefore claims **empirical compatibility** with, and _partial_ support for, path-dependent dynamics. It does not claim unique identification of the specific action-contingent active-inference learner.

## Appendix-to-repository crosswalk

- [Appendix A — Wichers N-of-1 transition](VALIDATION/appendices/A_wichers/NOTES_A.md)
- [Appendix B — Fisher avoidance/contentment](VALIDATION/appendices/B_fisher/NOTES_B.md)
- [Appendix C — Jang exercise/panic](VALIDATION/appendices/C_jang/NOTES_C.md)
- [Appendix D — Jin sparse symptom trajectories](VALIDATION/appendices/D_jin/NOTES_D.md)

## Canonical manuscript figures

| Figure | Analysis | Code |
|---|---|---|
| 1 | Learning–action loop and evidence hierarchy | `Python_code/export_figure1_concept.py` |
| 2 | Onset schedule, fixed points, kinetic lag, and ensemble variance | `Python_code/sim_psychopathology.py` |
| 3 | Paired recovery from verified withdrawal | `Python_code/recovery_bistable.py` |
| 4 | Fold region, finite-agent boundary, and live regime audit | `Python_code/recovery_bistable.py` |
| 5 | History dependence and retention-sensitive recovery thresholds | `Python_code/recovery_bistable.py` |
| 6 | Wichers density, boundary recovery, and penalised model comparison | `VALIDATION/generate_wichers_figure.R` |

Raster and vector copies are stored in `outputs/manuscript_figures/`.

## Reproduction routes

```bash
python3 -m pip install -r requirements.txt
python3 main.py
```

For analysis-specific reruns, see `README.md` and `VALIDATION/README.md`.

## Brief notes on the key claims

1. **Multistability:** established analytically only where the fixed-point audit returns two stable branches and one unstable branch.
2. **Finite-agent onset:** the live schedule loses the engaged quasistatic branch, but stochastic trajectories show substantial kinetic lag; fold timing is not treated as the observed transition time.
3. **Recovery and path dependence:** demonstrated from paired learned states inside the bistable regime. Differential retention raises recovery thresholds more strongly than symmetric or permanent evidence retention.
4. **Wichers transition:** the same weekly occasion is recovered with substantial bootstrap concentration, but the endpoint does not identify the learning mechanism.
5. **Secondary prediction:** cumulative history improves prediction mainly during persistent periods; transition-only performance is reported separately.
6. **Causality:** not identified because actions and the learning loop were not experimentally manipulated.
