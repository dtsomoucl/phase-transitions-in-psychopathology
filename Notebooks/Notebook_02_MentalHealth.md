# Notebook 02 — Clinical Translation and Bistable Recovery

**Series:** *Phase Transitions in Psychopathology*  
**Status:** Interpretation of the prespecified onset and recovery analyses  
**Author:** Dimitris I. Tsomokos

## 1. Clinical translation

The two policies are interpreted as broad courses of action: engagement and withdrawal. This is a worked clinical translation of a minimal model, not a fitted model of depression and not a claim that psychopathology has only two states.

| Model quantity | Working interpretation | Evidential status |
|---|---|---|
| $\pi_1$ | engagement | model policy |
| $\pi_2$ | withdrawal | model policy |
| $z$ | accumulated action-allocation balance | latent model quantity |
| $P(\pi_1)$ | probability of selecting engagement | latent model quantity |
| $\Delta c$ | relative preference field | model control |
| $\gamma$ | policy precision | model control |

Increasing $\gamma$ is not a generic improvement: it sharpens selection of whichever policy currently has lower expected free energy. Changing $\Delta c$ can reverse which policy is favoured.

## 2. Onset schedule and kinetic lag

The onset schedule uses $\alpha_0=40$, $p=.85$, an initial $\gamma=16$, and an initial $\Delta c=.30$. After 200 healthy trials, $\gamma$ and $\Delta c$ drift gradually toward 14 and −.15.

The live quasistatic schedule loses the engaged stable branch at trial 310. The finite stochastic process does not relax instantaneously to that moving fixed point. Of 100 prespecified agents, 60 meet the sustained-withdrawal criterion by trial 600; their median crossing is trial 469.5 (IQR 382.25–520), a median lag of 159.5 trials after fold loss. Among the 21 agents with an observable smoothed .80-to-.20 transition, the median duration is 162 trials (IQR 130–232), or 27% of the full simulation window.

The onset analysis therefore demonstrates a fold boundary and finite-agent transition under gradual drift, together with substantial kinetic lag. It does not support describing every simulated trajectory as an instantaneous jump. The late ensemble-variance peak is a retrospective heterogeneity diagnostic, not a prospective clinical early-warning signal.

## 3. Structurally valid recovery protocol

Recovery is evaluated under a separate reference protocol selected by structural criteria. It uses $\alpha_0=40$, $p=.85$, healthy controls $(\gamma,\Delta c)=(20,.15)$, and adverse floors $(18,-.05)$. The agent experiences 100 healthy trials, 100 trials of linear drift, and 100 trials at the adverse floors before intervention at trial 300. Thus $\tau=300/(2\alpha_0)=3.75$.

A withdrawn state is verified only when both floors have been reached and both instantaneous $P(engaged)$ and its trailing 30-trial mean remain below .50 for 30 consecutive trials. Two hundred verified agents are selected from deterministic candidate seeds. All intervention branches start from the same learned state and random-number-generator state.

At intervention, the control condition has two stable mean-field branches, $z^*=-.961$ and $z^*=.902$, separated by an unstable branch at $z^*=.613$. Restoring $\gamma$ alone retains three fixed points. Restoring $\Delta c$ collapses the system to one engaged fixed point ($z^*=.979$), and restoring both gives $z^*=.987$. The recovery analysis therefore takes place in a regime where bistability and hysteresis are possible.

## 4. Paired recovery result

Recovery is summarised as mean latent $P(engaged)$ over the final 100 of 250 post-intervention trials. The .50 line is a model criterion, not a clinical definition.

| Condition | Median | IQR | Proportion above .50 |
|---|---:|---:|---:|
| Control | .156 | [.077, .278] | .025 |
| Restore $\gamma$ only | .132 | [.060, .245] | .025 |
| Restore $\Delta c$ only | .829 | [.657, .893] | .860 |
| Restore both | .862 | [.693, .918] | .860 |

The field restoration removes the withdrawn branch and produces re-engagement for most agents. Precision restoration alone slightly strengthens selection of the already favoured withdrawn policy. This is a model-internal result, not an estimate of treatment efficacy.

## 5. History dependence and retention

Engaged and withdrawn histories were brought to identical post-intervention controls. Their late engagement probabilities remained separated across a range of restored fields, demonstrating path dependence inside the bistable region.

Recovery thresholds were then estimated on a $\Delta c$ grid with resolution .01 and 1,000 agent-bootstrap resamples per cell. With 0, 50, 100, and 150 additional adverse trials after verified withdrawal, the estimated thresholds were:

| Retention rule | 0 | 50 | 100 | 150 |
|---|---:|---:|---:|---:|
| Permanent evidence | .06 | .07 | .09 | .09 |
| Symmetric decay | .05 | .06 | .07 | .08 |
| Mild differential retention | .14 | .18 | .22 | .25 |
| Strong differential retention | .21 | .25 | .29 | >.30 |

Threshold uncertainty is reported in `outputs/tables/bistable_retention_thresholds.csv`. The result supports residence-time entrenchment in this implementation, particularly when engagement-consistent evidence decays faster than withdrawal-consistent evidence. Differential retention is imposed rather than estimated from clinical data.

## 6. Reproducibility map

- `Python_code/core_functions.py`: expected-free-energy and fixed-point utilities.
- `Python_code/sim_psychopathology.py`: onset schedule, operational transition metrics, and Figure 2.
- `Python_code/recovery_bistable.py`: verified-state selection, branch enumeration, recovery, history, retention, and Figures 3–5.
- `outputs/tables/onset_*`: kinetic-lag and abruptness diagnostics.
- `outputs/tables/bistable_*`: fixed points, regime audit, paired recovery, structural sensitivity, and threshold uncertainty.
