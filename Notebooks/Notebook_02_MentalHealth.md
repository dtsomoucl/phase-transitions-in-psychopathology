# Notebook 02 — Clinical Translation and Recovery Boundary Conditions

**Series:** *Phase Transitions in Psychopathology*  
**Status:** Interpretation of the simulations, including the state-anchored sensitivity  
**Author:** Dimitris I. Tsomokos

## 1. Clinical translation

The two policies are interpreted as broad courses of action: engagement and withdrawal. This is a worked _clinical translation_ of a minimal model, not a fitted model of depression and not a claim that psychopathology has only two states.

| Model quantity | Working interpretation | Evidential status |
|---|---|---|
| $\pi_1$ | engagement | model policy |
| $\pi_2$ | withdrawal | model policy |
| $z$ | accumulated allocation balance | latent model quantity |
| $P(\pi_1)$ | probability of selecting engagement | latent model quantity |
| $\Delta c$ | relative preference field | model control |
| $\gamma$ | policy precision | model control |
| cohort engagement composites | engagement-related proxies | observational measures, not $\Delta c$ |
| cognitive-task scores | precision-related proxies | observational measures, not $\gamma$ |

The clinical labels “motivational field” and “cognitive control” are interpretive shorthand. Increasing $\gamma$ is not a generic improvement: it sharpens selection of whichever policy currently has lower expected free energy.

## 2. Simulation schedules

Two parameter regimes answer different theoretical questions.

- `ONSET_REGIME` uses $\alpha_0=40$, slow adverse drift, and a high precision floor. It tests whether gradual drift can traverse a multistable mean-field region and generate abrupt finite-agent changes.
- `RECOVERY_REGIME` uses $\alpha_0=2$, faster precision erosion, and a slower field drift. It separates the consequences of restoring $\Delta c$ and $\gamma$.

These are theoretical probes. They are not fitted to MCS or ABCD and should not be combined into one calibrated developmental trajectory.

## 3. Onset result

Under the onset schedule, many agents show a relatively abrupt fall in latent $P(engaged)$ while the controls change gradually. Agents cross at different times, so the ensemble mean is smoother. The quasistatic mean-field analysis enters a multistable region before the engaged branch is lost.

The ensemble variance peak is a retrospective simulation diagnostic, not a prospective early-warning signal. It occurs after the quasistatic engaged-branch fold loss and partly reflects heterogeneous transition timing across agents. Event alignment additionally requires the transition time to be known in advance (these results should therefore not be presented as evidence for a clinical early-warning marker).

## 4. Fixed-schedule recovery result

The original recovery comparison intervenes at trial 400, which is 200 trials after the adverse schedule begins. At that trial, $\gamma=6$ and $\Delta c=0$; both adverse floors are not reached until trial 600.

Therefore, the legacy code variable `T_ill` means **elapsed adverse-schedule exposure after the healthy phase**. It does not mean verified clinical illness duration or verified residence time in withdrawal.

Within the fixed protocol, restoring $\Delta c$ produces more re-engagement than restoring $\gamma$ alone. This follows from their different mathematical roles: changing $\Delta c$ can reverse which policy is favoured, whereas increasing $\gamma$ sharpens the current ordering.

## 5. State-anchored recovery check

`step6_verified_withdrawal.py` diagnoses the original intervention state and then repeats the comparison under a stricter protocol.

### 5.1 Diagnostic of the original intervention time

Across the 100 seeds used for the main recovery comparison:

- median current $P(engaged)=.523$, IQR $[.446,.638]$;
- 41% had current $P(engaged)<.50$;
- 34% had a prior-15-trial mean below .50; and
- 18% had cumulative order parameter $z<0$.

The original intervention ensemble is therefore not a uniformly consolidated-withdrawal sample.

### 5.2 Verified-withdrawal protocol

The sensitivity waits until both adverse controls are at their floors, then requires $P(engaged)<.50$ for 30 consecutive trials. The intervention comparison occurs after at least 200 further adverse-schedule trials, with current and recent $P(engaged)$ again below .50. For each agent, all four intervention branches start from identical Dirichlet evidence and an identical random-number-generator state.

With 200 paired agents, median late $P(engaged)$ is:

| Condition | Median | IQR |
|---|---:|---:|
| control | .371 | [.347, .400] |
| restore $\gamma$ only | .157 | [.114, .212] |
| restore $\Delta c$ only | .630 | [.603, .654] |
| restore both | .849 | [.803, .893] |

The paired $\Delta c$-over-$\gamma$ advantage has median .464, IQR [.434, .488], and is positive for all 200 agents. Thus, the model-internal field-versus-precision ordering survives the state-anchored check.

## 6. Memory and duration: revised boundary condition

The fixed-schedule memory figures vary the interval from adversity onset to intervention. Because the controls are still drifting over part of this interval, the protocol combines:

1. elapsed exposure;
2. changing $\gamma$ and $\Delta c$;
3. evolving evidence counts; and
4. any imposed retention asymmetry.

The state-anchored sensitivity instead varies the number of additional adverse trials after a verified withdrawal criterion. Under strong retention asymmetry, the minimum restored-field thresholds across 0, 50, 100, 200, and 400 additional trials are .060, .060, .045, .040, and .040. Mild and symmetric conditions are also non-monotonic.

Accordingly:

- asymmetric retention is sufficient to generate an exposure-duration pattern in the original fixed schedule;
- it does not, in this implementation, generate a monotonic recovery cost with time elapsed after verified withdrawal; and
- the model does not establish that longer clinical illness or longer withdrawal residence progressively raises the recovery threshold.

This negative boundary-condition result should be retained. A verified residence-time effect would require a different or additional history-dependent mechanism and a protocol that separates elapsed state residence from control drift.

## 7. Empirical bridge

MCS and ABCD contain broad engagement and cognitive-task measures at widely spaced waves. They do not contain the intensive behavioural observations needed to estimate $\Delta c$, $\gamma$, $z$, or individual transition boundaries (that is, we leave such phenotyping for future work with more intensive longitudinal data).

The cohort analyses therefore ask only whether **engagement-related measures have stronger adjusted prospective associations with later internalising outcomes than cognitive-task measures**. This is a heuristic population-level _consistency check_. Differences can reflect reliability, method, informant, scaling, or construct proximity and cannot identify the active-inference mechanism.

## 8. Reproducibility map

- `Python_code/step2_hysteresis.py`: fixed-schedule recovery comparison and boundary.
- `Python_code/step4_decay.py`: symmetric-retention result.
- `Python_code/step5_asymmetric_memory.py`: asymmetric retention under the fixed schedule.
- `Python_code/step6_verified_withdrawal.py`: state diagnostic and event-anchored sensitivity.
- `outputs/tables/verified_withdrawal_*.csv`: numerical results.
- `Python_code/Figs_psychopathology/fig_S6_verified_withdrawal.*`: supplementary figure.

