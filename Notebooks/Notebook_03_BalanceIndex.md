# Notebook 03 — Policy Log Odds and Exploratory Questionnaire Ratios

**Series:** *Phase Transitions in Psychopathology*  
**Status:** Measurement note  
**Author:** Dimitris I. Tsomokos

## 1. Two distinct objects

This notebook separates an exact identity inside the generative model from an exploratory ratio constructed from cohort questionnaires.

They are not interchangeable:

1. the policy-posterior log odds are defined on model probabilities; and
2. distress-to-wellbeing ratios combine observed scales with unknown measurement functions.

## 2. Exact model-internal identity

For the two-policy softmax,

```math
P(\pi_1)=\frac{1}{1+\exp(\gamma\Delta G)},
\qquad
P(\pi_2)=1-P(\pi_1).
```

Therefore,

```math
\boxed{\log\frac{P(\pi_2)}{P(\pi_1)}=\gamma\Delta G.}
```

This identity is exact inside the specified model.

Using

```math
\Delta G=(1-a-b)\Delta c+\mathcal H(a)-\mathcal H(b),
```

and considering balanced learned experience, $a\approx b\approx\bar a$, gives the approximation

```math
\log\frac{P(\pi_2)}{P(\pi_1)}
\approx \gamma(1-2\bar a)\Delta c.
```

Thus, in this restricted regime, policy log odds are proportional to $\gamma\Delta c$ up to a state-dependent scaling factor. Away from balanced experience, the ambiguity difference need not vanish.

This is the main mathematical result that derives from the active-inference model used here. The question then is how can this log-ratio measure be approximated by measurable quantities (e.g. proxies from questionnaire items or structured interviews/observations). In a previous version of this work we used psychological distress and mental wellbeing as two such measures - however, these must be used with care, and in the next subsection we explain what such measures cannot do (before we move on to what we used them for).

## 3. Why questionnaire ratios do not inherit the identity

Psychological distress and wellbeing scales are neither policy probabilities nor complementary quantities. They also do not generally share a natural ratio scale. If

```math
D=f(P(\pi_2)), \qquad W=g(P(\pi_1)),
```

for unknown monotone measurement functions $f$ and $g$, then

```math
\log(D/W)
```

is not generally equal or proportional to

```math
\log\{P(\pi_2)/P(\pi_1)\}.
```

Unknown monotone transformations do not preserve ratios and need not preserve between-person ordering of ratios. Scale origins, floor corrections, and units can change the empirical value.

Consequently, a questionnaire distress-to-wellbeing log ratio does **not** recover $\gamma\Delta c$, does not locate an individual on a phase diagram, and is not a validated transition coordinate.

## 4. Role of the empirical balance analyses

The MCS and ABCD balance composites are retained as descriptive measurement analyses. They ask whether a prespecified nonlinear combination changes in-sample fit relative to its components.

The relevant conclusions are narrow:

- modest AIC differences can occur;
- changes in discrimination and calibration are very small;
- the analyses are not external validation; and
- no result identifies a latent active-inference parameter or phase coordinate.

Any future validation would require explicit measurement models, independent samples, preregistered transformations, and intensive repeated observations aligned to participant-level transition candidates.

In addition, the policy log odds may be a useful model-internal diagnostic. That does not make a questionnaire ratio a privileged early-warning variable. Early-warning behaviour depends on mechanism, sampling rate, noise, detrending, and the measured variable. Direct tests should compare candidate observables within fitted participant-level models rather than assuming a ratio is theoretically guaranteed to work.

## 5. Reproducibility map

- `Python_code/core_functions.py`: policy posterior and EFE difference.
- `R/balance_index.R`: exploratory cohort ratio construction and comparison.
- `SOM_Empirical_Results.Rmd`: full results and measurement caveats.
- `outputs/tables/mcs_balance_*.csv` and `outputs/tables/abcd_balance_*.csv`: aggregate results.

