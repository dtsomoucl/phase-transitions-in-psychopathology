# Notebook 03 — The Balance Index: Ratios as Phase-Space Observables

**Series:** *Quantifying Phase Transitions in Psychopathology*  
**Date:** April 2026  
**Status:** Analytical derivation and empirical operationalisation  
**Author:** Dr Dimitris I. Tsomokos  
Psychology & Human Development, Institute of Education, University College London

---

## 1. Purpose

Notebooks 01 and 02 established that the onset of and recovery from psychopathological states can be modelled as phase transitions governed by two control parameters: the motivational preference field $\Delta c$ and the cognitive precision $\gamma$. A natural question arises: **can we track an individual's position in this two-dimensional control space using observable measures?**

This notebook shows that the answer is a qualified *yes*. The **log-ratio** of the policy posterior — or equivalently, a log-ratio of observable proxies such as psychological distress and mental wellbeing — recovers the **product** $\gamma \Delta c$ up to a known scaling factor. This is informative: it tells us where the individual sits along a curve in the $(\gamma, \Delta c)$ phase space. But it cannot decompose the product into its two factors. We derive this result analytically, explain its empirical operationalisation, and discuss what the balance index can and cannot reveal.

---

## 2. The Log-Odds of the Policy Posterior

### 2.1 Starting point

From Notebook 01, the posterior probability of the engagement policy $\pi_1$ is:

```math
P(\pi_1) = \sigma(-\gamma\,\Delta G) = \frac{1}{1 + \exp(\gamma\,\Delta G)}
```

and the posterior probability of the withdrawal policy $\pi_2$ is:

```math
P(\pi_2) = 1 - P(\pi_1) = \frac{\exp(\gamma\,\Delta G)}{1 + \exp(\gamma\,\Delta G)}
```

### 2.2 The exact log-odds

The log-ratio of the two policy posteriors is:

```math
\log\frac{P(\pi_2)}{P(\pi_1)} = \gamma\,\Delta G
```

This is exact — no approximation is involved. Since $\Delta G = G(\pi_1) - G(\pi_2)$ is the difference in expected free energy between the two policies, the log-odds is simply $\gamma$ times the EFE difference. When $\Delta G > 0$ (engagement is more costly in free-energy terms), the log-ratio is positive and the agent favours withdrawal; when $\Delta G < 0$, the agent favours engagement.

### 2.3 Expanding the EFE difference

Recall from Notebook 01 that:

```math
\Delta G = (1 - a - b)\,\Delta c + \mathcal{H}(a) - \mathcal{H}(b)
```

where $a$ and $b$ are the learned discriminabilities and $\mathcal{H}$ is binary entropy. The first term is the **pragmatic** (preference-driven) component; the second is the **ambiguity** (information-driven) component.

---

## 3. The Balanced-Experience Approximation

### 3.1 Near the symmetric point

When the agent has roughly balanced experience with both strategies ($z \approx 0$, i.e. $\phi \approx 1/2$), the learned discriminabilities are approximately equal: $a \approx b \approx \bar{a}(\tau)$, where

```math
\bar{a}(\tau) = \frac{1}{2} + \left(p - \frac{1}{2}\right)\frac{\tau}{1+\tau}
```

is the mean discriminability at the symmetric point (Section 6 of Notebook 01). Under this condition, the ambiguity difference vanishes:

```math
\mathcal{H}(a) - \mathcal{H}(b) \approx 0
```

and the EFE difference reduces to its pragmatic component:

```math
\Delta G \approx (1 - 2\bar{a})\,\Delta c
```

### 3.2 The log-odds as a product

Substituting into the exact log-odds expression:

```math
\log\frac{P(\pi_2)}{P(\pi_1)} \approx \gamma\,(1 - 2\bar{a})\,\Delta c
```

Since $\bar{a} > 1/2$, the factor $(1 - 2\bar{a})$ is negative, so we can write:

```math
\boxed{\log\frac{P(\pi_2)}{P(\pi_1)} \approx -\underbrace{(2\bar{a} - 1)}_{\text{positive, depends on } \tau, p}\;\cdot\;\gamma\,\Delta c}
```

The key insight is that the log-odds is proportional to the **product** $\gamma\Delta c$, scaled by a factor $(2\bar{a} - 1)$ that depends only on the agent's developmental time $\tau$ and the environmental discriminability $p$, but **not** on $\gamma$ or $\Delta c$ themselves.

### 3.3 What the product tells us

The product $\gamma\Delta c$ has a clear geometric meaning: it locates the individual along a **hyperbolic curve** in the $(\gamma, \Delta c)$ control plane. Two individuals with the same value of $\gamma\Delta c$ — say, one with high precision and moderate motivation, the other with moderate precision and high motivation — will have the same log-odds and the same balance-index value. The ratio cannot distinguish between them.

This is an inherent limitation of any single observable that depends on both control parameters only through their product. It is not a failure of the measurement but a structural feature of the phase-space geometry.

---

## 4. Empirical Operationalisation

### 4.1 From latent posteriors to observable measures

In the psychopathology framing of Notebook 02:

- $P(\pi_1)$, the probability of engagement, should be reflected in **wellbeing-related observables**: positive mental health, life satisfaction, social participation, reward responsiveness.
- $P(\pi_2)$, the probability of withdrawal, should be reflected in **distress-related observables**: psychological distress, depressive symptoms, anxiety, social withdrawal.

Neither distress nor wellbeing is a direct measurement of $P(\pi_1)$ or $P(\pi_2)$. They are **noisy, monotone proxies**: higher distress indicates higher $P(\pi_2)$, and higher wellbeing indicates higher $P(\pi_1)$, but the functional form of the mapping is unknown. Despite this, the log-ratio inherits the structural property of reflecting the product $\gamma\Delta c$, because any monotone transformation of the policy posteriors preserves the sign and ordering of the log-ratio (even if it distorts the exact magnitude).

### 4.2 The balance index $\psi$

We define the empirical balance index as:

```math
\psi = \log\!\left(\frac{\text{distress}}{\text{wellbeing}}\right)
```

In the MCS data, the natural operationalisation at age 17 is:

```math
\psi_{17} = \log\!\left(\frac{\text{Kessler-6 psychological distress}}{\text{WEMWBS mental wellbeing}}\right)
```

Both instruments are measured on positive scales at the same sweep, so the ratio is well-defined for all individuals with non-zero scores on both measures.

### 4.3 Properties of the balance index

Several properties follow from the derivation:

1. **Sign convention.** Positive $\psi$ indicates that distress exceeds wellbeing (the agent is in the withdrawal-dominated regime); negative $\psi$ indicates the reverse.

2. **Theoretical grounding.** The balance index is not an arbitrary ratio. It is the empirical counterpart of the log-odds of the policy posterior, which the theory shows to be proportional to the product $\gamma\Delta c$ — the combined "location" in the phase space that governs transitions.

3. **Sensitivity to transitions.** Near the bifurcation boundary ($\gamma\mathcal{G} \approx 1$), the log-odds is maximally sensitive to small changes in $\Delta c$ (because the susceptibility diverges). The balance index should therefore be most informative — and most variable across individuals — at developmental moments near the critical window.

4. **Superiority over components alone.** Distress and wellbeing measured separately each carry partial information: distress alone cannot distinguish between an individual with low engagement and one with high distress tolerance, and wellbeing alone cannot distinguish between an individual with genuine positive engagement and one who simply avoids distressing contexts. The ratio captures the **relative balance** between the two, which is what the theory predicts to be the informative quantity.

---

## 5. What the Balance Index Cannot Tell Us

### 5.1 The decomposition problem

The balance index recovers $\gamma\Delta c$ but cannot decompose the product into its two factors. An individual with $\psi = 0$ (balanced distress and wellbeing) might be:

- in a genuinely healthy state with moderate $\gamma > 0$ and moderate $\Delta c > 0$, where the balance reflects stable engagement, or
- in an indeterminate state with high $\gamma$ but $\Delta c \approx 0$, where the agent is cognitively decisive but motivationally neutral, or
- at the critical point with $\gamma\mathcal{G} = 1$ and $\Delta c = 0$, poised for a bifurcation in either direction.

A single cross-sectional measurement of $\psi$ cannot distinguish among these possibilities. This is a fundamental identifiability limitation.

### 5.2 Longitudinal resolution

Longitudinal trajectories of $\psi$ carry additional information that is not available in a single snapshot. Near the critical point, the theory predicts **critical slowing down**: the system's return time to equilibrium after a perturbation increases, and the autocorrelation and variance of the order parameter grow. In empirical terms:

- If $\psi$ is measured repeatedly over time, an increase in its **temporal autocorrelation** (the balance index at time $t$ becoming more predictive of the balance index at time $t+1$) would signal proximity to the critical boundary.
- An increase in the **variance** of $\psi$ across individuals at a given developmental time point would signal that the population is near or within the bifurcation window.

These critical-slowing-down signatures are, in principle, exploitable as **early-warning indicators** of an impending transition. The balance index is the natural observable in which to look for them, because the theory predicts that $\psi$ tracks the relevant phase-space coordinate.

### 5.3 The role of separate predictors

The decomposition problem is precisely why the empirical analyses in the accompanying SOM do not rely on the balance index alone. The separate motivation/engagement and executive-control predictor families provide complementary information that the ratio cannot. The **field-dominance** result (Notebook 02) predicts that motivation/engagement predictors ($\Delta c$ proxies) should outperform executive-control predictors ($\gamma$ proxies) as prospective predictors of internalising trajectories. The balance index predicts that the **ratio** of distress to wellbeing should track these trajectories more effectively than either component alone. These are distinct, non-redundant predictions, and the empirical pipeline tests both.

---

## 6. Beyond the Symmetric Approximation

### 6.1 When is the approximation valid?

The balanced-experience approximation ($a \approx b$, ambiguity difference $\approx 0$) holds when the agent has not committed too strongly to one strategy — i.e., when $|z|$ is not too large. For individuals deep in the withdrawal regime ($z \ll 0$), the ambiguity difference can become non-negligible: the agent's evidence about the engagement state may be stale (low $n_1$, hence $a$ near 1/2 and $\mathcal{H}(a)$ near $\ln 2$), while evidence about the withdrawal state is rich ($b$ near $p$ and $\mathcal{H}(b)$ lower). In this regime, the full expression

```math
\log\frac{P(\pi_2)}{P(\pi_1)} = \gamma\left[(1 - a - b)\,\Delta c + \mathcal{H}(a) - \mathcal{H}(b)\right]
```

has a non-trivial ambiguity contribution that shifts the log-odds relative to the simple $\gamma\Delta c$ prediction.

### 6.2 Implications for measurement

In practice, the individuals for whom the approximation breaks down most severely are those who are most deeply withdrawn — precisely the clinical population of greatest interest. For these individuals, the balance index still tracks the correct *direction* (distress-dominated implies withdrawal-dominated), but the precise quantitative relationship to $\gamma\Delta c$ is distorted by the ambiguity asymmetry. This is a reason to use the balance index as an **ordinal** indicator of phase-space position rather than as a precise **cardinal** estimate of $\gamma\Delta c$.

---

## 7. Summary

| Claim | Derivation status |
|-------|-------------------|
| $\log(P(\pi_2)/P(\pi_1)) = \gamma\Delta G$ | Exact |
| Near balanced experience, $\gamma\Delta G \approx -(2\bar{a}-1)\gamma\Delta c$ | Approximation, valid when $\|z\| \ll 1$ |
| The balance index $\psi = \log(\text{distress}/\text{wellbeing})$ is a noisy proxy for $\gamma\Delta c$ | Empirical translation (monotone proxy assumption) |
| $\psi$ recovers the product $\gamma\Delta c$ but cannot decompose it | Structural identifiability limitation |
| Longitudinal $\psi$ trajectories carry additional information via critical-slowing-down signatures | Theoretical prediction, testable |
| $\psi$ should track internalising trajectories better than either distress or wellbeing alone | Testable prediction derived from the theory |

The balance index provides a principled, theory-derived observable for tracking an individual's position in the phase space that governs transitions into and out of psychopathological states. Its power lies in capturing the *relative balance* between engagement and withdrawal — the quantity that the active-inference framework identifies as the relevant dynamical variable. Its limitation is that it cannot, from a single measurement, distinguish between the two control parameters that jointly determine that balance. The empirical analyses in the SOM test whether $\psi$ constructed from MCS Kessler-6 and WEMWBS measures at age 17 predicts age-23 outcomes above and beyond what the individual components predict separately.

---

## References

Dakos, V., Scheffer, M., van Nes, E. H., Brovkin, V., Petoukhov, V., & Held, H. (2008). Slowing down as an early warning signal for abrupt climate change. *Proceedings of the National Academy of Sciences*, 105(38), 14308–14312.

Scheffer, M., Bascompte, J., Brock, W. A., Brovkin, V., Carpenter, S. R., Dakos, V., ... & Sugihara, G. (2009). Early-warning signals for critical transitions. *Nature*, 461(7260), 53–59.

Tsomokos, D. I. (2026). *From Incremental Learning to Developmental Stage Shifts: Phase transitions in active inference.* Code and notebooks: [https://github.com/dtsomoucl/phase-transitions-in-active-inference](https://github.com/dtsomoucl/phase-transitions-in-active-inference).

Wichers, M., Groot, P. C., & Psychosystems, ESM Group, & EWS Group. (2016). Critical slowing down as a personalized early warning signal for depression. *Psychotherapy and Psychosomatics*, 85(2), 114–116.
