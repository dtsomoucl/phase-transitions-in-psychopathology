# Notebook 02 — Application to Adolescent Psychopathology

**Series:** *Quantifying Phase Transitions in Psychopathology*  
**Date:** April 2026  
**Status:** Clinical translation and simulation results  
**Author:** Dr Dimitris I. Tsomokos  
Psychology & Human Development, Institute of Education, University College London

---

## 1. Purpose

Notebook 01 derived the core mathematical result in domain-general terms: a minimal (two-state) active-inference agent with Dirichlet learning exhibits phase transitions governed by a self-consistency equation that maps onto the Ising mean-field equation. This notebook translates that result into the specific domain of **mental health and psychopathology**, showing how the same mechanism can model the onset of, persistence in, and recovery from states resembling depressive disengagement and social withdrawal. Due to the subsequent application and testing of some of these ideas in large, longitudinal studies (Millennium Cohort Study and Adolescent Brain Cognitive Development Study), we will focus here on **adolescent** mental health.

The computational results summarised here are produced by the simulation codebase in the project repository, orchestrated via `main.py`, which pulls in `sim_psychopathology.py`, `step2_hysteresis.py`, `step2_robustness.py`, `step3_orthogonal.py`, `step4_decay.py`, `step5_asymmetric_memory.py`, and other scripts, and builds on the `core_functions.py` module reconstructed from the companion codebase (Tsomokos, 2026).

---

## 2. Clinical Reinterpretation of the Model

### 2.1 Mapping table

The two-state model from Notebook 01 is reinterpreted as follows:

| Model quantity | Clinical interpretation |
|---------------|------------------------|
| Policy $\pi_1$ | **Engagement**: active participation in social, academic, and goal-directed activities |
| Policy $\pi_2$ | **Withdrawal / avoidance**: reduced participation, social isolation, behavioural disengagement |
| Order parameter $z$ | Commitment balance: $z > 0$ = engagement-dominated; $z < 0$ = withdrawal-dominated |
| $P(\pi_1)$ | Latent probability of choosing engagement on a given occasion |
| Preference asymmetry $\Delta c$ | **Motivational field**: the subjective value the individual places on engagement relative to avoidance. Positive $\Delta c$ means engagement is experienced as more rewarding or meaningful |
| Policy precision $\gamma$ | **Cognitive self-efficacy / decisiveness**: the degree to which the individual can translate their beliefs into consistent, goal-directed action |
| Coupling function $\mathcal{G}(\tau; p)$ | Strength of the learning–action feedback at a given point in the individual's experience |
| Environmental discriminability $p$ | Clarity of environmental feedback: how distinguishable the consequences of engagement are from those of withdrawal |
| Prior strength $\alpha_0$ | Rigidity of pre-existing beliefs about what each behavioural strategy produces |

### 2.2 The learning–action coupling in clinical context

The positive-feedback mechanism identified in Notebook 01 acquires a direct clinical meaning. An adolescent who engages with their social and school environment samples observations that tend to reinforce the value of engagement — positive social feedback, academic achievement, reward. Conversely, an adolescent who withdraws samples a narrower, more impoverished set of observations, and the evidence gathered from that restricted sample makes withdrawal appear increasingly justified, because the individual no longer encounters the positive feedback that would support re-engagement.

This is the same self-reinforcing loop described in Notebook 01: the behavioural strategy determines the evidence encountered, and the evidence shapes the next behavioural choice. When this loop operates under favourable conditions (high $\Delta c$, adequate $\gamma$), it stabilises adaptive engagement. When external conditions deteriorate — through chronic stress, social rejection, peer victimisation, or loss — the same loop can lock the individual into a maladaptive withdrawn state.

### 2.3 Connection to the Behavioural Activation / Inhibition System framework

The reinterpretation of $\pi_1$ and $\pi_2$ as engagement versus withdrawal echoes a well-established tradition in personality and clinical psychology. Gray's **Reinforcement Sensitivity Theory** (RST) proposes that behaviour is governed by two broad motivational systems: the **Behavioural Activation System (BAS)**, sensitive to signals of reward and driving approach behaviour, and the **Behavioural Inhibition System (BIS)**, sensitive to signals of punishment and conflict and producing cautious withdrawal (Gray, 1982). (This is a brief description, simply to make the connection here - but obviously we are ignoring the revised RST and its third system, the Fight/Flight/Freeze that handles direct threat).

This framework maps onto our model in a natural way. The preference parameter $\Delta c$ reflects the relative strength of approach motivation over avoidance motivation — essentially, the balance of BAS over BIS influence on the individual's behavioural repertoire. High $\Delta c$ corresponds to strong approach motivation (dominant BAS); low or negative $\Delta c$ corresponds to a state in which avoidance and withdrawal are experienced as more valuable or safer than engagement. Bijttebier et al. (2009) reviewed the evidence linking BIS/BAS sensitivity to a range of psychiatric conditions and concluded that depression is characterised by a constellation of low BAS sensitivity (reduced reward responsiveness) combined with elevated BIS sensitivity, whilst anxiety disorders are more specifically associated with high BIS. Trew (2011) provided an integrative model of approach and avoidance processes in depression, arguing that depressive states involve a broad reduction in approach motivation alongside an increase in avoidance, with the resulting behavioural contraction maintaining the disorder. This pattern — low approach motivation plus high avoidance — is precisely the condition under which our model predicts a transition into the withdrawn state: low $\Delta c$ (reduced value of engagement) combined with the self-reinforcing withdrawal loop.

---

## 3. The Three-Phase Protocol

The simulations model an individual's trajectory through three distinct life-phases, implemented as a schedule of the two control parameters $\gamma$ and $\Delta c$.

### Phase A — Adaptive behaviour / healthy development (trials $0$ to $N_{\text{healthy}}$)

Both parameters are held at favourable values: $\gamma_0 = 16$ (high decisiveness) and $\Delta c_0 = 0.3$ (engagement is experienced as meaningfully more rewarding than withdrawal). The adolescent develops a stable engagement pattern, accumulating evidence that reinforces the value of participation. This corresponds to a period of adaptive functioning.

### Phase B — Adversity and prodromal deterioration (trials $N_{\text{healthy}}$ onward)

External stressors may begin to erode both parameters simultaneously. The motivational field $\Delta c$ drifts downward at a slow, continuous rate (default: $0.003$ per trial), reflecting the gradual loss of perceived value in engagement — perhaps through chronic stress from social rejection, family conflict, academic failure, etc. Cognitive precision $\gamma$ also declines (default rate: $0.005$ per trial), representing the erosion of self-efficacy and decisional consistency under sustained stress. Both parameters have floor values ($\Delta c_{\text{floor}} = -0.15$, $\gamma_{\text{floor}} = 14.0$), representing a chronic adverse state. (These are the default onset-regime values; the recovery analyses use a distinct parameter regime — see the note before Section 4.)

The critical feature of this phase is that the parameter drift is **slow and continuous**, but the resulting behavioural transition can be **abrupt**: once the system crosses the bifurcation boundary in the $(\gamma\mathcal{G}, \Delta c)$ control plane, the engaged attractor destabilises and the adolescent transitions rapidly into the withdrawn state. This is the onset of psychopathology as a phase transition.

### Phase C — Disordered state consolidation

Once $\gamma$ and $\Delta c$ reach their floor values, the adolescent is in a stable withdrawn state. The accumulated Dirichlet evidence now reflects extensive experience with withdrawal-related observations and relatively less recent engagement-related experience. This is the "stuck" state from which recovery must be attempted.

### Catastrophe flags in the clinical direction

The simulations confirm that the reverse transition (from engagement to withdrawal) reproduces all five catastrophe-flag signatures identified in Notebook 01, now in the clinically relevant direction: bimodality in the population distribution of $P(\pi_1)$ (some individuals remain engaged while others have transitioned); an inaccessible region at intermediate commitment levels; sudden jumps in individual trajectories; hysteresis when parameters are swept forward and then backward; and divergence across initial conditions.

### A note on parameter regimes

The onset simulations in Section 3 use the **onset regime** ($\alpha_0 = 40$, $\gamma_{\text{rate}} = 0.005$, $\Delta c_{\text{rate}} = 0.003$, $\Delta c_{\text{healthy}} = 0.3$, $\gamma_{\text{floor}} = 14.0$), which represents a strongly committed individual with rigid priors whose parameters erode slowly. The recovery and orthogonal-intervention analyses in Sections 4–5 use the **recovery regime** ($\alpha_0 = 2$, $\gamma_{\text{rate}} = 0.05$, $\Delta c_{\text{rate}} = 0.00075$, $\Delta c_{\text{healthy}} = 0.15$, $\gamma_{\text{floor}} = 5.0$), which represents a lightly committed individual with weak priors and faster parameter drift — a setting in which recovery is not trivially guaranteed and the relative contributions of $\gamma$ and $\Delta c$ can be cleanly distinguished. The robustness checks in Section 4.3 confirm that the core field-dominance result holds across both regimes as well as under a cross-regime condition that applies onset-level priors ($\alpha_0 = 40$) with recovery-level drift rates.

---

## 4. Recovery and Path Asymmetry

### 4.1 The central question

Wang et al. (2026) argued that most studies on psychopathology specify how individuals transition into disordered states but not out of them, implicitly assuming a symmetry between onset and recovery that complex-systems theory does not support. The present framework addresses this directly by asking: once the adolescent is in the consolidated "illness" state, which parameters must be restored, in what combination, and does the recovery path simply retrace the onset path?

### 4.2 The recovery boundary

The recovery analysis (`step2_hysteresis.py`) maps the outcome of attempted recovery across the two-dimensional control space $(\gamma_{\text{restored}}, \Delta c_{\text{restored}})$. At each point in this space, the adolescent's control parameters are set to the restored values after a period of illness consolidation, and the adolescent is allowed to re-learn under the new parameter regime. Recovery is defined as the adolescent's late-window $P(\pi_1)$ exceeding a threshold (default: $0.5$).

The recovery boundary — the curve in $(\gamma, \Delta c)$ space that separates recovery from persistent illness — reveals a stark **asymmetry between the two control parameters**:

- **Restoring $\Delta c$ alone** (returning the motivational field to its recovery-regime healthy value of $0.15$ while leaving $\gamma$ at its illness-phase floor) produces robust recovery. Even with low cognitive precision, a sufficiently positive motivational field drives the adolescent back toward engagement.
- **Restoring $\gamma$ alone** (returning cognitive precision to its healthy value while leaving $\Delta c$ at its adverse floor) **fails to produce recovery**. High decisiveness applied to an adverse motivational landscape simply makes the adolescent more decisively withdrawn.

This is the **dominance of motivational field** result ("field-dominance" result, for simplicity): recovery is driven primarily by the restoration of the motivational field $\Delta c$, not by the restoration of cognitive precision $\gamma$.

### 4.3 Robustness

The field-dominance result is not a fragile feature of one parameter setting. The robustness checks (`step2_robustness.py`) confirm that the ordering — $\Delta c$ matters more than $\gamma$ for recovery — holds across:

- Different recovery thresholds ($P(\pi_1) > 0.5, 0.6, 0.7$).
- Different environmental discriminabilities ($p = 0.75$ to $0.95$).
- Different prior strengths ($\alpha_0 = 1$ to $4$).

A cross-regime condition (`ONSET_AS_RECOVERY_REGIME`: onset-level priors $\alpha_0 = 40$ combined with recovery-level drift rates) further confirms the ordering. In every regime tested, the difference in late-window $P(\pi_1)$ between the $\Delta c$-restored and $\gamma$-restored conditions is positive: restoring the field outperforms restoring precision.

### 4.4 Weak evidence-driven hysteresis

The recovery analysis also documents a **substantive negative result**. In the minimal two-state model with permanent (non-decaying) Dirichlet evidence, the hysteresis driven by accumulated maladaptive evidence is **weak** — of order $O(0.001)$ in the ambiguity difference once both strategies have been explored for more than approximately 50 trials. The dominant barrier to recovery is the ongoing adverse preference field, not accumulated maladaptive memory. This finding constrains what the minimal model can and cannot explain: genuine "scar" effects from prolonged illness would require additional mechanisms beyond what this two-state symmetric model provides (see Section 6 below).

---

## 5. Orthogonal Interventions

### 5.1 Clinical motivation

Obviously, increasing motivation is not easy and psychotherapeutic interventions that do this take time. However, if $\Delta c$ can only be partially restored — to a neutral value ($\Delta c = 0$) rather than a "completely healthy" positive value (so to speak) — can other interventions close the remaining gap?

The orthogonal-intervention analysis (`step3_orthogonal.py`) tests two candidates:

1. **Count tempering** (parameter: $\eta$): rescaling all Dirichlet concentration parameters toward their prior baseline, effectively "softening" accumulated maladaptive evidence. This corresponds loosely to interventions that aim to weaken the hold of ingrained negative beliefs or expectations — along the lines of cognitive restructuring and schema therapy.

2. **Environmental enrichment** (parameter: $p_{\text{new}} > p$): increasing the discriminability of the environment during recovery, so that the consequences of engagement become more clearly distinguishable from those of withdrawal. This corresponds to structured therapeutic environments that provide clearer, more immediate feedback for participation.

### 5.2 Results

Under partial preference restoration ($\Delta c$ restored to $0$ rather than the recovery-regime healthy value of $0.15$):

- **Restoring $\gamma$** alone produces a negligible additional benefit beyond partial $\Delta c$ restoration.
- **Count tempering** ($\eta = 0.1$) provides a small additional benefit when $\Delta c$ is already favourable, but it can be counterproductive when $\Delta c$ remains adverse, because softening evidence removes the agent's memory of engagement-related observations along with its withdrawal-related ones.
- **Environmental enrichment** ($p = 0.95$) **amplifies whichever direction the field favours**: it helps when $\Delta c > 0$ (the clearer feedback reinforces re-engagement) but is counterproductive when $\Delta c < 0$ (the clearer feedback reinforces withdrawal more decisively).
- **Full preference restoration** ($\Delta c$ returned to the recovery-regime healthy value of $0.15$) remains the only condition that reliably produces robust recovery.

The key takeaway from all this is that it is possible that **no orthogonal intervention substitutes for restoring the motivational field** (i.e., further research is needed to explore this possibility further). This may a direct consequence of the field-dominance structure identified in Section 4: because the recovery boundary is organised primarily around $\Delta c$, interventions that operate on other parameters cannot compensate for an unremediated motivational deficit.

---

## 6. Symmetric Evidence Decay: A Structural Negative Result

### 6.1 Motivation

A natural extension of the model is to introduce **evidence decay**: at each trial, all Dirichlet concentration parameters are pulled back toward their prior baseline at rate $\eta$:

```math
\alpha_i \leftarrow \eta\,\alpha_i + (1 - \eta)\,\alpha_0/2
```

This makes recent evidence count more than distant evidence. The original hypothesis was that an agent who has been withdrawn for a long time would gradually forget the evidence supporting engagement, creating an asymmetry in learned discriminabilities that resists recovery even after full parameter reversal — a form of illness-residence-time ($T_{\text{ill}}$) dependent hysteresis.

### 6.2 Result

The hypothesis is not supported. Symmetric evidence decay (`step4_decay.py`) does **not** produce $T_{\text{ill}}$-dependent hysteresis under full same-variable reversal ($\gamma$ and $\Delta c$ both restored to healthy values). The reason is structural: because the decay applies identically to all concentration parameters, it erases illness-phase evidence during recovery just as effectively as it erases health-phase evidence during illness. Under full reversal with adequate recovery time, faster decay can actually *help* recovery by making the adolescent more responsive to the restored positive preference field.

### 6.3 Implications

Actually, this is a meaningful negative result. It establishes a **structural boundary** on what the minimal two-state model with symmetric Dirichlet learning can explain. Genuine asymmetric hysteresis — where prolonged disorder/illness makes recovery progressively harder even after the external conditions are restored — would require mechanisms that break the symmetry between the two strategies at the level of evidence accumulation. Candidates include strategy-dependent decay rates (withdrawal-related evidence is retained more strongly than engagement-related evidence), domain-specific forgetting, or structural asymmetries between the two behavioural modes. The first of these candidates — strategy-specific retention — is explored in Section 7.

---

## 7. Strategy-Specific Evidence Retention (Asymmetric Memory)

### 7.1 Motivation

Section 6 showed that **symmetric** evidence decay cannot produce illness-duration-dependent hysteresis. This section introduces the minimal asymmetry that can: **strategy-specific retention rates**. In each trial, the Dirichlet concentration parameters associated with each strategy are decayed toward their prior baseline at strategy-specific rates $\eta_{\text{eng}}$ and $\eta_{\text{wdr}}$, where $\eta = 1$ means perfect retention and lower values mean faster forgetting:

```math
\alpha_i^{(s)} \leftarrow \eta_s\,\alpha_i^{(s)} + (1 - \eta_s)\,\alpha_0/2
```

The hypothesis is that if withdrawal-related evidence is retained more strongly than engagement-related evidence ($\eta_{\text{wdr}} > \eta_{\text{eng}}$), then an individual who has been ill for longer will have accumulated a larger asymmetry in learned discriminabilities — and this asymmetry will resist recovery even after full parameter restoration.

### 7.2 Results

The simulations (`step5_asymmetric_memory.py`) confirm the hypothesis. When $\eta_{\text{eng}} < \eta_{\text{wdr}}$, recovery success under full parameter reversal declines monotonically with illness duration $T_{\text{ill}}$. The effect is absent when $\eta_{\text{eng}} = \eta_{\text{wdr}}$ (the symmetric case from Section 6) and strengthens as the gap between the two retention rates widens.

The per-trial retention rates used (e.g. $\eta_{\text{eng}} = 0.985$) may appear close to $1$, but their cumulative effect over long illness durations is substantial: $0.985^{200} \approx 0.05$, meaning that after 200 trials only about 5% of the original engagement-related evidence is retained. Lower per-trial values would create floor effects in which virtually all evidence is erased within a few dozen trials, removing the graduated relationship between $T_{\text{ill}}$ and recovery difficulty.

### 7.3 Implications

This result resolves the structural limitation identified in Section 6. The minimal two-state model can produce genuine illness-duration-dependent path asymmetry — where prolonged disorder makes recovery progressively harder — provided that the evidence-retention mechanism is strategy-specific rather than symmetric. Psychologically, this asymmetry is plausible: withdrawal restricts the range of experiences encountered, and the narrower behavioural repertoire during illness may be more internally consistent (and therefore more strongly consolidated) than the broader repertoire experienced during health. The result aligns with clinical observations that prolonged depressive episodes are associated with poorer treatment response and higher relapse risk.

---

## 8. Summary of Simulation Results

| Finding | Status | Implication |
|---------|--------|-------------|
| Onset reproduces all five catastrophe flags in the clinical direction | Confirmed across default and varied parameters | The phase-transition framework applies to the engagement → withdrawal transition |
| Recovery is field-dominated: restoring $\Delta c$ matters more than restoring $\gamma$ | Robust across thresholds, $p$, and $\alpha_0$ | Motivation/re-engagement is the primary driver of recovery; cognitive precision plays a secondary role |
| Orthogonal interventions cannot substitute for field restoration | Confirmed under partial $\Delta c$ restoration | Clinical interventions that do not address motivation directly have limited leverage |
| Symmetric evidence decay does not produce $T_{\text{ill}}$-dependent hysteresis | Structural negative result | The minimal model cannot explain illness-duration scarring without additional asymmetric mechanisms |
| Strategy-specific retention ($\eta_{\text{eng}} < \eta_{\text{wdr}}$) produces $T_{\text{ill}}$-dependent path asymmetry | Confirmed; effect scales with retention-rate gap | Asymmetric memory is the minimal mechanism needed for genuine illness-duration-dependent hysteresis |

---

## 9. Bridging to the Empirical Analyses

The dominance of the preference field result generates a testable prediction for longitudinal cohort data: **motivation and engagement indices should be stronger prospective predictors of internalising trajectories than cognitive-control indices**, once baseline distress and standard covariates are taken into account. In the language of the model:

- **Motivation/engagement measures** (reward responsiveness, school engagement, approach behaviour, sense of value in participation) are observational proxies for the preference field $\Delta c$.
- **Executive-control measures** (working memory, inhibitory control, decision-making consistency) are observational proxies for policy precision $\gamma$.

This prediction is tested directly in the MCS and ABCD cohort analyses documented in the Supplementary Online Material (SOM), which can be found alongside these notebooks and in the R code folder. As an additional check on the alignment between the simulation and the empirical results, a synthetic internal-consistency analysis (`empirical_bridge.py`) confirms that a field-dominant data-generating process reproduces the empirical coefficient ordering and the balance-index advantage documented in Notebook 03.

Finally, the next notebook (Notebook 03) develops an additional theoretical prediction concerning **balance indices** — ratios of distress to wellbeing measures — and their relationship to the product $\gamma\Delta c$.

---

## References

Bijttebier, P., Beck, I., Claes, L., & Vandereycken, W. (2009). Gray's Reinforcement Sensitivity Theory as a framework for research on personality–psychopathology associations. *Clinical Psychology Review*, 29(5), 421–430. https://doi.org/10.1016/j.cpr.2009.04.002

Gray, J. A. (1982). Précis of The neuropsychology of anxiety: An enquiry into the functions of the septo-hippocampal system. *Behavioral and Brain Sciences*, 5(3), 469-484. doi:10.1017/S0140525X00013066 

Trew, J. L. (2011). Exploring the roles of approach and avoidance in depression: An integrative model. *Clinical Psychology Review*, 31(7), 1156–1168. https://doi.org/10.1016/j.cpr.2011.07.007

Tsomokos, D. I. (2026). *From Incremental Learning to Developmental Stage Shifts: Phase transitions in active inference.* Code and notebooks: [https://github.com/dtsomoucl/phase-transitions-in-active-inference](https://github.com/dtsomoucl/phase-transitions-in-active-inference).

Wang SB, Blanken TF, van der Maas HLJ, Borsboom D. Path Asymmetry in Complex Dynamic Systems of Psychopathology. *JAMA Psychiatry.* 2026;83(1):99–100. doi:10.1001/jamapsychiatry.2025.3147
