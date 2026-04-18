# Notebook 01 — Phase Transitions in a Simple Active-Inference Model

**Series:** *Quantifying Phase Transitions in Psychopathology*  
**Date:** March-April 2026  
**Status:** Analytical derivation of the main result  
**Author:** Dr Dimitris I. Tsomokos  
Psychology & Human Development, Institute of Education, University College London

---

## 1. Purpose

This notebook derives the core mathematical result for a paper on *phase transitions in psychopathology and path asymmetry between disorder onset and recovery*. Starting from a minimal active-inference model, we show that the interplay between an agent's policy selection and its parameter learning produces a **self-consistency equation** that is formally identical to the **Curie–Weiss mean-field equation** from statistical physics. This equation admits bifurcations — abrupt, discontinuous transitions in behavioural commitment — and its local normal form is a **cusp catastrophe**. The model used here is absolutely minimal, and purposefully so as we aim to understand the nature of such transitions in the simplest possible case (with two states only).

The derivation presented is deliberately concise and domain-neutral: it makes no reference to any specific applied domain (that 'translation' into what it means practically is provided in Notebook 02 of this repository, and of course in the paper too). Note that, for the **full step-by-step derivation** with algebraic intermediates, one may consult the companion repository (https://github.com/dtsomoucl/phase-transitions-in-active-inference) and its [Notebook 04](https://github.com/dtsomoucl/phase-transitions-in-active-inference/blob/main/notebooks/Notebook_04_Step_by_Step.md), which reported the mathematical derivations initially in a different context within a developmental-psychology focus and framing (Tsomokos, 2026). The core algebraic steps have also been verified in [Lean 4](https://github.com/leanprover) (see the `lean/` directory of the companion repository). From herein on when we refer to notebooks, we mean notebooks in this folder (of the present repository), unless otherwise stated.

### Levels of argument

Three levels of argument are used below and should be kept distinct:

1. **Exact algebra** for the expected free energy (EFE) difference in the specified two-state model.
2. **Mean-field closure**, in which stochastic Dirichlet updates are replaced by their expected counts, reducing the system to a deterministic self-consistency equation.
3. **Local normal-form analysis** near the primary bifurcation, recovering the cusp catastrophe.

---

## 2. Model Specification

We define the simplest partially observable Markov decision process (POMDP) that can exhibit strategy competition coupled with parameter learning.

| Component | Specification |
|-----------|---------------|
| Hidden states | $s \in \{s_1, s_2\}$ |
| Observations | $o \in \{o_1, o_2\}$ |
| Policies | $\pi_1$ (leads to $s_1$), $\pi_2$ (leads to $s_2$) |
| Preferences | $\mathbf{C} = (c_1, c_2)^\top$, with asymmetry $\Delta c = c_1 - c_2$ |

### Observation likelihood (A matrix)

The agent's learned observation model is:

```math
\mathbf{A} = \begin{pmatrix} a & 1-b \\ 1-a & b \end{pmatrix}
```

where $a = \hat{P}(o_1 \mid s_1)$ and $b = \hat{P}(o_2 \mid s_2)$ are the learned discriminabilities. Both start near $1/2$ (maximum uncertainty) and, through Dirichlet learning, approach their true values $a^* = b^* = p \in (1/2, 1)$.

### Dirichlet parameterisation

Each column of **A** is parameterised by Dirichlet concentration parameters. Let $n_1$ and $n_2$ denote the number of observations accumulated from states $s_1$ and $s_2$ respectively, and let $\alpha_0$ denote the prior concentration sum per column. The expected point estimate for the first column is:

```math
a = \frac{1}{2} + \left(p - \frac{1}{2}\right) \frac{n_1}{\alpha_0 + n_1}
```

and symmetrically for $b$ with $n_2$.

### Transition model

Each policy deterministically selects its target state: $\pi_1$ transitions the agent to $s_1$, and $\pi_2$ transitions to $s_2$, regardless of the current state.

---

## 3. Expected Free Energy and Policy Selection

### EFE difference

Under the deterministic transition model, the EFE for policy $\pi_k$ decomposes into a **pragmatic** (preference-related) term and an **ambiguity** (entropy) term. The quantity governing policy selection is the difference $\Delta G = G(\pi_1) - G(\pi_2)$, which evaluates to:

```math
\boxed{\Delta G = (1 - a - b) \, \Delta c + \mathcal{H}(a) - \mathcal{H}(b)}
```

where $\mathcal{H}(x) = -x \ln x - (1-x)\ln(1-x)$ is the binary entropy function. The first term captures the **pragmatic difference** (which policy better satisfies preferences), and the second captures the **ambiguity difference** (which policy leads to more informative observations).

### Policy posterior

The posterior probability of policy $\pi_1$ is a softmax over negative EFE:

```math
P(\pi_1) = \sigma(-\gamma \, \Delta G) = \frac{1}{1 + \exp(\gamma \, \Delta G)}
```

where $\gamma > 0$ is the **policy precision** (inverse temperature), governing how sharply the agent commits to the better-scoring policy.

---

## 4. The Learning–Action Coupling

This is the critical mechanism for the phase transition. The agent's policy choice determines which state it visits, and therefore which column of **A** receives a Dirichlet update:

- Selecting $\pi_1$ → visiting $s_1$ → updating column 1 → $n_1$ increases, $a$ sharpens toward $p$.
- Selecting $\pi_2$ → visiting $s_2$ → updating column 2 → $n_2$ increases, $b$ sharpens toward $p$.

This creates **positive feedback**: an agent that favours $\pi_1$ accumulates more evidence about $s_1$, reducing the ambiguity of $\pi_1$, which makes $\pi_1$ even more attractive. The feedback between learning and action selection is the engine of the phase transition.

---

## 5. Mean-Field Reduction

### Order parameter

Let $N = n_1 + n_2$ denote total accumulated evidence and $\phi = n_1/N$ the fraction allocated to $\pi_1$. We define the order parameter:

```math
z = 2\phi - 1 \in [-1, 1]
```

which measures the degree of asymmetry in experience: $z = 0$ is balanced; $|z| \to 1$ is full commitment to one policy.

### Learned discriminabilities in terms of z

Defining the rescaled time $\tau = N/(2\alpha_0)$:

```math
a(z) = \frac{1}{2} + \left(p - \frac{1}{2}\right)\frac{(1+z)\tau}{1 + (1+z)\tau}, \qquad b(z) = \frac{1}{2} + \left(p - \frac{1}{2}\right)\frac{(1-z)\tau}{1 + (1-z)\tau}
```

### The self-consistency equation

Replacing the stochastic Dirichlet updates by their expectations (the mean-field approximation), the equilibrium condition $\phi^* = P(\pi_1)$ becomes:

```math
\boxed{z = \tanh\!\left(-\frac{\gamma}{2}\,\Delta G(z, \tau, \Delta c)\right)}
```

This is **formally identical to the Curie–Weiss mean-field equation** $m = \tanh(\beta J m + \beta h)$ from statistical mechanics, where $m$ is the magnetisation, $\beta$ is inverse temperature, $J$ is the spin coupling, and $h$ is the external field.

---

## 6. Bifurcation Analysis

### The coupling function

In the symmetric case ($\Delta c = 0$), $z = 0$ is always a fixed point. Its stability is determined by the **coupling function**:

```math
\mathcal{G}(\tau; p) = \left(p - \frac{1}{2}\right)\frac{\tau}{(1+\tau)^2}\left|\ln\frac{1-\bar{a}(\tau)}{\bar{a}(\tau)}\right|
```

where $\bar{a}(\tau) = \frac{1}{2} + (p - \frac{1}{2})\frac{\tau}{1+\tau}$ is the mean learned discriminability at $z = 0$. This function starts at zero, rises to a maximum $\mathcal{G}_{\max}$ at an intermediate time $\tau_{\max}$, and returns to zero as learning saturates.

### Critical condition

The symmetric fixed point $z = 0$ loses stability — and the system undergoes a **pitchfork bifurcation** — when:

```math
\boxed{\gamma \cdot \mathcal{G}_{\max}(p) = 1}
```

Equivalently, the **critical precision** is:

```math
\gamma_c = \frac{1}{\mathcal{G}_{\max}(p)}
```

If $\gamma > \gamma_c$, the learning–action feedback is strong enough to produce a genuine bifurcation in behavioural commitment. The critical precision depends only on the environmental discriminability $p$; the prior strength $\alpha_0$ affects only the absolute timing $N^* = 2\alpha_0\tau^*$ but not whether the bifurcation occurs.

### Key properties

1. **Finite developmental window.** The coupling function has a bounded support in time: before the window, the agent has too little evidence for feedback to matter; after it, learning has saturated and feedback weakens.
2. **Conditions favouring bifurcation.** Higher $p$ (clearer environment) and higher $\gamma$ (more decisive agent) make abrupt transitions more likely. Lower $\alpha_0$ (weaker prior) makes them occur earlier.
3. **Conditions favouring smooth change.** Low $p$, low $\gamma$, or a very strong prior produce gradual shifts without discontinuity.

---

## 7. The Ising Mapping

The identification between the self-consistency equation and the Curie–Weiss equation is not merely a formal analogy — it is a **mathematical identity** at the level of the mean-field equation:

| Ising model | Active inference | Interpretation |
|-------------|-----------------|----------------|
| Magnetisation $m$ | Order parameter $z$ | Degree of behavioural commitment |
| Inverse temperature $\beta$ | Policy precision $\gamma/2$ | Cognitive decisiveness |
| Coupling constant $J$ | Coupling function $\mathcal{G}(\tau; p)$ | Learning–action feedback strength |
| External field $h$ | $-(1-2\bar{a})\Delta c/2$ | Preference asymmetry / motivational bias |
| Curie temperature $T_c$ | Critical condition $\gamma \mathcal{G}_{\max} = 1$ | Threshold for symmetry breaking |

Because the mapping is exact at this level, the active-inference system inherits several well-known properties of the Ising model. Near the critical point, strategic commitment grows as $|z| \sim (\gamma\mathcal{G} - 1)^{1/2}$, population variance scales as $\text{Var}(z) \sim (\gamma\mathcal{G} - 1)^{-1}$, and the response to preference asymmetry follows $z \sim \Delta c^{1/3}$ at the critical point.

---

## 8. Local Cusp-Catastrophe Structure

### Breaking symmetry

When $\Delta c \neq 0$, the EFE difference at $z = 0$ becomes $\Delta G(0) = (1 - 2\bar{a})\Delta c$, which acts as an external field that breaks the $z \to -z$ symmetry.

### Normal-form expansion

Writing the self-consistency equation as $z = \tanh(f(z))$ where $f(z) = -(\gamma/2)\Delta G(z, \Delta c)$, and expanding both sides in a Taylor series around $z = 0$, the leading-order equation (after a standard coordinate shift to absorb the quadratic term) takes the **canonical cusp form**:

```math
0 = \tilde{y} + \tilde{x}\,\tilde{z} + \tilde{z}^3
```

with the following identifications:

| Cusp variable | Active-inference quantity | Role |
|--------------|--------------------------|------|
| $\tilde{z}$ (behavioural) | $\approx z$ (strategy allocation) | Degree of commitment to $\pi_1$ vs $\pi_2$ |
| $\tilde{x}$ (splitting) | Dominated by $(1 - \gamma\mathcal{G})/b$ | Controls whether the system has one or three fixed points |
| $\tilde{y}$ (normal) | Dominated by $f_0/b \propto \Delta c$ | Preference asymmetry breaking the symmetry |

The bifurcation set in the control plane — the region where multistability and sudden jumps can occur — is delimited by the standard cusp curve $4\tilde{x}^3 + 27\tilde{y}^2 = 0$.

### The five catastrophe flags

The cusp structure generates five empirically testable signatures in any population of agents near the critical region:

1. **Bimodality.** The distribution of $z$ across agents is bimodal: some are committed to $\pi_1$, others to $\pi_2$, with few in between.
2. **Inaccessible region.** Intermediate values of $z$ are unstable and therefore rarely observed.
3. **Sudden jump.** As a control parameter traverses the bifurcation set, the agent's commitment switches abruptly from one attractor to the other.
4. **Hysteresis.** If $\Delta c$ is swept forward and then backward, the transition thresholds differ in the two directions because the agent has accumulated asymmetric evidence.
5. **Divergence.** Two agents with slightly different initial conditions can end up in qualitatively different behavioural states.

While in this project we do not go deeper into this aspect of the theory, these catastrophe flags have been studied in more detail and we have provided simulations for them in the separate (yet theoretically related) project on "Phase Transitions in Active Inference" (early learning and developmental psychology context); for instance, for simulations details see Notebook 02 of that project: https://github.com/dtsomoucl/phase-transitions-in-active-inference/blob/main/notebooks/Notebook_02_Simulation_Verifying_Bifurcation.md

---

## 9. Summary of Key Results

The derivation establishes three principal results:

1. **Self-consistency equation.** The learning–action coupling in active inference reduces, under a mean-field closure, to $z = \tanh(-\gamma\Delta G(z)/2)$, which is the Curie–Weiss equation of the Ising model.

2. **Bifurcation condition.** A phase transition in policy selection occurs if and only if $\gamma > \gamma_c = 1/\mathcal{G}_{\max}(p)$, where $\mathcal{G}$ is the coupling function encoding the strength of the learning–action feedback at a given point in the agent's experience.

3. **Cusp-catastrophe normal form.** Near the primary bifurcation, the system has the local geometry of a cusp catastrophe, with policy precision $\times$ feedback strength as the splitting factor and preference asymmetry as the normal factor.

These results are domain-general. The next notebook (Notebook 02) translates them into the specific context of adolescent psychopathology, where $\pi_1$ and $\pi_2$ are reinterpreted as engagement and withdrawal, and the control parameters $\gamma$ and $\Delta c$ acquire clinical meaning.

---

## References

Friston, K., FitzGerald, T., Rigoli, F., Schwartenbeck, P., & Pezzulo, G. (2017). Active inference: a process theory. *Neural Computation*, 29(1), 1–49.

Nishimori, H. (2001). *Statistical physics of spin glasses and information processing: an introduction*. Clarendon Press.

Parr, T., Pezzulo, G., & Friston, K. J. (2022). *Active inference: the free energy principle in mind, brain, and behavior*. MIT Press.

Smith, R., Friston, K. J., & Whyte, C. J. (2022). A step-by-step tutorial on active inference and its application to empirical data. *Journal of Mathematical Psychology*, 107, 102632.

Tsomokos, D. I. (2026). *From Incremental Learning to Developmental Stage Shifts: Phase transitions in active inference.* Code and notebooks: [https://github.com/dtsomoucl/phase-transitions-in-active-inference](https://github.com/dtsomoucl/phase-transitions-in-active-inference).

Van der Maas, H. L., & Molenaar, P. C. (1992). Stagewise cognitive development: an application of catastrophe theory. *Psychological Review*, 99(3), 395–417.
