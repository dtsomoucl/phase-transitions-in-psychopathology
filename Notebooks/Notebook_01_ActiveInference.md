# Notebook 01 — Mean-Field Transition Structure in a Minimal Active-Inference Model

**Series:** *Phase Transitions in Psychopathology*  
**Status:** Analytical companion to the manuscript  
**Author:** Dimitris I. Tsomokos

## 1. Scope

This notebook records the derivation used in the manuscript. The model is a deliberately minimal two-policy POMDP. It is designed to isolate learning–action feedback, not to reproduce the heterogeneity of symptoms, environments, or treatments in clinical populations.

Three levels must be kept separate:

1. exact algebra inside the specified two-policy generative model;
2. a mean-field closure that replaces stochastic Dirichlet updates with expected evidence counts; and
3. a local normal-form analysis near the symmetric bifurcation.

_Note_ - The main result is a **local normal form**. We make this clear from the outset to avoid confusion, namely, the finite stochastic agent is not globally identical to an Ising model (nor, of course, is it isomorphic to a symptom network).

## 2. Model

There are two hidden states, two observations, and two policies. Policy $\pi_1$ selects state $s_1$ and policy $\pi_2$ selects state $s_2$. The deterministic policy-to-state mapping is an analytical simplification.

The learned likelihood matrix is

```math
\mathbf A = \begin{pmatrix} a & 1-b \\ 1-a & b \end{pmatrix},
```

where $a=\hat P(o_1\mid s_1)$ and $b=\hat P(o_2\mid s_2)$. With symmetric prior concentration $\alpha_0$ per column and true environmental discriminability $p>1/2$,

```math
a = \frac12 + \left(p-\frac12\right)\frac{n_1}{\alpha_0+n_1},
\qquad
b = \frac12 + \left(p-\frac12\right)\frac{n_2}{\alpha_0+n_2}.
```

The counts $n_1$ and $n_2$ are sufficient statistics of the sampled observation history. The base model retains them permanently; later extensions add symmetric or strategy-specific retention.

## 3. Expected-free-energy difference

Let $\Delta c=c_1-c_2$ denote the preference-field difference and let $\mathcal H(x)$ be binary entropy. The exact expected-free-energy difference governing the two-policy softmax is

```math
\boxed{\Delta G = G(\pi_1)-G(\pi_2)
= (1-a-b)\Delta c + \mathcal H(a)-\mathcal H(b).}
```

This is a difference between two policy-specific expected free energies.

The posterior probability of policy $\pi_1$ is

```math
P(\pi_1)=\frac{1}{1+\exp(\gamma\Delta G)},
```

where $\gamma>0$ is policy precision. In this binary architecture, $\Delta G$ determines which policy is favoured and $\gamma$ determines how sharply that advantage is expressed. This decomposition is not a claim that arbitrary multi-policy softmax models have exactly two independent controls.

## 4. Mean-field closure

Let $N=n_1+n_2$, $\phi=n_1/N$, and

```math
z=2\phi-1.
```

With $\tau=N/(2\alpha_0)$,

```math
a(z)=\frac12+\left(p-\frac12\right)
\frac{(1+z)\tau}{1+(1+z)\tau},
```

```math
b(z)=\frac12+\left(p-\frac12\right)
\frac{(1-z)\tau}{1+(1-z)\tau}.
```

Replacing stochastic count allocation by its expectation gives the equilibrium condition

```math
\boxed{z=\tanh\!\left[-\frac{\gamma}{2}\Delta G(z,\tau,\Delta c)\right].}
```

This has the same self-consistency form as a Curie–Weiss mean-field equation. The correspondence concerns expected order-parameter dynamics after closure. The effective coupling is state- and time-dependent (and the finite-agent stochastic POMDP is not thereby transformed into a physical spin system; spin-system language is only used as an analogy).

## 5. Bifurcation and local cusp

In the symmetric case $\Delta c=0$, $z=0$ is a fixed point. Its stability is governed by

```math
\mathcal G(\tau;p)=
\left(p-\frac12\right)\frac{\tau}{(1+\tau)^2}
\left|\log\frac{1-\bar a(\tau)}{\bar a(\tau)}\right|,
```

where

```math
\bar a(\tau)=\frac12+\left(p-\frac12\right)\frac{\tau}{1+\tau}.
```

The balanced fixed point loses stability when $\gamma\mathcal G(\tau;p)=1$. Because $\mathcal G$ rises and later falls as learning accumulates, multistability occurs only over a finite range of rescaled experience for suitable parameter values.

To characterise the nearby geometry, write the fixed-point residual as $F(z;\gamma,\Delta c)=0$ and expand around the symmetric critical point. Retaining terms through third order and applying the standard coordinate shift that removes the quadratic term gives the local canonical form

```math
0=\tilde y+\tilde x\tilde z+\tilde z^3.
```

The splitting control $\tilde x$ is governed primarily by departure from the critical coupling, and the normal control $\tilde y$ is governed primarily by the symmetry-breaking field $\Delta c$. The local fold set is

```math
4\tilde x^3+27\tilde y^2=0.
```

This is a local cusp normal-form equivalence. It is not a global identity between the full active-inference process and a cusp potential.

## 6. Finite-agent interpretation

The mean-field equation locates fixed-point structure. Finite simulated agents traverse that structure stochastically because policies and observations are sampled. Consequently:

- a gradual control schedule can produce abrupt individual changes;
- agents need not cross at the same trial;
- an ensemble mean can change smoothly even when individual traces jump; and
- variance can peak in the transition window without establishing a clinical early-warning signal.

Any claim about a clinical transition requires participant-level longitudinal measurement and model fitting. Population regressions cannot identify the fold geometry or the latent order parameter.

## 7. Reproducibility map

- `Python_code/core_functions.py`: entropy, EFE, posterior, and fixed-point functions.
- `Python_code/sim_psychopathology.py`: onset schedule and finite-agent simulations.
- `Python_code/psychopathology_regimes.py`: shared parameter regimes.
- `Python_code/step6_verified_withdrawal.py`: state-anchored boundary-condition checks.
- `main.py`: complete simulation sweep.

The companion repository, <https://github.com/dtsomoucl/phase-transitions-in-active-inference>, contains the earlier domain-general derivation and additional algebraic detail.

