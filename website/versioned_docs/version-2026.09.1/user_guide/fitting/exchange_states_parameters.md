---
sidebar_position: 6.5
title: Exchange States and Parameters
description: Interpret ChemEx state labels, populations, exchange rates, chemical-shift differences, and state-specific parameters.
---

# Exchange states and parameters

ChemEx uses letters such as A, B, C, and D as stable identifiers for exchange
states. The letters do not rank states by population or energy. In many analyses,
an analyst assigns A to the major or ground state and B to a minor or excited
state, but ChemEx does not enforce that convention or reorder labels after a fit.

This page defines the common parameter conventions. Individual kinetic models
can add constraints or use a different parameterization; see
[Kinetic Models](kinetic_models.md) for the available model families.

## State labels are identifiers

Selecting the `2st` model creates states A and B in that order. A three-state
model similarly uses A, B, and C. These labels remain attached to the same
populations, rates, shifts, relaxation rates, and magnetization components
throughout an analysis, regardless of their fitted populations.

Consequently, a result with $P_B > P_A$ is allowed. It means that the state
labeled B has the larger fitted population; it does not trigger a relabeling.
Assign physical conformations to the labels deliberately and keep the assignment
consistent across input files and experiments.

## Populations in the two-state model

`PB` is the equilibrium population fraction of state B. For `2st`, ChemEx fits
`PB` by default and derives `PA` from normalization:

$$
P_A = 1 - P_B.
$$

Both populations are fractions. The default bounds on `PB` are 0 to 1, but
nothing imposes $P_A > P_B$.

## Two-state exchange rates

`KEX_AB` is the total A↔B exchange rate, in s⁻¹:

$$
KEX_{AB} = K_{AB} + K_{BA}.
$$

The directional names describe the initial and final labels:

- `KAB` means A → B;
- `KBA` means B → A.

For the `2st` parameterization, ChemEx derives them as

$$
K_{AB} = KEX_{AB} P_B
$$

and

$$
K_{BA} = KEX_{AB} P_A.
$$

Thus `KEX_AB` and `PB` are optimizer-controlled by default, whereas `KAB`,
`KBA`, and `PA` are derived. A Method file can change whether eligible
independent parameters are fitted or fixed, but it cannot turn a derived
relationship into an unrelated independent parameter.

Multi-state models preserve the same directional naming rule and stable state
labels. Their allowed pathways and population/rate parameterizations vary by
model, however, so do not extrapolate the two-state formulas to every
multi-state topology.

## Chemical-shift differences

At the user-facing parameter level, chemical shifts and their differences are
in ppm. ChemEx defines

$$
DW_{AB} = CS_B - CS_A.
$$

Where the common shift construction applies, ChemEx therefore derives

$$
CS_B = CS_A + DW_{AB}.
$$

`DW_AB` is generally specific to a spin system or residue: different nuclei can
report different shift changes for the same shared exchange process. Internal
frequency conversions used by a pulse-sequence calculation do not change the
user-facing ppm convention.

## State-specific parameters

Names such as `R1_A`, `R2_A`, `R1_B`, and `R2_B` attach a quantity to a state.
The suffix identifies the state; it does not say how the quantity is obtained.
Depending on the experiment, kinetic model, and Method-file roles, a
state-specific quantity may be fitted, fixed from an input value, or derived
from another parameter or expression.

For example, many common parameterizations initially make state-B relaxation
rates equal to the corresponding state-A rates. Other experiments or Method
roles can require separate state-specific values. Do not infer that every
state-specific parameter is independently fitted merely because it appears in
an output file.

## What changes when A and B are swapped

Relabeling a two-state description is a coordinated transformation, not a
population sort. If old A becomes new B and old B becomes new A, then:

- the new `PB` is the old $P_A = 1 - P_B$;
- the new `KAB` represents the old B → A rate, and the new `KBA` represents the
  old A → B rate;
- the new `CS_A` is the old `CS_B`;
- `DW_AB` reverses sign.

All other state-specific quantities must be exchanged consistently. This simple
two-state relabeling example does not imply that an arbitrary relabeling leaves
every multi-state model or topology unchanged.

## Practical interpretation checklist

Before interpreting fitted state parameters, check that:

- the physical assignment of each letter is documented outside ChemEx;
- all experiment and parameter files use that assignment consistently;
- directional rates are read from the first state to the second;
- `DW_AB` is interpreted as B minus A;
- fitted, fixed, and derived parameter roles are distinguished;
- no conclusion depends only on assuming that A must be the most populated
  state.

For the next part of the fitting convention, see
[Scaling, Uncertainties, and Residuals](scaling_uncertainties_residuals.md).
