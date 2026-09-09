---
status: accepted
---

# Make generic N-state topology and domain authoritative

ChemEx will define each Generic N-state Model from one explicit topology record:
its public model name, ordered states, and Structural Exchange Edges. Unsuffixed
`Nst` names are complete graphs, `Nst_linear` names are consecutive chains, and
`Nst_fork` names are A-centered stars. `3st_triangle` remains an exact
compatibility name for complete `3st`. Structural absence is declared directly;
it is not represented by generating a complete graph and deleting or fixing an
edge.

The model owns two scientific functions. The population-complement function
validates the closed Population Simplex and derives `PA` from the independent
non-A populations. The pair-rate function validates a structural `KEX` and its
two endpoint populations, then derives both directional rates together. It
preserves `KEX = k_forward + k_reverse`, detailed balance, exact-zero endpoints,
and exact-zero `KEX` without a denominator floor. A positive `KEX` at two zero
endpoints, a negative input, or a mathematically positive but unrepresentable
directional rate is a scientific-domain failure.

Direct TRF compiles the model-owned Population Simplex into boundary-capable
Feasible Coordinates. Fixed populations and user bounds reduce the remaining
mass available to fitted populations; public continuation and output remain in
`PB`, `PC`, and higher-state coordinates. Deterministic resampling reuses Direct
TRF. MCMC remains in public coordinates and maps model-domain rejection through
its existing invalid-density contract. Deterministic Uncertainty uses registered
analytic function derivatives and retains its existing one-sided-boundary
warning policy for simplex faces.

Every complete-model Structural Exchange Edge is estimation-capable with the
same default. Consequently, bare `4st`, `5st`, and `6st` deliberately change
`KEX_AD` and `KEX_BD` from the historical zero/fixed behavior to `200 s^-1` and
fitted. An explicit parameter value of zero still disables either edge. To
reproduce the old role as well as the old value, a Method Step must explicitly
fix both parameters. Generated restart files preserve values and bounds but not
historical roles, and omitted parameter defaults have no reliable age signal;
ChemEx will not infer either intent.

## Considered options

- Recursive construction from smaller model dictionaries was rejected because
  topology depended on merge history and obscured the special AD/BD defaults.
- A generic graph DSL was rejected because three small edge constructors and an
  explicit public table expose the scientific contract more clearly.
- Complete-graph generation followed by edge deletion or dormant zero/fixed
  edges was rejected because structural absence and dynamic `KEX = 0` have
  different public meanings.
- A softmax population parameterization, clipping, renormalization, penalties,
  and epsilon floors were rejected because they remove exact boundaries or
  change supplied scientific values.
- A generic MCMC or optimizer rewrite was rejected because the existing
  model-owned feasible-coordinate and invalid-density seams already provide the
  required authority.
- Restart-age heuristics and hidden compatibility defaults were rejected because
  old artifacts do not preserve enough role provenance to infer user intent.

## Consequences

Scientists can verify all generic topologies by inspection, and all inputs,
constraints, restart values, deterministic trials, MCMC proposals, and
resampling fits cross the same model-owned Population Simplex and pair-rate
authority before reaching an NMR exchange matrix. Complete models can still
represent disconnected static mixtures by fixing selected structural `KEX`
values to zero; irreducibility is not required.

The public population and exchange names, condition scopes, valid ordinary
three-state behavior, directional orientation, NMR basis, and output coordinate
conventions remain unchanged. New 4/5/6-state linear and fork names are added.
The active AD/BD defaults for bare 4/5/6-state models, strict domain failures,
and newly qualified propagated directional-rate uncertainty are intentional
behavior changes.
