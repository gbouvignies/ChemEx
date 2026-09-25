---
sidebar_position: 7
---

# Kinetic Models

The kinetic model (specified with `-d` or `--model`) defines the exchange
system used for analysis. The following are **all** runtime-registered base
names. A compatibility spelling selects the indicated scientific model, but
ChemEx retains the spelling supplied by the user in model identity and run
provenance.

<!-- kinetic-model-names:start -->
| Model name | Physical system and topology | Canonical scientific equivalent |
| --- | --- | --- |
| `1st` | One non-exchanging state | — |
| `2st` | A ↔ B exchange (default) | — |
| `2st_binding` | Free protein ↔ ligand-bound protein | — |
| `2st_eyring` | A ↔ B with Eyring temperature dependence | — |
| `2st_hd` | Protonated ↔ deuterated state | — |
| `2st_monomer_dimer` | Monomer ↔ dimer | — |
| `2st_monomer_tetramer` | Monomer ↔ tetramer | — |
| `2st_monomer_trimer` | Monomer ↔ trimer | — |
| `2st_rs` | Historical residue-specific two-state exchange | `2st.rs` |
| `3st` | Complete A/B/C exchange graph | — |
| `3st_binding_cs` | Apo conformational selection, then ligand binding | — |
| `3st_binding_if` | Ligand binding, then induced fit | — |
| `3st_binding_partner_2st` | Protein binding two interconverting ligand forms | — |
| `3st_double_binding` | Two competing bound complexes | — |
| `3st_eyring` | Historical linear three-state Eyring exchange | `3st_eyring_linear` |
| `3st_eyring_fork` | A-centered Eyring fork | — |
| `3st_eyring_linear` | A ↔ B ↔ C Eyring exchange | — |
| `3st_fork` | A-centered three-state fork | — |
| `3st_linear` | A ↔ B ↔ C exchange | — |
| `3st_monomer_dimer_tetramer` | Sequential monomer/dimer/tetramer association | — |
| `3st_monomer_dimer_trimer` | Sequential monomer/dimer/trimer association | — |
| `3st_triangle` | Historical complete three-state graph | `3st` |
| `4st` | Complete four-state exchange graph | — |
| `4st_binding_3_bound_states` | Free protein and three interconverting bound forms | — |
| `4st_binding_partner_2st` | Protein binding two ligand forms, with a third bound form | — |
| `4st_eyring` | Complete four-state Eyring exchange graph | — |
| `4st_fork` | A-centered four-state fork | — |
| `4st_hd` | Two conformers, each with protonated/deuterated forms | — |
| `4st_linear` | Four-state chain | — |
| `5st` | Complete five-state exchange graph | — |
| `5st_fork` | A-centered five-state fork | — |
| `5st_linear` | Five-state chain | — |
| `6st` | Complete six-state exchange graph | — |
| `6st_fork` | A-centered six-state fork | — |
| `6st_linear` | Six-state chain | — |
<!-- kinetic-model-names:end -->

In these models, each state in the exchange process is represented with a unique
parameter suffix (`A`, `B`, `C`, `D`, etc.). For example, `R1_A` and `R2_B`
refer to the R<sub>1</sub> rate of state A and the R<sub>2</sub> rate of state B.
In many analyses A is assigned to the major or ground state and B to a minor or
excited state, but this is an analyst-defined convention rather than an ordering
enforced by ChemEx. See [Exchange States and Parameters](exchange_states_parameters.md)
for populations, directional rates, shift signs, and label swapping.

## Composable modifiers

The suffixes `.rs`, `.mf`, and `.tc` may be combined in any order; each changes
only its own parameter family. Accepted redundant combinations remain accepted.
For example, `2st_rs.rs` has the same kinetic scope as `2st.rs`, while the
selected spelling remains in provenance.

**`.rs`** adds residue/group identity to eligible kinetic quantities, including
their derived populations and rates. It leaves global sample conditions such as
`D2O` global and changes neither equations, domains, default roles, nor base
model physics. The compatibility base name `2st_rs` corresponds scientifically
to `2st.rs`.

**`.mf`** uses independent model-free `TAUC` (ns), `S2` (unitless, 0–1), and
`KHH` (s⁻¹, where applicable) in place of independently assigned relaxation
rates. Their defaults are respectively 4 ns, 0.9, and 0 s⁻¹; their bounds are
[0, 1000] ns, [0, 1], and [0, 1e6] s⁻¹. Relaxation and cross-correlation rates
are model-derived, with magnetic-field and nucleus dependence retained on the
rates. `KHH` stays with the physical proton even for reversed `HN`/`NH`
orientation; deuterated labels use the corresponding deuterated rate model.
Exact zero `TAUC`, `S2`, or `KHH` is retained, and the derived rates must still
satisfy ChemEx's relaxation feasibility conditions. For an example, see
`CEST_15N_TR/` under `Examples/Experiments/`.

## Linear chemical-shift temperature coefficients

For any kinetic model, add `.tc` to use reference-centered linear temperature
polynomials for both the state-A chemical shift and every A-to-X chemical-shift
difference. The coefficients are `CS0_A`, `CS1_A`, `DW0_AX`, and `DW1_AX`, and
the analysis-wide reference temperature is the protected `[GLOBAL]` constant
`TREF` (25 °C by default). Every `.tc` Experiment must provide `temperature` in
degrees Celsius.

See [Temperature-dependent chemical shifts](temperature_dependent_shifts.md) for
the equations, parameter roles, multistate and Eyring composition, fitting
limitations, and output behavior.

## Reading kinetic parameter settings

In the [runtime parameter inventory](#runtime-parameter-inventory), **input**
means an independent coordinate or fixed value that can be supplied through a
parameter file; **derived** means ChemEx calculates it, and **report-only**
means it is calculated for output but does not block dynamics when its value
cannot be represented. A default fitted role is marked `fit`; `fixed` values
can still be selected by an applicable Method Plan unless protected by their
model. `estimate` marks an explicit estimation capability. Bounds shown there
are the registered defaults; parameter files can override ordinary bounds.

Scopes use `T` for temperature, `P` for total protein, `L` for total ligand,
and `D` for D2O fraction. `group` identifies residue/group-specific kinetics;
`global` means no spin or group qualifier. A missing condition dimension means
the parameter is shared across that dimension. Independent binding and
oligomerization inputs are normally temperature-scoped and shared across spin
systems and concentrations; their derived populations and rates also depend on
total concentrations. A model may have further domain checks beyond the
registered default bounds; the family sections state the important ones.

### One and two states

`1st` has only state A, with `PA = 1` and no kinetic input or exchange edge.
For ordinary `2st`, A ↔ B has independent population `PB` (fraction, default
0.05, [0, 1]) and total exchange scale `KEX_AB` (s⁻¹, default 200,
[0, 1e6]). `PA = 1 − PB`; `KAB` and `KBA` are derived directional rates.
Both zero population endpoints and `KEX_AB = 0` are valid, with exact zero
rates where appropriate. These quantities are scoped by temperature and total
protein/ligand concentrations and are global across spins unless `.rs` is
selected. `2st_rs` has this same physics and residue/group scope.

### H/D solvent exchange

`2st_hd` uses A (protonated) ↔ B (deuterated). The independent physical inputs
are global solvent `D2O` (fraction), residue/group-specific `KDH` (s⁻¹), and
residue/group-specific fractionation factor `PHI` (unitless). `KDH` is fitted by
default; `PHI` is fixed by default. ChemEx derives both populations and
directional rates. `KDH = 0` freezes both directions without erasing the
solvent-defined equilibrium composition. `D2O` remains global under `.rs`.

`4st_hd` combines conformers A/B with their deuterated forms C/D. The square
has conformational A ↔ B and C ↔ D edges, and solvent-exchange A ↔ C and
B ↔ D edges; there is no A ↔ D or B ↔ C edge. `POP_B` selects the conformer
fraction, `KEX_AB` the conformational scale, and group-specific `KDH_A`,
`KDH_B`, and `PHI_A` the solvent kinetics. `PHI_B` is derived from `PHI_A`.
`D2O` is global. All four populations and eight directional rates are derived.
The conformational scale defaults to exact zero, leaving the two H/D pairs
without interconversion; zero solvent-exchange scales similarly freeze their
corresponding H/D edges without changing equilibrium composition.

### Ligand binding

`2st_binding` describes free protein A + ligand L ↔ complex B. Supply
`KD` (M) and `KOFF` (s⁻¹), together with positive total protein concentration
and nonnegative total ligand concentration. ChemEx derives free concentrations,
populations, and tagged directional rates. `KON` is report-only. `KD` must be
strictly positive even if bounds are overridden; zero ligand gives only free
protein, while `KOFF = 0` freezes tagged exchange without changing the
equilibrium species.

For all binding models below, `KD`/`KD_APP` are in M, `KOFF` and `KEX` in
s⁻¹, and `KEQ` ratios are dimensionless. Independent binding inputs are fitted
by default. Exact-zero `KOFF`/`KEX` freezes the corresponding tagged edge
without changing composition. Positive protein total and positive dissociation
constants are required; zero ligand is allowed. The derived concentrations,
populations, and directional rates follow the equilibrium species, while
`KON` quantities are report-only. The [runtime inventory](#runtime-parameter-inventory)
gives each model's exact public parameter names, defaults, bounds, and scope.

| Model | States and reaction scheme | Independent physical inputs | Important derived quantities and boundaries |
| --- | --- | --- | --- |
| `3st_double_binding` | A = free protein; A + L ↔ B and A + L ↔ C are competing bound complexes | `KD_AB`, `KOFF_AB`, `KD_AC`, `KOFF_AC` | Free species, `PA`–`PC`, tagged edge rates; no B ↔ C edge |
| `3st_binding_partner_2st` | Ligand L1 ↔ L2; A + L1 ↔ B and A + L2 ↔ C, with B ↔ C tagged exchange | `KD_AB`, `KOFF_AB`, `KD_AC`, `KOFF_AC`, `KEQ`, `KEX_BC` | `KEQ = 0` removes L2 and C exactly; zero ligand leaves A only, although the B ↔ C rate split remains defined by the input equilibrium ratios |
| `4st_binding_partner_2st` | Ligand L1 ↔ L2; A + L1 ↔ B and A + L2 ↔ C, then C ↔ D is a third bound form; B ↔ C also exchanges | `KD_AB`, `KOFF_AB`, `KD_AC`, `KOFF_AC`, `KEQ_L`, `KEQ_PL`, `KEX_BC`, `KEX_CD` | `KEQ_L = 0` removes L2, C, D; `KEQ_PL = 0` removes D exactly |
| `4st_binding_3_bound_states` | A = free protein; A + L ↔ B ↔ C ↔ D, with three bound forms | `KD_APP`, `KOFF_AB`, `KEQ_BC`, `KEX_BC`, `KEQ_CD`, `KEX_CD` | `KD_AB`, `KD_EFF`, free/bound concentrations and `PA`–`PD` are derived; zero `KEQ_BC` removes C/D and zero `KEQ_CD` removes D exactly |

In the partner models, `KEQ` or `KEQ_L` sets the free L2/L1 ratio;
`KEQ_PL` sets the D/C bound-state ratio. In the three-bound-state model,
`KEQ_BC` and `KEQ_CD` set C/B and D/C, and `KD_APP` refers to the total
bound protein concentration. Ratios of zero remove the named downstream
species exactly; a positive dissociation constant is still required.

### Oligomerization

For direct `2st_monomer_dimer`, `2st_monomer_trimer`, and
`2st_monomer_tetramer`, A is monomer and B is the named oligomer: respectively
2A ↔ B, 3A ↔ B, and 4A ↔ B. Independent `KD` (M to the power of oligomer size
minus one) and `KOFF` (s⁻¹) determine chemical equilibrium and tagged exchange.
ChemEx derives monomer/oligomer concentrations, protein-unit populations, and
directional rates. `KD` must be positive; `KOFF = 0` freezes tagged exchange
while preserving chemical equilibrium. Zero total protein leaves only state A
as the reference population and zero association.

The sequential models add C: `3st_monomer_dimer_trimer` has
2A ↔ B and A + B ↔ C, while `3st_monomer_dimer_tetramer` has
2A ↔ B and 2B ↔ C. Supply `KD1`, `KD2`, `KOFF1`, and `KOFF2`; the `KD`
units follow each reaction's stoichiometry. Monomer/dimer/higher-oligomer
concentrations, all three populations, and tagged rates are derived. The
trimer model includes A ↔ C, B ↔ C, and A ↔ B tagged edges; the tetramer
model includes only A ↔ B and B ↔ C. Zero `KOFF1` or `KOFF2` freezes the
associated tagged directions without erasing equilibrium species. Zero `KD1`
or `KD2` is invalid.

## Generic N-state models

The generic family supports three through six states. Its public model name is
the topology contract:

| Model | Topology |
| --- | --- |
| `3st`, `4st`, `5st`, `6st` | Complete graph: every unordered state pair exchanges |
| `3st_linear`, `4st_linear`, `5st_linear`, `6st_linear` | Chain A ↔ B ↔ C ↔ D ↔ E ↔ F, truncated at the selected state count |
| `3st_fork`, `4st_fork`, `5st_fork`, `6st_fork` | A-centered star: A exchanges directly with every other state |
| `3st_triangle` | Exact historical compatibility name for the canonical complete `3st` model |

`PB`, `PC`, and any higher non-A populations are independent fractions; `PA`
is derived so that all populations sum to one. ChemEx enforces the complete
closed simplex: every population must be finite and nonnegative and their sum
must be exactly representable as no greater than one. Exact zero populations
and `PA = 0` are valid. Values are not clipped, renormalized, or replaced by a
small positive floor.

For every structural edge `i` ↔ `j`, `KEX_ij` is the total exchange scale in
s⁻¹:

$$
KEX_{ij} = K_{ij} + K_{ji}.
$$

ChemEx derives the directional rates from the prescribed populations:

$$
K_{ij} = KEX_{ij}\frac{P_j}{P_i + P_j}, \qquad
K_{ji} = KEX_{ij}\frac{P_i}{P_i + P_j}.
$$

The first state in a directional name is its source, so `KAB` means A → B.
`KEX_ij = 0` produces exact zero in both directions and dynamically switches
off an existing structural edge. If exactly one endpoint population is zero,
the outgoing rate from that endpoint equals `KEX_ij` and the reciprocal rate is
zero. If both endpoint populations are zero, `KEX_ij = 0` remains valid, but a
positive `KEX_ij` is rejected because its directional split is undefined.
Positive rates too small to be represented as positive binary64 values are
also rejected instead of being changed into structural zero.

An absent edge in a `_linear` or `_fork` model is different from an edge fixed
to zero in a complete model: the absent edge has no `KEX` or directional-rate
parameters and cannot be selected by a Method Step. In a complete model,
setting selected structural `KEX` values to zero can create disconnected
kinetic components. Positive prescribed populations may remain in multiple
components; ChemEx treats them as static mixture weights for non-interconverting
subensembles and does not require an irreducible generator or a unique global
stationary distribution.

Direct TRF and its Monte Carlo/bootstrap reruns use closed, bounded simplex
coordinates internally while continuation and output retain public `PB`, `PC`,
and higher-state coordinates. MCMC remains in those public coordinates and
rejects proposals outside the simplex through its normal zero-density path.
Qualified interior deterministic covariance is propagated analytically to
`PA` and the directional rates. At a population boundary, symmetric errors are
reported only with ChemEx's boundary warning; the positive-`KEX` zero/zero split
has no valid central rate or derivative.

### Breaking default change for complete 4/5/6-state models

Bare `4st`, `5st`, and `6st` are now genuinely complete by default. Every
structural `KEX`, including `KEX_AD` and `KEX_BD`, starts at `200 s⁻¹`, is
bounded to `[0, 1e6] s⁻¹`, and is fitted by default. Earlier releases inherited
historical `KEX_AD = KEX_BD = 0` fixed defaults. This is an intentional breaking
change.

To reproduce the old sparse behavior for any of `4st`, `5st`, or `6st`, use a
normal parameter file:

```toml title="legacy-complete-parameters.toml"
[GLOBAL]
KEX_AD = 0.0
KEX_BD = 0.0
```

and explicitly preserve the roles in a version 2 Method Plan:

```toml title="legacy-complete-method.toml"
FORMAT_VERSION = 2

[FIT]
ROLES = [{ FIX = ["KEX_AD", "KEX_BD"] }]
```

Select the required state count normally:

```text
chemex fit ... -d 4st -p legacy-complete-parameters.toml -m legacy-complete-method.toml
chemex fit ... -d 5st -p legacy-complete-parameters.toml -m legacy-complete-method.toml
chemex fit ... -d 6st -p legacy-complete-parameters.toml -m legacy-complete-method.toml
```

Old generated `run_info/restart.toml` files preserve the numerical zero values
and their bounds, but do not preserve the historical fitted/fixed roles. When
such a restart is loaded under the new complete model, `KEX_AD` and `KEX_BD`
therefore become ordinary fitted coordinates unless the Method Step above fixes
them. Plain historical parameter files that omitted the two values carry no
reliable version or intent signal; ChemEx applies the new complete defaults
rather than guessing their age.

## Three-State Association Models

The `3st_binding_cs` and `3st_binding_if` models separate equilibrium
composition from tagged-magnetization exchange speed. In both models,
`KD_APP` and the equilibrium ratio determine `PA`, `PB`, and `PC`; changing a
`KEX` or `KOFF` value does not change those populations.

All independent association parameters are scoped by temperature. `KD_APP` is
in M, `KOFF_*` and `KEX_*` are in s⁻¹, and `KEQ_*` is dimensionless. The upper
bounds shown below are broad safety defaults rather than scientific priors;
parameter files can override them using the normal explicit-bound syntax.

### Conformational selection: `3st_binding_cs`

The reaction scheme is

```text
A ⇌ B
    B + L ⇌ C
```

where A is apo conformer 1, B is the apo binding-competent conformer, and C is
the bound complex. The independent parameters are:

| Parameter | Meaning | Default | Default bounds |
| --- | --- | ---: | ---: |
| `KD_APP` | apparent dissociation constant against total unbound protein | `1e-6` M | `(0, 1.0]` M |
| `KOFF_BC` | tagged B←C dissociation-rate scale | `100.0` s⁻¹ | `[0, 1e6]` s⁻¹ |
| `KEQ_AB` | apo equilibrium ratio `B / A` | `1.0` | `(0, 1e6]` |
| `KEX_AB` | total tagged A↔B exchange scale | `200.0` s⁻¹ | `[0, 1e6]` s⁻¹ |

Writing `q = KEQ_AB` and `U = A + B`, the apo distribution is
`A = U / (1 + q)` and `B = U q / (1 + q)`. The apparent equilibrium definition
is `KD_APP = U L / C`, so the intrinsic binding constant is
`KD_BC = KD_APP q / (1 + q)`. Directional conformational rates are
`KAB = KEX_AB q / (1 + q)` and `KBA = KEX_AB / (1 + q)`. The tagged binding
edge uses `KCB = KOFF_BC` and detailed balance, `PB KBC = PC KCB`.

`KEX_AB = 0` sets both A↔B tagged rates to zero without changing the apo
equilibrium. `KOFF_BC = 0` freezes both tagged binding directions without
changing composition. `L_TOTAL = 0` retains the A/B apo distribution and makes
the tagged association rate zero. `KD_APP = 0`, `KEQ_AB = 0`, and
`P_TOTAL = 0` are rejected; exact-zero `KEQ_AB` is incompatible with this
finite reversible apparent-binding parameterization.

`KAB`, `KBA`, `KD_BC`, `KON_BC`, `C_L`, `KBC`, `KCB`, `PA`, `PB`, and `PC`
remain public derived outputs. `KD_BC` and `KON_BC` are report-only, so an
unrepresentable intrinsic value cannot block otherwise finite tagged dynamics.
When representable, these report-only values are included in normal parameter
output with propagated deterministic uncertainty when its registered derivative
qualifies; otherwise the output gives the explicit uncertainty-unavailable reason.

### Induced fit: `3st_binding_if`

The reaction scheme is

```text
A + L ⇌ B ⇌ C
```

where A is free protein, B is the first bound complex, and C is the rearranged
bound complex. The independent parameters are:

| Parameter | Meaning | Default | Default bounds |
| --- | --- | ---: | ---: |
| `KD_APP` | apparent dissociation constant against total bound protein | `1e-3` M | `(0, 1.0]` M |
| `KOFF_AB` | tagged A←B dissociation-rate scale | `100.0` s⁻¹ | `[0, 1e6]` s⁻¹ |
| `KEQ_BC` | bound-state equilibrium ratio `C / B` | `1.0` | `[0, 1e6]` |
| `KEX_BC` | total tagged B↔C exchange scale | `200.0` s⁻¹ | `[0, 1e6]` s⁻¹ |

Writing `q = KEQ_BC`, the apparent equilibrium definition is
`KD_APP = A L / (B + C)`, with `B:C = 1:q`; therefore
`KD_AB = KD_APP (1 + q)`. Directional conformational rates are
`KBC = KEX_BC q / (1 + q)` and `KCB = KEX_BC / (1 + q)`. The tagged binding
edge uses `KBA = KOFF_AB` and detailed balance, `PA KAB = PB KBA`.

`KEX_BC = 0` and `KOFF_AB = 0` freeze their respective tagged edges without
changing composition. `KEQ_BC = 0` is supported exactly: C has zero equilibrium
weight, `KBC = 0`, `KCB = KEX_BC`, and `KD_AB = KD_APP`. `L_TOTAL = 0` puts all
protein in free A and makes the tagged association rate zero. `KD_APP = 0` and
`P_TOTAL = 0` are rejected.

`KBC`, `KCB`, `KD_AB`, `KON_AB`, `C_L`, `KAB`, `KBA`, `PA`, `PB`, and `PC`
remain public derived outputs. `KD_AB` and `KON_AB` are report-only, with the same
representability and uncertainty-reporting policy as `KD_BC` and `KON_BC` above.

### Migrating directional-rate inputs

Directional rates remain available as outputs but can no longer be supplied as
independent parameter values, Method Plan roles, constraints, or grid targets.
For `3st_binding_cs`, replace a positive legacy pair using
`KEQ_AB = KAB / KBA` and `KEX_AB = KAB + KBA`. For `3st_binding_if`, use
`KEQ_BC = KBC / KCB` and `KEX_BC = KBC + KCB`.

A legacy `(0, 0)` pair determines only `KEX = 0`; its equilibrium ratio is
ambiguous and must be chosen explicitly. For induced fit, `KBC = 0` with
`KCB > 0` maps exactly to `KEQ_BC = 0` and `KEX_BC = KCB`. A positive forward
rate with a zero reverse rate would require an infinite equilibrium ratio and
cannot be migrated to the finite parameter domain. Conformational selection
also rejects a legacy `KAB = 0` mapping because `KEQ_AB` must be positive.

The migration diagnostic prints a numerical replacement only when the exact
positive ratio and sum are representable in binary64. If either result is outside
that domain, it asks for manual model/parameter reconsideration instead of
printing a false zero or infinity. A representable replacement can still exceed
the new broad default upper bound of `1e6`. In that case keep the exact migrated
value and override the bound explicitly with the normal three-value parameter
syntax, for example `KEX_BC = [1200000.0, 0.0, 1200000.0]`; for positive CS
`KEQ_AB`, use a positive lower bound such as `5e-324`.

## Temperature-Dependent Eyring Models

The `2st_eyring`, three-state Eyring variants, and `4st_eyring` models
implement exchange systems with temperature-dependent rate constants calculated
using Eyring transition state theory. These models are particularly useful for
studying exchange processes where thermodynamic parameters govern the
temperature dependence of exchange rates.

### Thermodynamic convention and temperature

ChemEx input temperatures are in degrees Celsius. Eyring calculations convert
them to Kelvin internally with `T = temperature + 273.15`. The public domain is
every finite temperature strictly above -273.15 °C; absolute zero, lower
temperatures, NaN, and infinities are rejected. There is no artificial upper
temperature limit.

State A is the thermodynamic reference, so `H_A = S_A = 0`. Every `DH_I` and
`DS_I` parameter is the coordinate of state I relative to A. Likewise,
`DH_IJ` and `DS_IJ` are the shared IJ transition-state coordinates relative to
A, not separate forward activation parameters. Enthalpies use J mol⁻¹ and
entropies use J mol⁻¹ K⁻¹.

For the directional transition i → j, ChemEx therefore uses

```
Delta H double dagger (i -> j) = H_TS_ij - H_i
Delta S double dagger (i -> j) = S_TS_ij - S_i

log(k_ij) = log(k_B / h) + log(T)
            + (S_TS_ij - S_i) / R
            - (H_TS_ij - H_i) / (R * T)
```

Here `T` is in Kelvin, `k_ij` is in s⁻¹, and the transmission coefficient is
one. A shared transition state makes the forward/reverse pair obey

```
log(k_ij / k_ji) = -((H_j - H_i) - T * (S_j - S_i)) / (R * T)
```

The equilibrium populations are calculated directly from state coordinates,
not from rate magnitudes:

```
log(w_i) = S_i / R - H_i / (R * T)
p_i = w_i / sum(w)
```

The topology still selects the kinetic edges. This population authority ensures
detailed balance for every present edge and thermodynamic cycle consistency in
`4st_eyring`, even when all transition states are shifted to make kinetics much
slower.

Rates are evaluated in the log domain and are not clipped or saturated. If a
mathematically positive rate lies below the minimum positive binary64 value or
above the maximum finite binary64 value, parameter evaluation fails explicitly
instead of producing zero, infinity, or a capped rate.

Finite Eyring state coordinates likewise imply strictly positive populations.
If a normalized state population lies below the minimum positive binary64 value,
ChemEx rejects the evaluation instead of converting that state to structural
zero.

### 2st_eyring Model Parameters

The `2st_eyring` model uses the following thermodynamic parameters:

**State Energies (relative to state A):**

- `DH_B`: Enthalpy difference (J/mol) for state B relative to A
- `DS_B`: Entropy difference (J/mol/K) for state B relative to A

**Transition-state coordinates (relative to state A):**

- `DH_AB`: Enthalpy coordinate (J/mol) of the shared AB transition state
- `DS_AB`: Entropy coordinate (J/mol/K) of the shared AB transition state

The model automatically calculates both forward (k_AB) and reverse (k_BA) rate constants from these parameters.

### Three-State Eyring Model Parameters

`3st_eyring_linear` has the linear topology A ↔ B ↔ C.
`3st_eyring_fork` has the fork topology B ↔ A ↔ C. The historical
`3st_eyring` name remains supported as a compatibility name for
`3st_eyring_linear`; it does not select a triangular topology.

Both topologies use the following state parameters:

**State Energies (relative to state A):**

- `DH_B`, `DH_C`: Enthalpy differences (J/mol) for states B, C
- `DS_B`, `DS_C`: Entropy differences (J/mol/K) for states B, C

**Linear transition-state coordinates (`3st_eyring`, `3st_eyring_linear`):**

- `DH_AB`, `DH_BC`: Enthalpy coordinates (J/mol) of the AB and BC transition states
- `DS_AB`, `DS_BC`: Entropy coordinates (J/mol/K) of the AB and BC transition states

The linear models calculate k_AB, k_BA, k_BC, and k_CB. There is no direct A ↔ C
pathway.

**Fork transition-state coordinates (`3st_eyring_fork`):**

- `DH_AB`, `DH_AC`: Enthalpy coordinates (J/mol) of the AB and AC transition states
- `DS_AB`, `DS_AC`: Entropy coordinates (J/mol/K) of the AB and AC transition states

The fork model calculates k_AB, k_BA, k_AC, and k_CA. There is no direct B ↔ C
pathway.

Parameter files used with `3st_eyring` or `3st_eyring_linear` should not contain
`DH_AC` or `DS_AC`; those parameters are not part of the linear topology.

### 4st_eyring Model Parameters

The `4st_eyring` model implements a full 4-state system:

**State Energies (relative to state A):**

- `DH_B`, `DH_C`, `DH_D`: Enthalpy differences (J/mol) for states B, C, D
- `DS_B`, `DS_C`, `DS_D`: Entropy differences (J/mol/K) for states B, C, D

**Transition-state coordinates (relative to state A):**

- `DH_AB`, `DH_AC`, `DH_AD`: Enthalpy coordinates (J/mol) of the AB, AC, and AD transition states
- `DH_BC`, `DH_BD`, `DH_CD`: Enthalpy coordinates (J/mol) of the BC, BD, and CD transition states
- `DS_AB`, `DS_AC`, `DS_AD`: Entropy coordinates (J/mol/K) of the AB, AC, and AD transition states
- `DS_BC`, `DS_BD`, `DS_CD`: Entropy coordinates (J/mol/K) of the BC, BD, and CD transition states

The model automatically calculates all 12 rate constants (k_AB, k_BA, k_AC, k_CA, k_AD, k_DA, k_BC, k_CB, k_BD, k_DB, k_CD, k_DC).

All Eyring state and transition enthalpy parameters have default bounds of
[-2×10⁵, 2×10⁵] J/mol. The corresponding entropy parameters have default
bounds of [-5×10², 5×10²] J/mol/K. Explicit bounds in parameter files override
these defaults.

### Identifiability, uncertainty, and parameter scope

At one temperature, varying both H and S for the same state or transition state
is non-identifiable: the data determine the corresponding free-energy
combination, not H and S separately. Different H/S pairs can therefore give the
same rate or population, and deterministic uncertainty can correctly be
unavailable because the covariance is rank deficient. Measurements at multiple
temperatures can restore mathematical rank, but H/S correlation remains strong
over a narrow temperature interval. Use a meaningful temperature span when H
and S are both fitted.

When the covariance is qualified, ChemEx propagates deterministic uncertainty
to derived directional rates and to the coupled, normalized populations. Eyring
models use the generic MCMC and resampling workflows; each proposal or replicate
recalculates rates and populations through the same thermodynamic authority.

The current ChemEx scoping contract shares Eyring H/S parameters across
temperature, magnetic field, and spin system, but splits them when `p_total` or
`l_total` changes. This concentration scoping is a compatibility contract, not
a universal thermodynamic requirement. Derived rates and populations are scoped
by temperature, `p_total`, and `l_total`; changing magnetic field alone does not
create another chemical kinetic rate. Public parameter names, units, defaults,
and TOML syntax are unchanged.

## Runtime parameter inventory

This compact inventory is projected from the registered settings under
representative conditions (25 °C, 1 mM total protein, 2 mM total ligand,
20% D2O). It lists only kinetic quantities, not experiment-specific shift or
relaxation parameters. Family sections above explain scientific meaning and
units. The exact table is checked against runtime settings and representative
sealed parameter construction by the test suite.

<details>
<summary>Show all kinetic inputs, defaults, bounds, scopes, and derived names</summary>

<!-- kinetic-parameters:start -->
| Model | Independent input: default [bounds] (role; scope) | Derived | Report-only |
| --- | --- | --- | --- |
| `1st` | — | `PA` (all; global) | — |
| `2st` | `KEX_AB` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `PB` 0.05 [0.0, 1.0] (fit; T,P,L; global) | `PA`, `KAB`, `KBA` (T,P,L; global) | — |
| `2st_binding` | `KOFF` 100.0 [0.0, 1000000.0] (fit; T; global); `KD` 0.001 [5e-324, 1.0] (fit; T; global) | `P_FREE`, `L_FREE`, `PL`, `KAB`, `KBA`, `PA`, `PB` (T,P,L; global) | `KON` (T; global) |
| `2st_eyring` | `DH_B` 8000.0 [-200000.0, 200000.0] (fit; P,L; global); `DS_B` 0.0 [-500.0, 500.0] (fit; P,L; global); `DH_AB` 65000.0 [-200000.0, 200000.0] (fit; P,L; global); `DS_AB` 0.0 [-500.0, 500.0] (fit; P,L; global) | `KAB`, `KBA`, `PA`, `PB` (T,P,L; global) | — |
| `2st_hd` | `D2O` 0.2 [0.0, 1.0] (fixed, estimate; D; global); `KDH` 1.0 [0.0, 1000000.0] (fit; T; group); `PHI` 1.1 [0.75, 1.5] (fixed; T; group) | `KAB`, `KBA`, `PA`, `PB` (T,D; group) | — |
| `2st_monomer_dimer` | `KOFF` 100.0 [0.0, 1000000.0] (fit; T; global); `KD` 1e-06 [5e-324, 1.0] (fit; T; global) | `C_MONOMER`, `C_DIMER`, `KAB`, `KBA`, `PA`, `PB` (T,P; global) | — |
| `2st_monomer_tetramer` | `KOFF` 100.0 [0.0, 1000000.0] (fit; T; global); `KD` 1e-06 [5e-324, 1.0] (fit; T; global) | `C_MONOMER`, `C_TETRAMER`, `KAB`, `KBA`, `PA`, `PB` (T,P; global) | — |
| `2st_monomer_trimer` | `KOFF` 100.0 [0.0, 1000000.0] (fit; T; global); `KD` 1e-06 [5e-324, 1.0] (fit; T; global) | `C_MONOMER`, `C_TRIMER`, `KAB`, `KBA`, `PA`, `PB` (T,P; global) | — |
| `2st_rs` | `KEX_AB` 200.0 [0.0, 1000000.0] (fit; T,P,L; group); `PB` 0.05 [0.0, 1.0] (fit; T,P,L; group) | `PA`, `KAB`, `KBA` (T,P,L; group) | — |
| `3st` | `PB` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PC` 0.02 [0.0, 1.0] (fit; T,P,L; global); `KEX_AB` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_AC` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_BC` 200.0 [0.0, 1000000.0] (fit; T,P,L; global) | `PA`, `KAB`, `KBA`, `KAC`, `KCA`, `KBC`, `KCB` (T,P,L; global) | — |
| `3st_binding_cs` | `KD_APP` 1e-06 [5e-324, 1.0] (fit; T; global); `KOFF_BC` 100.0 [0.0, 1000000.0] (fit; T; global); `KEQ_AB` 1.0 [5e-324, 1000000.0] (fit; T; global); `KEX_AB` 200.0 [0.0, 1000000.0] (fit; T; global) | `KAB`, `KBA` (T; global); `C_L`, `KBC`, `KCB`, `PA`, `PB`, `PC` (T,P,L; global) | `KD_BC`, `KON_BC` (T; global) |
| `3st_binding_if` | `KD_APP` 0.001 [5e-324, 1.0] (fit; T; global); `KOFF_AB` 100.0 [0.0, 1000000.0] (fit; T; global); `KEQ_BC` 1.0 [0.0, 1000000.0] (fit; T; global); `KEX_BC` 200.0 [0.0, 1000000.0] (fit; T; global) | `KBC`, `KCB` (T; global); `C_L`, `KAB`, `KBA`, `PA`, `PB`, `PC` (T,P,L; global) | `KD_AB`, `KON_AB` (T; global) |
| `3st_binding_partner_2st` | `KOFF_AB` 100.0 [0.0, 1000000.0] (fit; T; global); `KD_AB` 0.001 [5e-324, 1.0] (fit; T; global); `KOFF_AC` 100.0 [0.0, 1000000.0] (fit; T; global); `KD_AC` 0.001 [5e-324, 1.0] (fit; T; global); `KEQ` 1.0 [0.0, 1000000.0] (fit; T; global); `KEX_BC` 1000.0 [0.0, 1000000.0] (fit; T; global) | `L1_FREE`, `L2_FREE`, `PL1`, `PL2`, `KAB`, `KBA`, `KAC`, `KCA`, `KBC`, `KCB`, `PA`, `PB`, `PC` (T,P,L; global) | `KON_AB`, `KON_AC` (T; global) |
| `3st_double_binding` | `KOFF_AB` 100.0 [0.0, 1000000.0] (fit; T; global); `KD_AB` 0.001 [5e-324, 1.0] (fit; T; global); `KOFF_AC` 100.0 [0.0, 1000000.0] (fit; T; global); `KD_AC` 0.001 [5e-324, 1.0] (fit; T; global) | `PFREE`, `L_FREE`, `KAB`, `KBA`, `KAC`, `KCA`, `PA`, `PB`, `PC` (T,P,L; global) | `KON_AB`, `KON_AC` (T; global) |
| `3st_eyring` | `DH_B` 8000.0 [-200000.0, 200000.0] (fit; P,L; global); `DH_C` 8000.0 [-200000.0, 200000.0] (fit; P,L; global); `DS_B` 0.0 [-500.0, 500.0] (fixed; P,L; global); `DS_C` 0.0 [-500.0, 500.0] (fixed; P,L; global); `DH_AB` 65000.0 [-200000.0, 200000.0] (fit; P,L; global); `DS_AB` 0.0 [-500.0, 500.0] (fixed; P,L; global); `DH_BC` 65000.0 [-200000.0, 200000.0] (fit; P,L; global); `DS_BC` 0.0 [-500.0, 500.0] (fixed; P,L; global) | `KAB`, `KBA`, `KBC`, `KCB`, `PA`, `PB`, `PC` (T,P,L; global) | — |
| `3st_eyring_fork` | `DH_B` 8000.0 [-200000.0, 200000.0] (fit; P,L; global); `DH_C` 8000.0 [-200000.0, 200000.0] (fit; P,L; global); `DS_B` 0.0 [-500.0, 500.0] (fixed; P,L; global); `DS_C` 0.0 [-500.0, 500.0] (fixed; P,L; global); `DH_AB` 65000.0 [-200000.0, 200000.0] (fit; P,L; global); `DS_AB` 0.0 [-500.0, 500.0] (fixed; P,L; global); `DH_AC` 65000.0 [-200000.0, 200000.0] (fit; P,L; global); `DS_AC` 0.0 [-500.0, 500.0] (fixed; P,L; global) | `KAB`, `KBA`, `KAC`, `KCA`, `PA`, `PB`, `PC` (T,P,L; global) | — |
| `3st_eyring_linear` | `DH_B` 8000.0 [-200000.0, 200000.0] (fit; P,L; global); `DH_C` 8000.0 [-200000.0, 200000.0] (fit; P,L; global); `DS_B` 0.0 [-500.0, 500.0] (fixed; P,L; global); `DS_C` 0.0 [-500.0, 500.0] (fixed; P,L; global); `DH_AB` 65000.0 [-200000.0, 200000.0] (fit; P,L; global); `DS_AB` 0.0 [-500.0, 500.0] (fixed; P,L; global); `DH_BC` 65000.0 [-200000.0, 200000.0] (fit; P,L; global); `DS_BC` 0.0 [-500.0, 500.0] (fixed; P,L; global) | `KAB`, `KBA`, `KBC`, `KCB`, `PA`, `PB`, `PC` (T,P,L; global) | — |
| `3st_fork` | `PB` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PC` 0.02 [0.0, 1.0] (fit; T,P,L; global); `KEX_AB` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_AC` 200.0 [0.0, 1000000.0] (fit; T,P,L; global) | `PA`, `KAB`, `KBA`, `KAC`, `KCA` (T,P,L; global) | — |
| `3st_linear` | `PB` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PC` 0.02 [0.0, 1.0] (fit; T,P,L; global); `KEX_AB` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_BC` 200.0 [0.0, 1000000.0] (fit; T,P,L; global) | `PA`, `KAB`, `KBA`, `KBC`, `KCB` (T,P,L; global) | — |
| `3st_monomer_dimer_tetramer` | `KD1` 1e-06 [5e-324, 1.0] (fit; T; global); `KD2` 1e-06 [5e-324, 1.0] (fit; T; global); `KOFF1` 100.0 [0.0, 1000000.0] (fit; T; global); `KOFF2` 100.0 [0.0, 1000000.0] (fit; T; global) | `C_MONOMER`, `C_DIMER`, `C_TETRAMER`, `KAB`, `KBA`, `KBC`, `KCB`, `PA`, `PB`, `PC` (T,P; global) | — |
| `3st_monomer_dimer_trimer` | `KD1` 1e-06 [5e-324, 1.0] (fit; T; global); `KD2` 1e-06 [5e-324, 1.0] (fit; T; global); `KOFF1` 100.0 [0.0, 1000000.0] (fit; T; global); `KOFF2` 100.0 [0.0, 1000000.0] (fit; T; global) | `C_MONOMER`, `C_DIMER`, `C_TRIMER`, `KAB`, `KBA`, `KAC`, `KCA`, `KBC`, `KCB`, `PA`, `PB`, `PC` (T,P; global) | — |
| `3st_triangle` | `PB` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PC` 0.02 [0.0, 1.0] (fit; T,P,L; global); `KEX_AB` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_AC` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_BC` 200.0 [0.0, 1000000.0] (fit; T,P,L; global) | `PA`, `KAB`, `KBA`, `KAC`, `KCA`, `KBC`, `KCB` (T,P,L; global) | — |
| `4st` | `PB` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PC` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PD` 0.02 [0.0, 1.0] (fit; T,P,L; global); `KEX_AB` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_AC` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_AD` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_BC` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_BD` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_CD` 200.0 [0.0, 1000000.0] (fit; T,P,L; global) | `PA`, `KAB`, `KBA`, `KAC`, `KCA`, `KAD`, `KDA`, `KBC`, `KCB`, `KBD`, `KDB`, `KCD`, `KDC` (T,P,L; global) | — |
| `4st_binding_3_bound_states` | `KD_APP` 1e-06 [5e-324, 1.0] (fit; T; global); `KOFF_AB` 100.0 [0.0, 1000000.0] (fit; T; global); `KEX_BC` 1000.0 [0.0, 1000000.0] (fit; T; global); `KEQ_BC` 1.0 [0.0, 100.0] (fit; T; global); `KEX_CD` 1000.0 [0.0, 1000000.0] (fit; T; global); `KEQ_CD` 1.0 [0.0, 100.0] (fit; T; global) | `KD_AB` (T; global); `C_P`, `C_L`, `C_PL1`, `C_PL2`, `C_PL3`, `C_PL`, `KD_EFF`, `KAB`, `KBA`, `KBC`, `KCB`, `KCD`, `KDC`, `PA`, `PB`, `PC`, `PD` (T,P,L; global) | `KON_AB` (T; global) |
| `4st_binding_partner_2st` | `KOFF_AB` 100.0 [0.0, 1000000.0] (fit; T; global); `KD_AB` 0.001 [5e-324, 1.0] (fit; T; global); `KOFF_AC` 100.0 [0.0, 1000000.0] (fit; T; global); `KD_AC` 0.001 [5e-324, 1.0] (fit; T; global); `KEQ_L` 1.0 [0.0, 1000000.0] (fit; T; global); `KEQ_PL` 1.0 [0.0, 1000000.0] (fit; T; global); `KEX_BC` 1000.0 [0.0, 1000000.0] (fit; T; global); `KEX_CD` 1000.0 [0.0, 1000000.0] (fit; T; global) | `P_FREE`, `L1_FREE`, `L2_FREE`, `PL1`, `PL2`, `PL3`, `KAB`, `KBA`, `KAC`, `KCA`, `KBC`, `KCB`, `KCD`, `KDC`, `PA`, `PB`, `PC`, `PD` (T,P,L; global) | `KON_AB`, `KON_AC` (T; global) |
| `4st_eyring` | `DH_B` 8000.0 [-200000.0, 200000.0] (fit; P,L; global); `DH_C` 8000.0 [-200000.0, 200000.0] (fit; P,L; global); `DH_D` 8000.0 [-200000.0, 200000.0] (fit; P,L; global); `DH_AB` 75000.0 [-200000.0, 200000.0] (fit; P,L; global); `DH_AC` 75000.0 [-200000.0, 200000.0] (fixed; P,L; global); `DH_AD` 75000.0 [-200000.0, 200000.0] (fixed; P,L; global); `DH_BC` 75000.0 [-200000.0, 200000.0] (fit; P,L; global); `DH_BD` 75000.0 [-200000.0, 200000.0] (fit; P,L; global); `DH_CD` 75000.0 [-200000.0, 200000.0] (fit; P,L; global); `DS_B` 0.0 [-500.0, 500.0] (fixed; P,L; global); `DS_C` 0.0 [-500.0, 500.0] (fixed; P,L; global); `DS_D` 0.0 [-500.0, 500.0] (fixed; P,L; global); `DS_AB` 0.0 [-500.0, 500.0] (fixed; P,L; global); `DS_AC` 0.0 [-500.0, 500.0] (fixed; P,L; global); `DS_AD` 0.0 [-500.0, 500.0] (fixed; P,L; global); `DS_BC` 0.0 [-500.0, 500.0] (fixed; P,L; global); `DS_BD` 0.0 [-500.0, 500.0] (fixed; P,L; global); `DS_CD` 0.0 [-500.0, 500.0] (fixed; P,L; global) | `KAB`, `KBA`, `KAC`, `KCA`, `KAD`, `KDA`, `KBC`, `KCB`, `KBD`, `KDB`, `KCD`, `KDC`, `PA`, `PB`, `PC`, `PD` (T,P,L; global) | — |
| `4st_fork` | `PB` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PC` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PD` 0.02 [0.0, 1.0] (fit; T,P,L; global); `KEX_AB` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_AC` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_AD` 200.0 [0.0, 1000000.0] (fit; T,P,L; global) | `PA`, `KAB`, `KBA`, `KAC`, `KCA`, `KAD`, `KDA` (T,P,L; global) | — |
| `4st_hd` | `D2O` 0.2 [0.0, 1.0] (fixed; D; global); `POP_B` 0.02 [0.0, 1.0] (fit; T,P,L; global); `KEX_AB` 0.0 [0.0, 1000000.0] (fit; T,P,L; global); `KDH_A` 1.0 [0.0, 1000000.0] (fit; T; group); `KDH_B` 1.0 [0.0, 1000000.0] (fit; T; group); `PHI_A` 1.1 [0.75, 1.5] (fixed; T; group) | `PHI_B` (T; group); `KAB`, `KBA`, `KCD`, `KDC` (T,P,L; global); `KAC`, `KCA`, `KBD`, `KDB` (T,D; group); `PA`, `PB`, `PC`, `PD` (T,P,L,D; group) | — |
| `4st_linear` | `PB` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PC` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PD` 0.02 [0.0, 1.0] (fit; T,P,L; global); `KEX_AB` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_BC` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_CD` 200.0 [0.0, 1000000.0] (fit; T,P,L; global) | `PA`, `KAB`, `KBA`, `KBC`, `KCB`, `KCD`, `KDC` (T,P,L; global) | — |
| `5st` | `PB` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PC` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PD` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PE` 0.02 [0.0, 1.0] (fit; T,P,L; global); `KEX_AB` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_AC` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_AD` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_AE` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_BC` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_BD` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_BE` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_CD` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_CE` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_DE` 200.0 [0.0, 1000000.0] (fit; T,P,L; global) | `PA`, `KAB`, `KBA`, `KAC`, `KCA`, `KAD`, `KDA`, `KAE`, `KEA`, `KBC`, `KCB`, `KBD`, `KDB`, `KBE`, `KEB`, `KCD`, `KDC`, `KCE`, `KEC`, `KDE`, `KED` (T,P,L; global) | — |
| `5st_fork` | `PB` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PC` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PD` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PE` 0.02 [0.0, 1.0] (fit; T,P,L; global); `KEX_AB` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_AC` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_AD` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_AE` 200.0 [0.0, 1000000.0] (fit; T,P,L; global) | `PA`, `KAB`, `KBA`, `KAC`, `KCA`, `KAD`, `KDA`, `KAE`, `KEA` (T,P,L; global) | — |
| `5st_linear` | `PB` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PC` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PD` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PE` 0.02 [0.0, 1.0] (fit; T,P,L; global); `KEX_AB` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_BC` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_CD` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_DE` 200.0 [0.0, 1000000.0] (fit; T,P,L; global) | `PA`, `KAB`, `KBA`, `KBC`, `KCB`, `KCD`, `KDC`, `KDE`, `KED` (T,P,L; global) | — |
| `6st` | `PB` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PC` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PD` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PE` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PF` 0.02 [0.0, 1.0] (fit; T,P,L; global); `KEX_AB` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_AC` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_AD` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_AE` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_AF` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_BC` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_BD` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_BE` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_BF` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_CD` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_CE` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_CF` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_DE` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_DF` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_EF` 200.0 [0.0, 1000000.0] (fit; T,P,L; global) | `PA`, `KAB`, `KBA`, `KAC`, `KCA`, `KAD`, `KDA`, `KAE`, `KEA`, `KAF`, `KFA`, `KBC`, `KCB`, `KBD`, `KDB`, `KBE`, `KEB`, `KBF`, `KFB`, `KCD`, `KDC`, `KCE`, `KEC`, `KCF`, `KFC`, `KDE`, `KED`, `KDF`, `KFD`, `KEF`, `KFE` (T,P,L; global) | — |
| `6st_fork` | `PB` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PC` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PD` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PE` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PF` 0.02 [0.0, 1.0] (fit; T,P,L; global); `KEX_AB` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_AC` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_AD` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_AE` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_AF` 200.0 [0.0, 1000000.0] (fit; T,P,L; global) | `PA`, `KAB`, `KBA`, `KAC`, `KCA`, `KAD`, `KDA`, `KAE`, `KEA`, `KAF`, `KFA` (T,P,L; global) | — |
| `6st_linear` | `PB` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PC` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PD` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PE` 0.02 [0.0, 1.0] (fit; T,P,L; global); `PF` 0.02 [0.0, 1.0] (fit; T,P,L; global); `KEX_AB` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_BC` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_CD` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_DE` 200.0 [0.0, 1000000.0] (fit; T,P,L; global); `KEX_EF` 200.0 [0.0, 1000000.0] (fit; T,P,L; global) | `PA`, `KAB`, `KBA`, `KBC`, `KCB`, `KCD`, `KDC`, `KDE`, `KED`, `KEF`, `KFE` (T,P,L; global) | — |
<!-- kinetic-parameters:end -->

</details>
