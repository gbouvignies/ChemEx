---
sidebar_position: 7
---

# Kinetic Models

The kinetic model (specified with the `-d` or `--model` option) defines the type of exchange model to be used for data analysis. Available models include:

| Model Name          | Description                                                                     |
| ------------------- | ------------------------------------------------------------------------------- |
| `2st`               | 2-state exchange model (default)                                                |
| `3st`               | Complete 3-state exchange model                                                  |
| `3st_triangle`      | Historical explicit compatibility name for `3st`                                |
| `3st_linear`        | Linear A ↔ B ↔ C exchange model                                                  |
| `3st_fork`          | A-centered B ↔ A ↔ C exchange model                                              |
| `4st`               | Complete 4-state exchange model                                                  |
| `4st_linear`        | Linear 4-state exchange model                                                    |
| `4st_fork`          | A-centered 4-state exchange model                                                |
| `5st`               | Complete 5-state exchange model                                                  |
| `5st_linear`        | Linear 5-state exchange model                                                    |
| `5st_fork`          | A-centered 5-state exchange model                                                |
| `6st`               | Complete 6-state exchange model                                                  |
| `6st_linear`        | Linear 6-state exchange model                                                    |
| `6st_fork`          | A-centered 6-state exchange model                                                |
| `2st_hd`            | 2-state exchange model for H/D solvent exchange studies                         |
| `2st_eyring`        | 2-state exchange model for temperature-dependent studies                        |
| `3st_eyring`        | Compatibility name for the linear 3-state Eyring model                          |
| `3st_eyring_linear` | Linear A ↔ B ↔ C Eyring model for temperature-dependent studies                 |
| `3st_eyring_fork`   | Fork B ↔ A ↔ C Eyring model for temperature-dependent studies                   |
| `4st_eyring`        | 4-state exchange model for temperature-dependent studies                        |
| `2st_binding`       | 2-state exchange model for ligand binding studies                               |
| `3st_binding_cs`    | 3-state conformational-selection ligand binding model                           |
| `3st_binding_if`    | 3-state induced-fit ligand binding model                                         |
| `4st_hd`            | 4-state exchange model for simultaneous normal and H/D solvent exchange studies |

In these models, each state in the exchange process is represented with a unique
parameter suffix (`A`, `B`, `C`, `D`, etc.). For example, `R1_A` and `R2_B`
refer to the R<sub>1</sub> rate of state A and the R<sub>2</sub> rate of state B.
In many analyses A is assigned to the major or ground state and B to a minor or
excited state, but this is an analyst-defined convention rather than an ordering
enforced by ChemEx. See [Exchange States and Parameters](exchange_states_parameters.md)
for populations, directional rates, shift signs, and label swapping.

:::note
For any kinetic model, you can add the `.rs` suffix to make the kinetic parameters residue-specific (for example, `2st.rs` or `3st_eyring.rs`). Suffixes can be combined, such as `2st.rs.mf`. The legacy `2st_rs` name remains supported as an alias for `2st.rs`.
:::

:::note
For any kinetic model, you can add the `.mf` suffix to create a model that fits model-free parameters directly (e.g., `TAUC_A`, `S2_A`), rather than individual relaxation parameters (e.g., `R1_A`, `R2_A`). For an example, see `CEST_15N_TR/` under `Examples/Experiments/`.
:::

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
