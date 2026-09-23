---
sidebar_position: 7
---

# Temperature-dependent chemical shifts

The `.tc` model suffix gives chemical shifts a shared linear temperature
dependence. It is independent of the kinetic model: `2st.tc`, `4st.tc`,
`2st_eyring.tc`, and compositions with `.mf` or `.rs` all use the same shift
parameterization. Every Experiment in a `.tc` analysis must specify
`temperature` under `[conditions]` in degrees Celsius.

Select it through the ordinary model option, for example:

```bash
chemex fit -e Experiments/*.toml -p Parameters.toml -m Method.toml -d 2st.tc
```

`.tc` changes chemical shifts only; `2st.tc` leaves the two-state exchange
kinetics and populations on their ordinary model. To combine shift-temperature
coefficients with Eyring temperature-dependent kinetics, compose the features
explicitly:

```bash
chemex fit -e Experiments/*.toml -p Parameters.toml -m Method.toml -d 2st_eyring.tc
```

## Model and units

With `dT = T - TREF`, ChemEx resolves

$$
CS_A(T) = CS0_A + CS1_A\,dT,
$$

$$
DW_{AX}(T) = DW0_{AX} + DW1_{AX}\,dT,
$$

$$
CS_X(T) = CS_A(T) + DW_{AX}(T).
$$

`CS0_A` and `DW0_AX` are in ppm and are the respective values at `TREF`.
`CS1_A` and `DW1_AX` are in ppm/°C. A Celsius difference is numerically equal
to a kelvin difference, so `.tc` composes consistently with Eyring kinetic
models: Eyring calculations convert absolute temperatures to kelvin, while the
shift polynomial uses the Celsius difference.

The coefficient names are order-indexed, but `.tc` implements only the linear
orders 0 and 1.

| Parameter | Meaning | Default | Bounds | Default role |
| --- | --- | ---: | ---: | --- |
| `TREF` | common reference temperature | 25 °C | finite and above −273.15 °C | protected fixed constant |
| `CS0_A` | state-A shift at `TREF` | 0 ppm | [−100, 300] ppm | Experiment Type authority |
| `CS1_A` | state-A linear coefficient | 0 ppm/°C | [−1, 1] ppm/°C | Experiment Type authority |
| `DW0_AX` | A-to-X shift difference at `TREF` | 0 ppm | [−100, 100] ppm | fitted |
| `DW1_AX` | A-to-X linear coefficient | 0 ppm/°C | [−1, 1] ppm/°C | fitted |

Direct chemical-shift Experiment Types already declare the absolute reference
shift fit-capable. Under `.tc`, that declaration expands to both `CS0_A` and
`CS1_A`. CEST and CPMG Experiment Types do not start fitting absolute shifts
merely because `.tc` is selected; a Method Plan can explicitly select those
coefficients. The shift-difference coefficients keep the established `.tc`
default fitted behavior.

## Configuring the reference temperature

An analysis has exactly one unqualified model constant:

```toml
[GLOBAL]
TREF = 20.0
```

`TREF` is shared across all residues, nuclei, fields, concentrations, Experiment
Types, and temperatures in that analysis. It cannot be residue-, nucleus-, or
condition-scoped, and it cannot have fitting bounds or a grid step. Method Plan
`FIT`, `FIX`, `CONSTRAIN`, `GRID`, and DE-coordinate operations cannot target it.

The coefficients retain their normal spin/nucleus and state-pair scope but omit
temperature and magnetic field. Consequently one coefficient line is shared by
all temperatures and fields unless the user applies an explicit supported
constraint or sharing rule.

## Multistate semantics

ChemEx keeps state A as the chemical-shift reference. A multistate model has one
`CS0_A`/`CS1_A` pair for state A and a `DW0_AX`/`DW1_AX` pair for each non-A
state X. State X is always formed as `CS_A + DW_AX`; there are no independent
`CS0_X` or `CS1_X` coordinates. This prevents the reference-state contribution
from being counted twice.

The temperature-specific `CS_A`, `DW_AX`, and `CS_X` entries remain derived
model quantities. Select, initialize, fix, fit, share, or constrain the
order-indexed coefficients instead of overriding those derived values.

## Reference-temperature invariance

Changing the reference from `r` to `r'` preserves the physical line when each
coefficient pair is transformed as

$$
Q0' = Q0 + Q1(r' - r), \qquad Q1' = Q1,
$$

where Q is either `CS` or `DW`. ChemEx does not automatically rewrite arbitrary
fitted parameter files between reference temperatures. Independent coefficient
bounds, priors, grid ranges, and constraints generally do **not** retain the same
meaning under this coordinate transformation. Change `TREF` only together with
an explicit scientific review of those settings.

## Fitting limitation

At one temperature, a free order-0 coefficient and its free order-1 coefficient
enter observables only through one linear combination and cannot be estimated
independently without an additional constraint. Use data at multiple
temperatures, or explicitly fix or constrain one coefficient. ChemEx does not
silently change the fitting roles for a single-temperature analysis.

## Breaking change from the previous `.tc` interface

The released `.tc` implementation used an uncentered line,
`DW_AX(T) = DWP_AX + DWM_AX T`, with `T` in °C. This release intentionally
replaces `DWP_AX` and `DWM_AX` with `DW0_AX` and `DW1_AX`; the old names are not
aliases and are not migrated automatically. Existing `.tc` parameter and Method
files must be rewritten in the centered coordinate system:

$$
DW0_{AX} = DWP_{AX} + TREF\,DWM_{AX}, \qquad DW1_{AX} = DWM_{AX}.
$$

Review bounds and constraints when converting because changing the reference
temperature changes their coordinate meaning.

## Output, provenance, and restart

New parameter output contains `TREF`, `CS0`/`CS1`, and `DW0`/`DW1` only. `TREF`
is written as a fixed global value whenever its associated shift polynomial is
active. Temperature-specific resolved shifts may appear as derived output, but
they are not independent restart coordinates. Run provenance archives the input
files as supplied; restart from ChemEx's canonical output so the reference and
coefficients remain together.
