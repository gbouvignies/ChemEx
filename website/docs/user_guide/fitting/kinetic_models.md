---
sidebar_position: 7
---

# Kinetic Models

The kinetic model (specified with the `-d` or `--model` option) defines the type of exchange model to be used for data analysis. Available models include:

| Model Name    | Description                                                                     |
| ------------- | ------------------------------------------------------------------------------- |
| `2st`         | 2-state exchange model (default)                                                |
| `3st`         | 3-state exchange model                                                          |
| `4st`         | 4-state exchange model                                                          |
| `2st_hd`      | 2-state exchange model for H/D solvent exchange studies                         |
| `2st_eyring`  | 2-state exchange model for temperature-dependent studies                        |
| `3st_eyring`  | Compatibility name for the linear 3-state Eyring model                          |
| `3st_eyring_linear` | Linear A ↔ B ↔ C Eyring model for temperature-dependent studies          |
| `3st_eyring_fork` | Fork B ↔ A ↔ C Eyring model for temperature-dependent studies                |
| `4st_eyring`  | 4-state exchange model for temperature-dependent studies                        |
| `2st_binding` | 2-state exchange model for ligand binding studies                               |
| `4st_hd`      | 4-state exchange model for simultaneous normal and H/D solvent exchange studies |

In these models, each state in the exchange process is represented with a unique parameter suffix (`A`, `B`, `C`, `D`, etc.). For example, `R1_A` denotes the R<sub>1</sub> relaxation rate of the major (ground) state, while `R2_B` refers to the R<sub>2</sub> rate of the first minor state, and so forth.

:::note
For any kinetic model, you can add the `.rs` suffix to make the kinetic parameters residue-specific (for example, `2st.rs` or `3st_eyring.rs`). Suffixes can be combined, such as `2st.rs.mf`. The legacy `2st_rs` name remains supported as an alias for `2st.rs`.
:::

:::note
For any kinetic model, you can add the `.mf` suffix to create a model that fits model-free parameters directly (e.g., `TAUC_A`, `S2_A`), rather than individual relaxation parameters (e.g., `R1_A`, `R2_A`). For an example, see `CEST_15N_TR/` under `Examples/Experiments/`.
:::

## Temperature-Dependent Eyring Models

The `2st_eyring`, three-state Eyring variants, and `4st_eyring` models
implement exchange systems with temperature-dependent rate constants calculated
using Eyring transition state theory. These models are particularly useful for
studying exchange processes where thermodynamic parameters govern the
temperature dependence of exchange rates.

### Theoretical Background

The Eyring equation relates the rate constant to the activation free energy:

```
k_ij = (k_B * T / h) * exp(-ΔG‡_ij / RT)
```

where:
- `k_ij` is the rate constant for transition from state i to j (s⁻¹)
- `k_B` is Boltzmann's constant (1.380649×10⁻²³ J/K)
- `T` is temperature in Kelvin
- `h` is Planck's constant (6.62607015×10⁻³⁴ J·s)
- `ΔG‡_ij` is the activation free energy (J/mol)
- `R` is the gas constant (8.314462618 J/mol/K)

The activation free energy is calculated from enthalpic and entropic contributions:

```
ΔG‡_ij = ΔH‡_ij - T * ΔS‡_ij
```

### 2st_eyring Model Parameters

The `2st_eyring` model uses the following thermodynamic parameters:

**State Energies (relative to state A):**
- `DH_B`: Enthalpy difference (J/mol) for state B relative to A
- `DS_B`: Entropy difference (J/mol/K) for state B relative to A

**Transition Barriers:**
- `DH_AB`: Activation enthalpy (J/mol) for A → B transition
- `DS_AB`: Activation entropy (J/mol/K) for A → B transition

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

**Linear Transition Barriers (`3st_eyring`, `3st_eyring_linear`):**
- `DH_AB`, `DH_BC`: Activation enthalpies (J/mol) for the AB and BC transitions
- `DS_AB`, `DS_BC`: Activation entropies (J/mol/K) for the AB and BC transitions

The linear models calculate k_AB, k_BA, k_BC, and k_CB. There is no direct A ↔ C
pathway.

**Fork Transition Barriers (`3st_eyring_fork`):**
- `DH_AB`, `DH_AC`: Activation enthalpies (J/mol) for the AB and AC transitions
- `DS_AB`, `DS_AC`: Activation entropies (J/mol/K) for the AB and AC transitions

The fork model calculates k_AB, k_BA, k_AC, and k_CA. There is no direct B ↔ C
pathway.

Parameter files used with `3st_eyring` or `3st_eyring_linear` should not contain
`DH_AC` or `DS_AC`; those parameters are not part of the linear topology.

### 4st_eyring Model Parameters

The `4st_eyring` model implements a full 4-state system:

**State Energies (relative to state A):**
- `DH_B`, `DH_C`, `DH_D`: Enthalpy differences (J/mol) for states B, C, D
- `DS_B`, `DS_C`, `DS_D`: Entropy differences (J/mol/K) for states B, C, D

**Transition Barriers:**
- `DH_AB`, `DH_AC`, `DH_AD`: Activation enthalpies (J/mol) for transitions from A
- `DH_BC`, `DH_BD`, `DH_CD`: Activation enthalpies (J/mol) for transitions between B, C, D
- `DS_AB`, `DS_AC`, `DS_AD`: Activation entropies (J/mol/K) for transitions from A
- `DS_BC`, `DS_BD`, `DS_CD`: Activation entropies (J/mol/K) for transitions between B, C, D

The model automatically calculates all 12 rate constants (k_AB, k_BA, k_AC, k_CA, k_AD, k_DA, k_BC, k_CB, k_BD, k_DB, k_CD, k_DC).

:::note
State A serves as the reference state with ΔH_A = ΔS_A = 0 for all Eyring models. Rate constants are automatically clipped to [0, 1×10¹⁶ s⁻¹] for numerical stability.
:::

All Eyring state and transition enthalpy parameters have default bounds of
[-2×10⁵, 2×10⁵] J/mol. The corresponding entropy parameters have default
bounds of [-5×10², 5×10²] J/mol/K. Explicit bounds in parameter files override
these defaults.
