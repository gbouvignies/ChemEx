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
| `2st_rs`      | 2-state exchange model for residue-specific studies                             |
| `2st_hd`      | 2-state exchange model for H/D solvent exchange studies                         |
| `2st_eyring`  | 2-state exchange model for temperature-dependent studies                        |
| `3st_eyring`  | 3-state exchange model for temperature-dependent studies                        |
| `4st_eyring`  | 4-state exchange model for temperature-dependent studies                        |
| `2st_binding` | 2-state exchange model for ligand binding studies                               |
| `4st_hd`      | 4-state exchange model for simultaneous normal and H/D solvent exchange studies |

In these models, each state in the exchange process is represented with a unique parameter suffix (`A`, `B`, `C`, `D`, etc.). For example, `R1_A` denotes the R<sub>1</sub> relaxation rate of the major (ground) state, while `R2_B` refers to the R<sub>2</sub> rate of the first minor state, and so forth.

:::note
For any kinetic model, you can add the `.mf` suffix to create a model that fits model-free parameters directly (e.g., `TAUC_A`, `S2_A`), rather than individual relaxation parameters (e.g., `R1_A`, `R2_A`). For an example, see `CEST_15N_TR/` under `Examples/Experiments/`.
:::

## Temperature-Dependent Eyring Models

The `2st_eyring`, `3st_eyring`, and `4st_eyring` models implement 2-state, 3-state, and 4-state exchange systems (respectively) with temperature-dependent rate constants calculated using Eyring transition state theory. These models are particularly useful for studying exchange processes where thermodynamic parameters govern the temperature dependence of exchange rates.

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

### 3st_eyring Model Parameters

The `3st_eyring` model extends the 2-state model with an additional state:

**State Energies (relative to state A):**
- `DH_B`, `DH_C`: Enthalpy differences (J/mol) for states B, C
- `DS_B`, `DS_C`: Entropy differences (J/mol/K) for states B, C

**Transition Barriers:**
- `DH_AB`, `DH_AC`, `DH_BC`: Activation enthalpies (J/mol) for transitions
- `DS_AB`, `DS_AC`, `DS_BC`: Activation entropies (J/mol/K) for transitions

The model automatically calculates all 6 rate constants (k_AB, k_BA, k_AC, k_CA, k_BC, k_CB).

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

