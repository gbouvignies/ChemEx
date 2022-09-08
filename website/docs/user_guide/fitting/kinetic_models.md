---
sidebar_position: 7
---

# Kinetic models

The kinetic model (indicated with `-d` or `--model`) indicates the type of
exchange model to be used for data analysis, which can be one of the following:

| Model Name    | Description                                                                   |
| ------------- | ----------------------------------------------------------------------------- |
| `2st`         | 2-state exchange model (default)                                              |
| `3st`         | 3-state exchange model                                                        |
| `4st`         | 4-state exchange model                                                        |
| `2st_rs`      | 2-state exchange model for residue-specific study                             |
| `2st_hd`      | 2-state exchange model for H/D solvent exchange study                         |
| `2st_eyring`  | 2-state exchange model for temperature-dependent study                        |
| `3st_eyring`  | 3-state exchange model for temperature-dependent study                        |
| `2st_binding` | 2-state exchange model for ligand binding study                               |
| `4st_hd`      | 4-state exchange model for simultaneous normal and H/D solvent exchange study |

Multiple states involved in an exchange process are distinguished with different
parameter suffix `A`, `B`, `C`, `D`, etc. For example, `R1_A` represents
R<sub>1</sub> rate of the major state (i.e., ground state), `R2_B` represents
R<sub>2</sub> rate of the first minor state, and so on.

:::note

For each kinetic model, it is possible to add `.mf` suffix to create a new
kinetic model (e.g., `2st.mf`). In such kinetic models, model-free parameters
(e.g., `TAUC_A`, `S2_A`, etc.) are directly fitted instead of each individual
relaxation parameter (e.g., `R1_A`, `R2_A`, etc.). See `CEST_15N_TR/` under
`Examples/Experiments/` as an example.

:::
