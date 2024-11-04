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
| `2st_binding` | 2-state exchange model for ligand binding studies                               |
| `4st_hd`      | 4-state exchange model for simultaneous normal and H/D solvent exchange studies |

In these models, each state in the exchange process is represented with a unique parameter suffix (`A`, `B`, `C`, `D`, etc.). For example, `R1_A` denotes the R<sub>1</sub> relaxation rate of the major (ground) state, while `R2_B` refers to the R<sub>2</sub> rate of the first minor state, and so forth.

:::note
For any kinetic model, you can add the `.mf` suffix to create a model that fits model-free parameters directly (e.g., `TAUC_A`, `S2_A`), rather than individual relaxation parameters (e.g., `R1_A`, `R2_A`). For an example, see `CEST_15N_TR/` under `Examples/Experiments/`.
:::
