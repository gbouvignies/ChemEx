---
sidebar_position: 1
---

# Running ChemEx

ChemEx is a command-line tool. Once installed, you can execute it directly from the shell prompt using the `chemex` command:

```shell
chemex <subcommand> <options>
```

ChemEx provides four primary **sub-commands**, each dedicated to a specific function:

- [`fit`](fitting/chemex_fit.md): Initiates the fitting of experimental data. This is the main sub-command for most ChemEx users.
- `simulate`: Generates synthetic profiles based on specified experiments and models.
- `pick_cest`: Opens a small GUI for plotting CEST, cos-CEST, and D-CEST profiles, enabling interactive peak picking and initial estimation of minor state positions.
- `plot_param`: Creates plots for selected parameters derived from fitting results.

:::tip

To view a list of all sub-commands, use the `--help` option:

```bash
chemex --help
```

For detailed options specific to each sub-command, append `--help` after the sub-command:

```bash
chemex <subcommand> --help
```

:::
