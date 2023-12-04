---
sidebar_position: 1
---

# Running ChemEx

ChemEx operates as a command-line tool. Once installed, execute it from the shell prompt with the `chemex` command:

```shell
chemex <subcommand> <options>
```

It offers four **sub-commands**, each dedicated to a distinct functionality:

- [`fit`](fitting/chemex_fit.md): Initiates the fitting of experimental data. This is the primary sub-command for most ChemEx users.
- `simulate`: Runs the module to generate synthetic profiles using specified experiments and models.
- `pick_cest`: Provides a small GUI for plotting CEST, cos-CEST, and D-CEST profiles. It facilitates interactive dip picking and initial estimation of minor state positions.
- `plot_param`: Generates plots for selected parameters derived from fitting results.

:::tip

Retrieve a list of all sub-commands by using the `--help` option:

```bash
chemex --help
```

For details on each sub-command's options, append `--help` after the sub-command:

```bash
chemex <subcommand> --help
```

:::