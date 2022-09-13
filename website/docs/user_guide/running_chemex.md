---
sidebar_position: 2
lastmod: 2022-04-19T20:14:49.735Z
---

# Running ChemEx

ChemEx is a command-line application. After the installation, you can run the
program directly from the shell prompt using the command `chemex`:

```shell
chemex <subcommand> <options>
```

Four **sub-commands** are available, each corresponding to a specific task:

- [`fit`](fitting/chemex_fit.md) starts the fits of experimental datasets: if
  you are interested in ChemEx, this is likely the sub-command you are looking
  for.
- `simulate` executes the module, which allows to calculate synthetic profiles
  from a given experiment and a given model.
- `pick_cest` is a small GUI application, which plots CEST, cos-CEST and D-CEST
  profiles for interactive dip picking and obtaining starting minor state
  position values.
- `plot_param` produces plots of selected parameters based on fitting output
  results.

:::tip

You can obtain the list of available sub-commands using the `--help` option.

```bash
chemex --help
```

You can display the list of options associated with each sub-command using the
sub-command name followed by the `--help` option.

```bash
chemex <subcommand> --help
```

:::
