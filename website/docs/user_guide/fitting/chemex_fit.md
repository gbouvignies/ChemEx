---
sidebar_position: 1
---

# Starting a Fit

## Command

To fit NMR chemical exchange datasets, use the `chemex fit` command from the command line. This command accepts various options to customize the fitting process:

```shell
chemex fit <options>
```

## Options

The table below lists the available options for `chemex fit`. Each option links to a detailed explanation, which you’ll find in subsequent sections.

| Name | Description |
| ---- | ----------- |
| [`-e`](experiment_files.md) or [`--experiments`](experiment_files.md) | Specify files containing experimental setup and data. |
| [`-p`](parameter_files.md) or [`--parameters`](parameter_files.md) | Specify files containing initial fitting parameters. |
| [`-m`](method_files.md) or [`--methods`](method_files.md) | Indicate the fitting method file (optional). |
| [`-d`](kinetic_models.md) or [`--model`](kinetic_models.md) | Specify the kinetic model for fitting (optional, default: `2st`). |
| [`-o`](outputs.mdx) or [`--output`](outputs.mdx) | Set the output directory (optional, default: `./Output`). |
| `--plot {nothing, normal, all}` | Select the plotting level (optional, default: `normal`). |
| `--include` | Define residues to include in the fit (optional). |
| `--exclude` | Define residues to exclude from the fit (optional). |

:::note
The `--experiments` and `--parameters` options are required.
:::

:::tip
For options that accept files (e.g., `-p`), you can use the wildcard character (`*`) to input multiple files collectively, rather than listing each file individually.
:::

:::important TOML File Formats

ChemEx uses the [TOML](https://toml.io/) format for input and output files. For more information about this format, visit the [TOML website](https://en.wikipedia.org/wiki/TOML).
:::

## Combining Multiple Experiments

ChemEx supports combined analysis of multiple experiments. To include various experiments in a fit, list the corresponding [experiment files](experiment_files.md) after the `--experiments` (or `-e`) option. This feature is particularly useful for fitting different types of CEST or CPMG experiments or a combination of both.

For an example of protein-ligand binding analysis using both CPMG and CEST experiments, see [this example](examples/binding.md).

## Example

Here is a typical command for running ChemEx:

```shell
chemex fit -e Experiments/*.toml \
           -p Parameters/parameters.toml \
           -m Methods/method.toml \
           -o Output
```

Output files are saved in the directory specified by `-o`. By default, plots illustrating best-fit lines are generated; you can adjust this behavior with the `--plot` option.

:::tip
For convenience, consider saving this command in a shell script, typically named `run.sh`:

```shell title=run.sh
#!/bin/sh

chemex fit -e Experiments/*.toml \
           -p Parameters/parameters.toml \
           -m Methods/method.toml \
           -o Output
```

To make the script executable, run:

```shell
chmod +x run.sh
```
:::
