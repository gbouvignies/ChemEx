---
sidebar_position: 1
---

# Starting a Fit

## Command

To initiate a fit of NMR chemical exchange datasets, use the `chemex fit` command in the command line. This command can be accompanied by various options to customize the fitting process:

```shell
chemex fit <options>
```

## Options

The following table lists the options for `chemex fit`. Each option's name is a link to its detailed explanation, and these options are further elaborated in subsequent sections.

| Name | Description |
| ---- | ----------- |
| [`-e`](experiment_files.md) or [`--experiments`](experiment_files.md) | Specify files containing experimental setup and data. |
| [`-p`](parameter_files.md) or [`--parameters`](parameter_files.md) | Specify files containing initial fitting parameters. |
| [`-m`](method_files.md) or [`--methods`](method_files.md) | Indicate the fitting method file (optional). |
| [`-d`](kinetic_models.md) or [`--model`](kinetic_models.md) | Specify the kinetic model for fitting (optional, default: `2st`). |
| [`-o`](outputs.mdx) or [`--output`](outputs.mdx) | Set the output directory (optional, default: `./Output`). |
| [`--plot {nothing,normal,all}`](method_files.md#plotting) | Select the plotting level (optional, default: `normal`). |
| `--include` | Define residues to include in the fit (optional). |
| `--exclude` | Define residues to exclude from the fit (optional). |

:::note

The `--experiments` and `--parameters` options are mandatory.

:::

:::tip

For file-input arguments (e.g., `-p`), you can use the wildcard character (`*`) to collectively input multiple files instead of listing each file individually.

:::

:::important TOML File Formats

ChemEx uses the [TOML](https://toml.io/) file format for its input and output files. You can find detailed information about this format on the [TOML website](https://en.wikipedia.org/wiki/TOML).

:::

## Combining Multiple Experiments

ChemEx facilitates the combined analysis of multiple experiments. To include various experiments in a fit, add their corresponding [experiment files](experiment_files.md) after the `--experiments` (or `-e`) option in the command line. This is particularly useful for fitting different types of CEST or CPMG experiments, or a combination thereof.

See an [example of protein-ligand binding analysis](https://github.com/gbouvignies/chemex/tree/master/examples/Combinations/2stBinding/) using both CPMG and CEST experiments [here](examples/binding.md).

## Example

A typical command for invoking ChemEx is as follows:

```shell
chemex fit -e Experiments/*.toml \
           -p Parameters/parameters.toml \
           -m Methods/method.toml \
           -o Output
```

The output files are saved in the directory specified by `-o`. The default setting includes the generation of plots to illustrate the best fit lines, which can be adjusted using the `--plot` option.

:::tip

For convenience, consider saving this command in a shell script, commonly named `run.sh` in examples:

```shell title=run.sh
#!/bin/sh

chemex fit -e Experiments/*.toml \
           -p Parameters/parameters.toml \
           -m Methods/method.toml \
           -o Output
```

To make the script executable, use:

```shell
chmod +x run.sh
```

:::
