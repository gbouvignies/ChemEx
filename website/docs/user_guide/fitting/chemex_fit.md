---
sidebar_position: 1
---

# Starting a fit

## Command

To initiate a fit of NMR chemical exchange dataset(s), run the command
`chemex fit` from the command-line followed by a set of options to set and
control the fitting process.

```shell
chemex fit <options>
```

## Options

The list of available options is given in the table below. More details about
the options and their associated files are given in following sections. Some of
the option names of the table are clickable, so that you can reach the related
section directly. Let's dive in.

| Name                                                                  | Description                                                                    |
| --------------------------------------------------------------------- | ------------------------------------------------------------------------------ |
| [`-e`](experiment_files.md) or [`--experiments`](experiment_files.md) | Specify the files containing the experimental setup and data location.         |
| [`-p`](parameter_files.md) or [`--parameters`](parameter_files.md)    | Specify the files containing the initial values of fitting parameters.         |
| [`-m`](method_files.md) or [`--methods`](method_files.md)             | Specify the file containing the fitting method (optional).                     |
| [`-d`](kinetic_models.md) or [`--model`](kinetic_models.md)           | Specify the kinetic model used to fit the datasets (optional, default: `2st`). |
| [`-o`](outputs.md) or [`--output`](outputs.md)                        | Specify the output directory (optional, default: `./Output`).                  |
| [`--plot {nothing,normal,all}`](method_files.md#plotting)             | Plotting level (optional, default: `normal`).                                  |
| `--include`                                                           | Residue(s) to include in the fit (optional).                                   |
| `--exclude`                                                           | Residue(s) to exclude from the fit (optional).                                 |

:::note

`--experiments` and `--parameters` options are mandatory.

:::

:::tip

For arguments requiring input files (e.g., `-p`), it is possible to provide
multiple files by using the wildcard character (i.e., `*`), instead of
specifying the name of each file individually.

:::

:::important TOML File formats

The input and output files of ChemEx use the
[TOML](https://en.wikipedia.org/wiki/TOML) file format. You can find a detailed
description of this file format on the [TOML website](https://toml.io/).

:::

## Combining multiple experiments

ChemEx offers a very simple way to jointly analyse multiple experiments. For
each experiment, simply add the corresponding
[experiment file](experiment_files.md) right after the `--experiments` (or `-e`)
option in the command-line.

This is useful, for example, when you want to fit CEST experiments with
different B<sub>1</sub> field together, or CPMG relaxation dispersion
experiments recorded at different B<sub>0</sub> field together, or even a
combination of the two.

An example illustrating
[the study of protein-ligand binding with CPMG and CEST experiments](examples/binding.md)
is available
[here](https://github.com/gbouvignies/chemex/tree/master/examples/Combinations/2stBinding/).

## Example

A typical command line invocation would be:

```shell
chemex fit -e Experiments/*.toml \
           -p Parameters/parameters.toml \
           -m Methods/method.toml \
           -o Output
```

Output files are located in the directory specified with the option `-o`. Plots
showing the line of best fit are also produced and saved to the disk by default
– you can modulate this behaviour with the option `--plot`.

:::tip

You can save the command in a shell script (typically called `run.sh` in most
examples) to save some typing efforts.

```shell title=run.sh
#!/bin/sh

chemex fit -e Experiments/*.toml \
           -p Parameters/parameters.toml \
           -m Methods/method.toml \
           -o Output
```

Make it executable with the following command:

```shell
chmod +x run.sh
```

:::
