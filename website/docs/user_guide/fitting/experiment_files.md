---
sidebar_position: 3
---

# Experiment Files

## Overview

In ChemEx, experiment files specify details about the datasets and file locations for analysis. These files are provided to ChemEx through the `-e` or `--experiments` command-line option:

```shell
chemex fit -e <experiment_file1> <experiment_file2> [...]
```

## Structure of Experiment Files

Experiment files consist of three main sections:

1. **`[experiment]`**: Contains information about the experiment setup.
2. **`[conditions]`**: Specifies sample conditions.
3. **`[data]`**: Details data location, spin system assignments, and methods for estimating uncertainties or filtering data.

### Example of an Experiment File

```toml title="experiment.toml"
[experiment]
name         = "cest_15n"
time_t1      = 0.5
carrier      = 118.987
b1_frq       = 26.3

[conditions]
h_larmor_frq = 499.243
# p_total     = 2.0e-3
# sample      = "A39G FF domain"
# temperature = 1.0

[data]
path           = "../Data/26Hz/"
error          = "scatter"
filter_offsets = [[0.0, 26.0]]

  [data.profiles]
  13N = "13N-HN.out"
  26N = "26N-HN.out"
  28N = "28N-HN.out"
  29N = "29N-HN.out"
```

:::note
Each experiment available in ChemEx includes a sample configuration file, which details available key-value pairs. These samples are found in the [Experiments](/docs/experiments) section.
:::

## Detailed Sections

### `[experiment]`

The `[experiment]` section defines the type and settings of the pulse sequence. Common keys include:

| Key               | Description                                                                                           |
| ----------------- | ----------------------------------------------------------------------------------------------------- |
| `name`            | Pulse sequence name.                                                                                  |
| `carrier`         | Carrier position during the experiment, in ppm.                                                       |
| `time_t2`, `time_t1` | Relaxation delays in seconds (e.g., for CPMG relaxation dispersion experiments).                   |
| `pw90`            | 90-degree pulse width, in seconds.                                                                    |
| `b1_frq`          | B1 radio-frequency field strength, in Hz.                                                             |
| `observed_state`  | ID of the observed state (e.g., `a`, `b`, `c`, or `d`).                                              |

### `[conditions]`

The `[conditions]` section provides experimental and sample conditions, such as Larmor frequency, temperature, concentration, and labeling. `h_larmor_frq` (Larmor frequency) is required, while `temperature`, `p_total`, and `l_total` depend on the kinetic model used.

| Key                 | Description                                                                                           |
| ------------------- | ----------------------------------------------------------------------------------------------------- |
| `h_larmor_frq`      | Larmor frequency in MHz.                                                                              |
| `temperature`       | Sample temperature in °C (optional, required by certain kinetic models).                              |
| `p_total`, `l_total`| Protein and ligand concentrations in M (optional, required by certain kinetic models).                |
| `label`             | Labeling scheme for the sample.                                                                       |

#### `label`

The `label` key indicates the sample’s isotopic labeling scheme, used to account for isotopic effects in certain experiments. Examples:

- `"2H"` for deuterated samples.
- `"13C"` for uniformly <sup>13</sup>C-labeled samples.

For a uniformly <sup>13</sup>C-labeled, perdeuterated sample:

```toml
label = ["13C", "2H"]
```

### `[data]`

The `[data]` section specifies data file locations, spin system assignments, and methods for error estimation and data filtering.

| Key              | Description                                                                                             |
| ---------------- | ------------------------------------------------------------------------------------------------------- |
| `path`           | Path to the directory containing data files.                                                            |
| `error`          | Method for error estimation.                                                                            |
| `filter_offsets` | List of offsets to exclude from calculations (e.g., in CEST experiments).                               |
| `filter_planes`  | List of planes to exclude from calculations (e.g., in CEST and CPMG experiments).                       |
| `[data.profiles]`| Subsection listing experimental profile file names with their spin-system assignments.                  |

#### `error`

The `error` key specifies the method for estimating uncertainties:

| Value         | Description                                                                                                                               |
| ------------- | ----------------------------------------------------------------------------------------------------------------------------------------- |
| `"file"`      | Use uncertainties directly from the data file.                                                                                            |
| `"duplicates"`| Calculate uncertainty using pooled standard deviation, or return average error if duplicates are not available.                           |
| `"scatter"`   | Estimate uncertainty by assuming additive Gaussian noise on the profile, suitable for CEST data.                                          |

#### `filter_offsets`

This key filters out specified offsets from CEST profiles, useful for removing artifacts like sidebands. Provide a list of offset and bandwidth pairs, where each pair specifies an offset relative to the main resonance and a bandwidth around the offset to exclude. Values are in Hz.

Example:

```toml
filter_offsets = [
   [0.0, 20.0],
   [-300.0, 20.0],
   [+300.0, 20.0],
 ]
```

#### `[data.profiles]`

The `[data.profiles]` subsection lists filenames of experimental profiles along with their spin-system assignments. Spin-system names follow Sparky-NMR conventions, using a group name (e.g., amino acid and position, like ALA3 or A3) and atom name (e.g., N, CA, CG1).

For multi-spin systems, use the `-` sign to separate spins. For example, `G23N-G23H` or `G23N-H` are both valid.

Example:

```toml
  [data.profiles]
  13N = "13N-HN.out"
  26N = "26N-HN.out"
  28N = "28N-HN.out"
  29N = "29N-HN.out"
```

:::note
Select spin-system names that accurately reflect the spin system of interest in each experiment. While ChemEx can auto-correct minor name inconsistencies, precise naming is recommended.
:::
