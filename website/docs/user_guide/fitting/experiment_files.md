---
sidebar_position: 3
---

# Experiment files

## Description

The experiment files containing the experimental details about the datasets to
be fitted as well as the location of the files of the profiles. These files are
provided to ChemEx using the `-e` or `--experiments` option.

```shell
chemex fit -e <experiment_file1> <experiment_file2> [...]
```

Experiment files are divided in 3 different sections:

-   `[experiment]` includes the information about the experiment.
-   `[conditions]` provides the sample conditions.
-   `[data]` contains the details about the data location, spin system assignment
    as well as methods to estimate uncertainties or filter out some measurements.

## Example

Here is an example experiment file:

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

The list of available key-value pairs depends on each experiment. For each
experiment available in ChemEx, a sample configuration file is provided in the
section [Experiments](/docs/experiments). Detailed descriptions about the
meaning of each key are included as comments in sample config files.

:::

## Sections

### `[experiment]`

The `[experiment]` section contains the details about the pulse sequence type
and settings. Below are some of the keys that are commonly found in the
`[experiment]` section:

| Name                 | Description                                                                                           |
| -------------------- | ----------------------------------------------------------------------------------------------------- |
| `name`               | Pulse sequence name.                                                                                  |
| `carrier`            | Position of the carrier during the experiment, in ppm.                                                |
| `time_t2`, `time_t1` | Transverse and longitudinal relaxation delay in second (e.g. CPMG relaxation dispersion experiments). |
| `pw90`               | 90 degree pulse width, in second.                                                                     |
| `b1_frq`             | B1 radio-frequency field strength, in Hz.                                                             |
| `observed_state`     | The ID of the state that is observed (one of `a`, `b`, `c` or `d`).                                   |

### `[conditions]`

The section `[conditions]` provides information about the experimental and
sample conditions (Larmor frequency, sample temperature, protein concentration
and labeling):

| Name                 | Description                                                                                           |
| -------------------- | ----------------------------------------------------------------------------------------------------- |
| `h_larmor_frq`       | Larmor frequency, in MHz.                                                                             |
| `temperature`        | Sample temperature, in °C (optional, required by some the kinetic models).                            |
| `p_total`, `l_total` | Protein and ligand concentration, respectively, in M (optional, required by some the kinetic models). |
| `label`              | Labeling scheme of the sample.                                                                        |

While you should always provide the Larmor frequency (`h_larmor_frq`), the
`temperature`, `p_total` and `l_total` keys are only required depending on which
[kinetic model](kinetic_models.md) is employed.

#### `label`

The `label` key is a list providing information about the labeling scheme of the
sample. The labeling information is used in some experiments to take into
account some isotopic effects that may affect relaxation rates or the presence
of scalar couplings.

-   Use `"2H"` for deuterated samples to obtain accurate initial estimates of
    relaxation rates based on model-free parameters.
-   Use `"13C"` for uniformly <sup>13</sup>C-labeled samples to account for scalar
    couplings in CEST experiments properly.

For example, for a uniformly <sup>13</sup>C-labeled, perdeuterated sample, use
the following configuration:

```toml
label = ["13C", "2H"]
```

### `[data]`

The section `[data]` contains the details about the data location, the
spin-system assignments, and methods to estimate uncertainties and/or filter out
some data points.

| Name              | Description                                                                                                                 |
| ----------------- | --------------------------------------------------------------------------------------------------------------------------- |
| `path`            | Path to the directory containing the data files.                                                                            |
| `error`           | The method for error estimation.                                                                                            |
| `filter_offsets`  | The list of offsets to exclude from the calculation (in CEST experiments).                                                  |
| `filter_planes`   | The list of planes from which data points should be excluded from the calculation (in CEST and CPMG experiments).           |
| `[data.profiles]` | The subsection containing the list of the experimental profile file names are given along with their spin-system assignment |

#### `error`

The `error` key is a string providing the method for error estimation.

| Key value      | Description                                                                                                                                                                                                                      |
| -------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `"file"`       | Uncertainty values directly taken from the file, no estimation are run.                                                                                                                                                          |
| `"duplicates"` | Estimate the uncertainty using the [pooled standard deviation](https://goldbook.iupac.org/html/P/P04758.html). If no duplicates are found, the average of the experimental error values is returned.                             |
| `"scatter"`    | Estimate the uncertainty from CEST profiles. Assume the profile has additive, Gaussian noise on top of a smoothly varying signal. Adapted from [there](https://fr.mathworks.com/matlabcentral/fileexchange/16683-estimatenoise). |

#### `filter_offsets`

This is to filter out certain offsets from the CEST profiles. This can be
useful, for example, when decoupling sideband artifacts are present in the
profiles, as their locations are in general predictable. The filter is applied
to all the profiles of the experiments.

You should provide a list of pairs of values, where the first value of the pair
corresponds to the offset relative to the main resonance position (ν) and the
second value to the bandwidth around the offset (Δν) where points are excluded
from the calculation (ν ± 0.5 \* Δν). All values should be given in Hz. This key
is optional.

```toml
filter_offsets = [
   [0.0, 20.0],
   [-300.0, 20.0],
   [+300.0, 20.0],
 ]
```

#### `[data.profiles]`

This subsection of the `[data]` section contains the list of the filenames of
the experimental profiles and their spin-system assignments.

The name of the spin-systems follows the Sparky-NMR peak label conventions. Each
nucleus is designated by a **group name**, usually a 1- or 3-letter amino acid
with the position number in the protein sequence (e.g. ALA3 or A3), followed by
an **atom name** (e.g. N, CA, CG1, etc.).

For spin systems with more than one spin, use the "-" sign to separate the
different spin names. If successive spins of the system belong to the same
residue, the group ID (e.g. A12) can be omitted after the first spin. For
example, G23N-G23H and G23N-H are both valid and equivalent.

Here is an example:

```toml
  [data.profiles]
  13N = "13N-HN.out"
  26N = "26N-HN.out"
  28N = "28N-HN.out"
  29N = "29N-HN.out"
```

:::note

Choose the name of the spin-systems to reflect the spin-system of interest in
each experiment. For example, "G2N" is suitable for experiments based on a
single-spin system, whereas "G2N-H" is suitable for experiments based on a
two-spin system, etc. Although ChemEx has a built-in mechanism to correct
user-supplied spin system names automatically, we recommend that you be as
accurate as possible when choosing them.

:::

---
sidebar_position: 3
---

# Experiment Files

## Overview

Experiment files in ChemEx are used to provide detailed information about datasets and file locations for analysis. These files are accessible through the `-e` or `--experiments` command line option as shown below:

```shell
chemex fit -e <experiment_file1> <experiment_file2> [...]
```

## Structure of Experiment Files

Experiment files consist of three primary sections:

1. **`[experiment]`**: Contains experiment-related information.
2. **`[conditions]`**: Details sample conditions.
3. **`[data]`**: Includes data location, spin system assignments, and methods for uncertainty estimation or filtering measurements.

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
Each experiment in ChemEx comes with a sample configuration file. These files, accessible in the [Experiments](/docs/experiments) section, include detailed key-value pair descriptions.
:::

## Detailed Sections

### `[experiment]`

This section includes pulse sequence type and settings. Common keys found here are:

- `name`: Name of the pulse sequence.
- `carrier`: Carrier position in ppm.
- `time_t1`, `time_t2`: Relaxation delays in seconds.
- `pw90`: 90 degree pulse width in seconds.
- `b1_frq`: B1 field strength in Hz.
- `observed_state`: Observed state ID (a, b, c, d).

### `[conditions]`

Provides experimental and sample conditions like Larmor frequency, sample temperature, protein concentration, and labeling. While `h_larmor_frq` (Larmor frequency) is always required, `temperature`, `p_total`, and `l_total` depend on the kinetic model used.

#### `label`
The `label` key, specifying the sample's labeling scheme, is crucial for certain experiments. Examples include `"2H"` for deuterated samples and `"13C"` for <sup>13</sup>C-labeled samples. 

### `[data]`

Details data location, spin system assignments, and methods for error estimation and data filtering.

#### `error`
Methods for error estimation include:
- `"file"`: Direct from file.
- `"duplicates"`: Pooled standard deviation or average error values.
- `"scatter"`: Gaussian noise estimation on CEST profiles.

#### `filter_offsets`
Used to exclude specific offsets from CEST profiles, typically provided as pairs of values indicating offset and bandwidth.

#### `[data.profiles]`
Lists filenames of experimental profiles and their spin-system assignments. Names should follow Sparky-NMR peak label conventions.

:::note
Ensure accuracy in naming spin-systems, as it impacts experiment specificity.
:::