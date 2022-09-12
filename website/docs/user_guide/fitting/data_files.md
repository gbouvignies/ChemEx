---
sidebar_position: 4
---

# Data files

## Description

The location of data files is specified within
[experiment files](experiment_files.md). Data files typically have three columns
containing the following information:

| Experiment type          | First column     | Second column | Third column |
| ------------------------ | ---------------- | ------------- | ------------ |
| CEST / D-CEST / cos-CEST | Offset (Hz)      | Intensity     | Uncertainty  |
| CPMG                     | Number of cycles | Intensity     | Uncertainty  |
| Relaxation               | Delay (s)        | Intensity     | Uncertainty  |

Data files are peak intensity tables that can be obtained from the output of
many peak fitting programs, such as the `autoFit` subroutine in the
[NMRPipe](https://www.ibbr.umd.edu/nmrpipe/index.html) program.

## Example

import Tabs from '@theme/Tabs'; import TabItem from '@theme/TabItem';

<Tabs>
<TabItem value="cest" label="CEST" default>

```python title="./Data/30Hz/F52N-H.out"
# Offset (Hz)        Intensity    Uncertainty
   -2.000e+04    2.5352430e+07  7.5500000e+04
   -4.000e+02    1.7378230e+07  7.5500000e+04
   -3.500e+02    1.7043830e+07  7.5500000e+04
   -3.000e+02    1.5819080e+07  7.5500000e+04
   -2.500e+02    8.8238340e+06  7.5500000e+04
   -2.000e+02    3.0651600e+06  7.5500000e+04
   -1.500e+02    1.4991090e+07  7.5500000e+04
   -1.000e+02    1.6768950e+07  7.5500000e+04
   -5.000e+01    1.7215280e+07  7.5500000e+04
    0.000e+00    1.7636170e+07  7.5500000e+04
    5.000e+01    1.7612980e+07  7.5500000e+04
    1.000e+02    1.7666720e+07  7.5500000e+04
    1.500e+02    1.7694720e+07  7.5500000e+04
    2.000e+02    1.7692340e+07  7.5500000e+04
    2.500e+02    1.7747680e+07  7.5500000e+04
    3.000e+02    1.7633460e+07  7.5500000e+04
    3.500e+02    1.7461530e+07  7.5500000e+04
    4.000e+02    1.7377560e+07  7.5500000e+04
```

</TabItem>
<TabItem value="cpmg" label="CPMG">

```python title="./Data/800MHz/G39N-H.out"
#    ncyc_cp        Intensity      Esd(Int.)
   0.000e+00    1.2584470e+07  9.8068130e+03
   4.000e+01    9.1841462e+06  9.8068130e+03
   1.000e+00    6.7679280e+06  9.8068130e+03
   3.600e+01    9.0545262e+06  9.8068130e+03
   2.000e+00    6.7616357e+06  9.8068130e+03
   3.200e+01    8.8758267e+06  9.8068130e+03
   3.000e+00    6.7377252e+06  9.8068130e+03
   3.000e+01    8.7499820e+06  9.8068130e+03
   4.000e+00    6.8950311e+06  9.8068130e+03
   2.800e+01    8.6316880e+06  9.8068130e+03
   5.000e+00    6.8094567e+06  9.8068130e+03
   2.600e+01    8.5410798e+06  9.8068130e+03
   6.000e+00    6.6798367e+06  9.8068130e+03
   2.400e+01    8.4227858e+06  9.8068130e+03
   7.000e+00    6.6458586e+06  9.8068130e+03
   2.200e+01    8.2566708e+06  9.8068130e+03
   8.000e+00    6.7125563e+06  9.8068130e+03
   2.000e+01    8.0012060e+06  9.8068130e+03
   9.000e+00    6.6949380e+06  9.8068130e+03
   1.800e+01    7.8099221e+06  9.8068130e+03
   1.000e+01    6.7440175e+06  9.8068130e+03
   1.600e+01    7.5330637e+06  9.8068130e+03
   1.100e+01    6.8925142e+06  9.8068130e+03
   1.400e+01    7.2436209e+06  9.8068130e+03
   1.200e+01    6.9504028e+06  9.8068130e+03
   1.000e+00    6.8107152e+06  9.8068130e+03
   4.000e+01    9.1904384e+06  9.8068130e+03
   2.000e+00    6.7855462e+06  9.8068130e+03
   3.600e+01    9.0104805e+06  9.8068130e+03
```

</TabItem>
<TabItem value="relaxation" label="Relaxation">

```python title="./Data/800MHz/S5N-H.out"
#  Times (s)        Intensity      Esd(Int.)
   5.000e-02    9.2822624e+05  1.0000000e+04
   1.000e-01    8.6879557e+05  1.0000000e+04
   1.500e-01    7.8922229e+05  1.0000000e+04
   2.000e-01    7.7218645e+05  1.0000000e+04
   2.500e-01    6.8675373e+05  1.0000000e+04
   3.000e-01    6.4862935e+05  1.0000000e+04
   3.500e-01    6.1765521e+05  1.0000000e+04
   4.000e-01    5.5875618e+05  1.0000000e+04
   4.500e-01    5.2201228e+05  1.0000000e+04
   5.000e-01    4.8749523e+05  1.0000000e+04
```

</TabItem>
</Tabs>

:::note

The input data files for `shift_experiments` have special format, refer to the
`Shifts/` example under `examples/Combinations/` to learn how to create data
files for `shift_experiments`.

:::
