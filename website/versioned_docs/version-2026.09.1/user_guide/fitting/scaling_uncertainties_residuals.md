---
sidebar_position: 6.6
title: Scaling, Uncertainties, Residuals, and Fit Statistics
description: Follow ChemEx data from profile scaling and uncertainty estimation to normalized residuals and fit statistics.
---

# Scaling, uncertainties, residuals, and fit statistics

ChemEx fits a residual vector built from selected profiles and their retained
observations. Profile scaling, uncertainty estimation, and masks determine the
weight of every contribution, so they are part of the scientific model rather
than presentation-only settings.

## From experimental points to the fit objective

For each selected profile, ChemEx:

1. obtains point uncertainties from the data file or the configured estimator;
2. calculates the profile for the current parameter values;
3. identifies the observations retained by the active mask;
4. analytically determines an amplitude scale from those observations when
   `scaled = true`; and
5. concatenates the normalized residuals from all retained observations.

The optimizer minimizes the sum of squares of that residual vector.

## Profile scaling

`scaled = true` is the common default. Shift experiments override the default to
`false`. Check the reference page for an experiment before deciding whether an
amplitude scale is scientifically appropriate.

For every scaled profile, ChemEx independently recomputes the weighted
least-squares scale during evaluation. Using only retained observations, the
scale is

$$
s =
\frac{\sum_i y_i f_i / \sigma_i^2}
{\sum_i f_i^2 / \sigma_i^2},
$$

where $y_i$ is the experimental value, $f_i$ is the unscaled calculated value,
and $\sigma_i$ is the experimental uncertainty.

The scaled calculation is $y_i^{\mathrm{calc}} = s f_i$. This lets an overall
profile amplitude adjust analytically, so the fit is driven primarily by the
profile shape. With `scaled = false`, ChemEx uses the unscaled calculation and
retains absolute-amplitude information in the objective.

## Experimental uncertainty modes

The `[data]` table in an experiment file selects `error = "file"`,
`"scatter"`, or `"duplicates"`.

### `file`

ChemEx uses the pointwise uncertainty supplied in the third column of each data
file. The points in one profile may therefore have different weights.

### `scatter`

ChemEx estimates one noise variance per profile from local point-to-point
scatter using a finite-difference estimator. The resulting scalar standard
uncertainty replaces the file uncertainties for every point in that profile,
unless experiment-level pooling is enabled as described below. This estimator
is intended for suitable sampled profiles such as CEST data; it is not a
universal uncertainty model.

### `duplicates`

Within a profile that contains repeated metadata values, ChemEx calculates the
sample variance of each duplicate group and combines those variances with their
within-group degrees of freedom. The resulting profile variance gives one
standard uncertainty for that profile unless experiment-level pooling is
enabled.

Duplicate availability is handled as follows:

- If every profile has duplicates, every profile contributes its
  duplicate-derived variance.
- If no profile in the constructed experiment has duplicates, ChemEx preserves
  the original pointwise file uncertainties and emits a fallback notice. This
  happens regardless of `global_error`.
- If duplicate and duplicate-free profiles are mixed, a duplicate-free profile
  contributes the mean of its pointwise file variances:

$$
v_{\mathrm{fallback}} = \frac{1}{n}\sum_i \sigma_i^2.
$$

With profile-specific estimation, the scalar uncertainty assigned to that
duplicate-free profile is therefore the root-mean-square file uncertainty:

$$
\sigma_{\mathrm{fallback}} =
\sqrt{\frac{1}{n}\sum_i \sigma_i^2}.
$$

For example, file uncertainties of 0.03 and 0.05 give
$\sqrt{(0.03^2 + 0.05^2)/2} \approx 0.0412$.

## Global versus profile-specific error estimation

`global_error = true` is the default. It affects the `scatter` and `duplicates`
estimated-error modes, but not `file` mode.

- With `global_error = false`, each profile keeps its own estimated scalar
  uncertainty.
- With `global_error = true`, ChemEx takes the equal-profile mean of the
  per-profile variances, takes its square root, and assigns that common scalar
  uncertainty to every profile in the constructed experiment.

For `duplicates`, this means each profile contributes one variance to the global
pool. A profile with duplicates contributes its pooled duplicate variance; in a
mixed experiment, a duplicate-free profile contributes its mean file variance.

Here, **global means all profiles built from one experiment TOML file**. ChemEx
does not pool error estimates across separate experiment files in the complete
fit.

## Selection, filtering, and active observations

Only retained active observations contribute to the residual vector and χ².
There are two distinct levels of selection:

- data-file and Method-file selection determine which profiles are constructed
  or active in a fit step;
- experiment-specific filters set individual observations inactive within an
  active profile.

Common filters can exclude plane indices or offset regions. In CEST-style
configuration, `filter_ref_planes = true` excludes reference planes from the
objective; the common CEST default is `false`, while experiment types without
usable reference planes can override it. A point is not excluded merely because
it is marked as a reference plane.

Masked observations are omitted from scaling, the residual vector, and χ².
Calculated values can remain available in machine-readable `Data/` output when
the relevant mechanism retains them. Plot visibility depends on the plotter and
filter path: ordinary non-reference filtered observations may appear as lighter
excluded points where supported, whereas CEST reference planes removed by
`filter_ref_planes` may be omitted from the plotted profile. See the relevant
[experiment reference](../../experiments/index.mdx) for supported filter keys
and their physical meaning.

## Normalized residuals

For every retained observation, the native fitting residual is

$$
r_i =
\frac{y_i^{\mathrm{calc}} - y_i^{\mathrm{exp}}}{\sigma_i}.
$$

The sign is calculated minus experimental. For scaled profiles,
$y_i^{\mathrm{calc}}$ already includes the analytically optimized profile scale.
Masking is applied before residuals are concatenated.

Residuals are expressed in units of the supplied or estimated experimental
uncertainty. Thus $|r_i| \approx 1$ means a discrepancy of about one standard
uncertainty **if the uncertainty model itself is meaningful**. Normalizing a
residual does not make it statistically normal or remove systematic structure.

## χ² and residual degrees of freedom

ChemEx uses these counts:

- $N$: retained residual observations;
- $P$: optimizer-controlled fit coordinates;
- $G$: analytically fitted profile normalizations, one for each scaled active
  profile;
- $K = P + G$: total fitted model dimension used for AIC and BIC;
- $\nu = N - P - G$: residual degrees of freedom.

The objective is

$$
\chi^2 = \sum_i r_i^2.
$$

When $\nu > 0$, reduced χ² is

$$
\chi^2_\nu = \frac{\chi^2}{\nu}.
$$

The χ² goodness-of-fit calculation uses the same $\nu$. When $\nu \le 0$, the
fit can still complete: χ², AIC, and BIC remain defined, while reduced χ² and the
χ² goodness-of-fit result are reported as `nan`.

:::note Artifact naming

`"number of variables"` in `statistics.toml` reports $P$, the
optimizer-controlled coordinates. It does not include the analytically profiled
normalization count $G$. Nevertheless, $G$ contributes to residual degrees of
freedom, the χ² goodness-of-fit calculation, AIC, and BIC. No separate $G$ field
is currently serialized.

:::

For example, the first CPMG tutorial has $N = 230$, $P = 17$, and $G = 10$, so
$\nu = 203$. See [Run Your First ChemEx Analysis](../../first_analysis.md) and
the [full CPMG workflow](../../full_cpmg_analysis.md).

## AIC and BIC

Within ChemEx's weighted least-squares convention, the current statistics use
$K = P + G$:

$$
AIC = \chi^2 + 2K
$$

and

$$
BIC = \chi^2 + K \ln N.
$$

This is why the model dimension used by AIC/BIC can exceed the serialized
`"number of variables"`. See [Outputs](outputs.mdx#statisticstoml) for the
artifact format.

## Interpretation cautions

A reduced χ² near 1 is consistent with residual scatter on the scale of the
stated uncertainties, under the fitted model. It is not proof that the kinetic
model is correct, parameters are identifiable, residuals are structure-free, or
the uncertainty estimates are trustworthy.

Likewise, optimizer success means that the requested numerical procedure
completed under its convergence rules. It does not establish scientific model
adequacy. Inspect residual patterns, parameter constraints and boundaries,
covariance evidence, and sensitivity to scientifically plausible alternatives
before drawing conclusions.
