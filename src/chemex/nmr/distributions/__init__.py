"""B1 field inhomogeneity distribution plugins.

This package provides a plugin-based system for B1 field inhomogeneity
distributions used in NMR relaxation experiments. Each distribution is
self-contained in its own module and is automatically discovered and
registered at import time.

Available Distributions
-----------------------
The following distributions are currently available:

- **gaussian**: Original ChemEx distribution using a linear grid from -2sigma
  to +2sigma with Gaussian PDF weights (default, backward compatible)

- **hermite**: Gauss-Hermite quadrature distribution using Hermite polynomial
  roots for optimal integration of Gaussian-weighted functions

- **skewed**: Skew-normal distribution with configurable skewness parameter,
  useful for asymmetric B1 inhomogeneity profiles centered around nominal B1.
  Note: unbounded, can extend beyond nominal B1.

- **truncated_skewed**: Truncated skew-normal distribution with hard upper
  bound at nominal B1. Combines asymmetry with upper limit enforcement.
  Ideal for asymmetric profiles where B1 cannot exceed the RF power limit.

- **beta**: Upper-bounded distribution on [0, nominal_B1] for modeling B1
  degradation (coil inhomogeneity, sample loading). Use when the nominal B1
  is a physical upper limit and the field can only degrade below that value,
  never exceed it. Ideal for matching nutation experiment profiles.

- **dephasing**: Dephasing mode for extreme B1 inhomogeneity, replacing the
  old b1_inh_scale = inf hack

- **custom**: User-defined distribution with explicit B1 scaling factors and
  weights specified directly in the configuration

Distribution Selection Guide
-----------------------------
Choose the appropriate distribution based on your physical setup:

- **Symmetric inhomogeneity** (field varies equally above/below nominal):
  Use gaussian or hermite

- **Asymmetric inhomogeneity** (field has skewed profile but can exceed nominal):
  Use skewed with negative/positive skewness parameter

- **Asymmetric + upper-bounded** (skewed profile, cannot exceed nominal):
  Use truncated_skewed - combines asymmetry with hard upper limit

- **Upper-bounded inhomogeneity** (field cannot exceed nominal RF power):
  Use beta - common for:
  * Coil inhomogeneity (field weaker away from coil center)
  * Sample loading effects (RF penetration limits)
  * Cryoprobes with power-limited amplifiers
  * When nutation experiments show sharp upper cutoff

- **Complete dephasing** (extreme inhomogeneity, phase incoherence):
  Use dephasing

- **Known empirical profile** (from nutation or other measurement):
  Use custom with explicit B1 values and weights

Usage in TOML Configuration
----------------------------
B1 distributions can be specified in experiment TOML files using inline tables:

Example 1: Beta distribution (upper-bounded, default mean at 95% of nominal)
    [experiment]
    b1_frq = 20.0
    b1_distribution = { type = "beta", scale = 0.10, res = 11 }

Example 2: Beta with custom mean position (e.g., 90% of nominal)
    [experiment]
    b1_frq = 20.0
    b1_distribution = { type = "beta", scale = 0.15, res = 11, mean_frac = 0.90 }

Example 3: Skewed distribution (asymmetric but unbounded)
    [experiment]
    b1_frq = 15.95
    b1_distribution = { type = "skewed", scale = 0.15, skewness = -0.5, res = 11 }

Example 4: Custom user-defined distribution
    [experiment]
    b1_frq = 15.95
    b1_distribution = { type = "custom", scales = [0.95, 1.0, 1.03], weights = [0.25, 0.5, 0.25] }

Example 5: Default Gaussian (backward compatible)
    [experiment]
    b1_frq = 15.0
    b1_inh_scale = 0.1
    b1_inh_res = 11

Adding a New Distribution
--------------------------
To add a new distribution, create a new Python module in this package following
this template:

1. Create a file: src/chemex/nmr/distributions/my_distribution.py

2. Implement three components:
   - A `generate()` function that creates the distribution
   - A Pydantic config class that defines the TOML schema
   - A `register()` function that registers the distribution

3. Template::

    from __future__ import annotations
    from typing import Literal
    import numpy as np
    from pydantic import BaseModel, Field
    from chemex.nmr.constants import Distribution
    from chemex.nmr.distributions.registry import registry

    def generate(value: float, scale: float, res: int) -> Distribution:
        '''Generate my custom B1 distribution.'''
        values = np.linspace(value * (1 - scale), value * (1 + scale), res)
        weights = np.ones(res) / res
        return Distribution(values, weights)

    class MyDistributionConfig(BaseModel):
        '''Configuration for my distribution.'''
        type: Literal["my-distribution"] = "my-distribution"
        scale: float = Field(default=0.1, ge=0.0, le=1.0)
        res: int = Field(default=11, ge=1, le=51)

        def get_distribution(self, value: float) -> Distribution:
            return generate(value, self.scale, self.res)

    def register() -> None:
        registry.register("my-distribution", generate, MyDistributionConfig)

4. The distribution is automatically discovered and available for use

Key Requirements
----------------
- Use kebab-case for the distribution type name (e.g., "my-distribution")
- The config class must have a `type` field with a Literal type annotation
- The config class must implement `get_distribution(self, value: float)`
- The `generate()` function must return a Distribution namedtuple
- Weights must be normalized to sum to 1.0

Notes
-----
- All distributions are auto-loaded when this package is imported
- The registry is built dynamically from discovered plugins
- Adding a distribution requires no modifications to any central configuration
- Each distribution is fully self-contained and independently testable

"""

from __future__ import annotations

# Auto-load all distribution plugins
from chemex.nmr.distributions.loader import register_distributions

# Import registry and core functionality
from chemex.nmr.distributions.registry import (
    DistributionGenerator,
    get_b1_distribution,
    registry,
)

# Trigger auto-registration
register_distributions()

# Public API
__all__ = [
    "DistributionGenerator",
    "get_b1_distribution",
    "registry",
]
