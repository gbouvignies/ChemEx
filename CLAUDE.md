# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ChemEx is a Python-based NMR chemical exchange analysis tool designed for analyzing NMR experimental data to characterize chemical exchange processes. It specializes in CPMG relaxation dispersion and Chemical Exchange Saturation Transfer (CEST) experiments.

## Development Commands

### Installation and Setup
```bash
# Install in development mode with uv
pip install -e .

# Or with pip from source
pip install git+https://github.com/gbouvignies/ChemEx.git
```

### Running ChemEx
```bash
# Run a fit
chemex fit -e path/to/experiments.toml -p path/to/parameters.toml [-m path/to/method.toml] [-o Output]

# Run a simulation
chemex simulate -e path/to/experiments.toml -p path/to/parameters.toml [-o OutputSim]

# Pick CEST profiles
chemex pick_cest -e path/to/experiments.toml [-o Sandbox]

# Plot fitted parameters
chemex plot_param -p path/to/parameters.fit -n PARAMETER_NAME
```

### Linting and Code Quality
```bash
# Run ruff linter (configured in pyproject.toml)
ruff check .

# Auto-fix issues
ruff check --fix .

# Format code
ruff format .
```

### Type Checking
```bash
# Run pyright (configured in pyproject.toml)
pyright
```

## High-Level Architecture

### Core Data Flow
1. **Entry Point**: `src/chemex/chemex.py::main()` registers experiments and kinetic models, then parses CLI arguments
2. **Experiment Loading**: `src/chemex/experiments/builder.py::build_experiments()` reads TOML files and creates experiment objects
3. **Parameter Management**: `src/chemex/parameters/database.py` maintains a global parameter database
4. **Fitting/Simulation**: `src/chemex/optimize/fitting.py` orchestrates the optimization process
5. **Output**: Results are written to output directories with fitted parameters, plots, and statistics

### Key Components

#### Experiment System
- **Catalog**: `src/chemex/experiments/catalog/` contains all NMR experiment implementations (CPMG, CEST, DCEST, etc.)
- **Factories**: `src/chemex/experiments/factories.py` provides factory pattern for creating experiment-specific components
- **Builder**: `src/chemex/experiments/builder.py` constructs experiments from TOML configuration files
- **Loader**: `src/chemex/experiments/loader.py` registers all available experiment types via `register_experiments()`

Each experiment type defines:
- Configuration schema (Pydantic models)
- Dataset creation from data files
- Spectrometer/sequence setup
- Printer and plotter for output

#### Kinetic Models
- **Models**: `src/chemex/models/kinetic/` contains kinetic exchange models (2-state, 3-state, 4-state, binding models, Eyring, etc.)
- **Model Selection**: `src/chemex/models/model.py::model` is a singleton that tracks the current model (e.g., "2st", "3st_binding")
- **Factory**: `src/chemex/models/factory.py` creates model instances from names
- **Extensions**: Models can have extensions like `.mf` (model-free) or `.tc` (temperature coefficient)

#### Liouvillian Calculations
- **Core Module**: `src/chemex/nmr/liouvillian.py` handles quantum mechanical calculations
- **Basis**: `src/chemex/nmr/basis.py` defines operator basis sets for spin systems
- Supports various spin operators (Ix, Iy, Iz, etc.) and multi-spin terms (2IxSz, 2IzSy, etc.)

#### Configuration System
- **Base**: `src/chemex/configuration/base.py` provides Pydantic base models
- **Experiment**: `src/chemex/configuration/experiment.py` defines experiment configuration structure
- **Conditions**: `src/chemex/configuration/conditions.py` handles experimental conditions (temperature, field strength, labels, etc.)
- **Methods**: `src/chemex/configuration/methods.py` defines fitting methods and constraints
- All configurations use TOML format with Pydantic validation

#### Parameters
- **Database**: `src/chemex/parameters/database.py` is a global singleton managing all fitting parameters
- **Spin Systems**: `src/chemex/parameters/spin_system/` defines nuclei and spin system naming (e.g., "42N-HN")
- **Factory**: `src/chemex/parameters/factory.py` creates parameter objects from configurations
- Parameters include chemical shifts, exchange rates, relaxation rates, etc.

#### Optimization
- **Fitting**: `src/chemex/optimize/fitting.py` manages the fitting workflow (steps, groups, methods)
- **Minimizer**: `src/chemex/optimize/minimizer.py` wraps lmfit for parameter optimization
- **Gridding**: `src/chemex/optimize/gridding.py` performs grid searches over parameter space
- **Statistics**: Supports Monte Carlo, bootstrap, and nucleus-based bootstrap analysis

### Data Containers
- **Profile**: `src/chemex/containers/profile.py` represents data for a single spin system in an experiment
- **Experiment**: `src/chemex/containers/experiment.py` contains all profiles for a single experiment file
- **Experiments**: `src/chemex/containers/experiments.py` manages multiple experiments across files
- **Dataset**: `src/chemex/containers/dataset.py` is a list of (SpinSystem, data) tuples

### Input/Output
- **TOML Reader**: `src/chemex/toml.py` handles TOML file parsing
- **Messages**: `src/chemex/messages.py` centralizes all console output formatting (uses Rich library)
- **Printers**: `src/chemex/printers/` formats experimental data for output
- **Plotters**: `src/chemex/plotters/` creates matplotlib visualizations

## Code Style and Conventions

- **Linting**: Uses Ruff with strict settings (select = ["ALL"]) except specific ignored rules in pyproject.toml
- **Code Style**: Follows black-compatible formatting
- **Type Hints**: Uses type hints throughout; checked with Pyright
- **Python Version**: Requires Python >=3.12
- **Docstrings**: Some modules have docstrings, but many functions/classes lack them (ignored by Ruff: D100-D105, D107)

## Important Notes

- **Global State**: The `model` singleton in `src/chemex/models/model.py` and `database` in `src/chemex/parameters/database.py` maintain global state
- **Registration Pattern**: Experiments and kinetic models must be registered via loader modules before use
- **Factory Pattern**: Heavy use of factory pattern for creating experiment-specific components
- **Rich Output**: All console output uses Rich library for formatting and progress bars
- **TOML Configuration**: All inputs (experiments, parameters, methods) use TOML format with Pydantic validation

## Examples

The `examples/` directory contains:
- **Experiments/**: Individual experiment type examples (CPMG_15N_IP, CEST_15N, DCEST_15N, etc.)
- **Combinations/**: Multi-experiment fitting examples (N15_NH_RDC, 2stBinding, etc.)

Each example includes:
- Experiment TOML files (define experimental setup and data paths)
- Parameters TOML files (initial parameter values)
- Methods TOML files (fitting strategy, constraints, statistics)
- Data files (raw NMR measurements)

## Documentation

Full documentation available at: https://gbouvignies.github.io/ChemEx/
