# Copilot Instructions for ChemEx

This guide enables AI coding agents to be productive in the ChemEx codebase. It summarizes architecture, workflows, and conventions specific to ChemEx.

## Project Overview
- **ChemEx** is a Python tool for analyzing NMR chemical exchange data, focusing on CPMG and CEST experiments.
- Main logic is in `src/chemex/` (core, CLI, models, experiments, optimization, plotting).
- Example scripts and data are in `examples/`.
- Website docs are in `website/` (Docusaurus).

## Architecture & Patterns
- **CLI Entrypoint:** `src/chemex/cli.py` and `src/chemex/__main__.py` implement command-line interface.
- **Modular Design:** Submodules for configuration, containers, experiments, models, nmr, optimize, parameters, plotters, printers, tools.
- **Configuration:** Uses TOML files for experiments, parameters, and methods. See `examples/` for templates.
- **Model Registration:** Models are registered in `src/chemex/models/` and used via CLI commands.
- **Testing:** Tests are in `tests/`, with model-specific tests in `tests/models/kinetic/`.

## Developer Workflows
- **Install (dev):**
  ```bash
  pip install -e .
  # or
  pip install git+https://github.com/gbouvignies/ChemEx.git
  ```
- **Run CLI:**
  ```bash
  chemex fit -e <experiments.toml> -p <parameters.toml> [-m <method.toml>] [-o Output]
  chemex simulate -e <experiments.toml> -p <parameters.toml> [-o OutputSim]
  chemex pick_cest -e <experiments.toml> [-o Sandbox]
  chemex plot_param -p <parameters.fit> -n PARAMETER_NAME
  ```
- **Lint:**
  ```bash
  ruff check .
  # Auto-fix: ruff check . --fix
  ```
- **Test:**
  ```bash
  python -m pytest tests/ -v
  # Coverage:
  python -m pytest tests/models/kinetic/ --cov=src/chemex/models/kinetic/settings_4st_eyring --cov-report=html
  ```
- **Website:**
  ```bash
  yarn start   # dev server
  yarn build   # static build
  yarn deploy  # deploy docs
  ```

## Conventions & Integration
- **Python 3.13 required** (see README for details).
- **Code style:** Black, Ruff (configured in `pyproject.toml`).
- **TOML config:** All experiment/parameter/method files use TOML format.
- **Examples:** Use scripts in `examples/Combinations/` and `examples/Experiments/` for typical workflows.
- **External dependencies:** See `pyproject.toml` for required packages.

## Key Files & Directories
- `src/chemex/cli.py`, `src/chemex/__main__.py`: CLI logic
- `src/chemex/models/`: Model implementations
- `examples/`: Example data and scripts
- `tests/`: Unit and integration tests
- `website/`: Documentation site

## Tips for AI Agents
- Reference TOML templates in `examples/` when generating configs.
- Follow CLI usage patterns from `CLAUDE.md` and `README.md`.
- Use model registration and test patterns from `tests/models/kinetic/`.
- Prefer Pythonic, modular code; keep CLI and core logic separate.

---
If any section is unclear or missing, please ask for feedback to improve these instructions.
