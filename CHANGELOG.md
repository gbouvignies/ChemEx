# Changelog

All notable changes to ChemEx will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project uses [Calendar Versioning](https://calver.org/) (YYYY.MM.MICRO).

## [2025.10.0] - 2025-10-27

### Added
- **New Experiments**: Added fitting module for 13CH3 13C CPMG [0013] experiment (#374)
- **New Features**:
  - Implemented 4-state Eyring model for chemical exchange kinetics (#384)
  - Added `detect_all_states` option for CEST, DCEST, and CPMG experiments (#385)
- **Documentation**: Added CLAUDE.md for project guidance and development instructions

### Changed
- **Build System**:
  - Migrated to `src/` layout for better package structure
  - Updated build system configuration to use `uv_build`
  - Switched to static versioning (removed hatch-vcs dependency)
- **Type System**: Major refactor of array typing across the codebase
  - Updated type annotations to use `Array` for intensities and magnetization
  - Improved type consistency and clarity across parameters
- **Configuration**: Refactored configuration merging to use `deep_update` for improved readability
- **Documentation**:
  - Updated installation instructions for uv in README and welcome guide
  - Updated Docusaurus config to enable new experimental features

### Fixed
- Fixed initial conditions for two-spin systems in ChemEx (#392)
- Updated FIX parameters in `method.toml` for CPMG_CH3_MQ to fix DW_AB configurations
- Updated suffix handling in CpmgCh313CH2c0013Settings for start_terms and detection properties

### Dependencies
- **Python**: Bumped minimum version to 3.13
- **Core Dependencies**:
  - matplotlib: 3.10.0 → 3.10.1
  - numpy: 2.2.3 → 2.2.5
  - pydantic: 2.10.6 → 2.11.4
  - rapidfuzz: 3.12.1 → 3.14.1
  - lmfit: 1.3.2 → 1.3.3
  - cachetools: 5.5.1 → 5.5.2
  - rich: 13.9.4 → 14.0.0

### Infrastructure
- **CI/CD**:
  - Refactored workflows and fixed dependabot group conflicts
  - Updated GitHub Actions: setup-python (5→6), setup-node (4→6), checkout (4→5)
  - Updated sigstore action (3.0.0→3.1.0)
  - Added Claude Code GitHub Workflow (#383)
  - Updated Node.js version to 20 in deployment workflows
- **Dependabot**: Updated configuration for pip, uv, github-actions, and npm ecosystems

### Website
- Multiple dependency updates for Docusaurus and related packages
- Security updates for various npm packages

## [2025.4.0] - 2025-04-XX

### Added
- Support for 5-state and 6-state kinetic models
- Enhanced CEST plotter with updated line styles

---

For a complete list of changes, see the [GitHub Releases](https://github.com/gbouvignies/chemex/releases) page.
