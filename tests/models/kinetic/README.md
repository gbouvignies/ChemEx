# Tests for ChemEx 4-State Eyring Model

This directory contains comprehensive tests for the `4st_eyring` kinetic model implementation.

## Test Structure

- `test_settings_4st_eyring.py`: Comprehensive unit tests for the core functionality
- `test_4st_eyring_integration.py`: Integration tests for model registration and framework integration

## Running Tests

### Prerequisites

Install pytest if not already available:
```bash
pip install pytest pytest-cov
```

### Run All Tests

```bash
# From the repository root
python -m pytest tests/models/kinetic/ -v

# Run with coverage
python -m pytest tests/models/kinetic/ --cov=src/chemex/models/kinetic/settings_4st_eyring --cov-report=html
```

### Run Specific Test Classes

```bash
# Unit tests only
python -m pytest tests/models/kinetic/test_settings_4st_eyring.py::TestCalculateKij4stEyring -v

# Integration tests only
python -m pytest tests/models/kinetic/test_4st_eyring_integration.py -v
```

## Test Coverage

The test suite covers:

### Core Function Tests (`TestCalculateKij4stEyring`)
- Basic functionality and return value validation
- Temperature range validation and boundary conditions
- Temperature dependence (Arrhenius behavior)
- Correct implementation of Eyring equation
- Detailed balance consistency for forward/reverse reactions
- Entropy effects on equilibrium ratios
- Numerical stability with extreme parameters
- Function caching behavior
- Parameter type handling (int, float, numpy types)

### Settings Creation Tests (`TestMakeSettings4stEyring`)
- Parameter creation and default values
- Temperature validation error handling
- Expression generation with temperature embedding
- Parameter bounds and variation flags

### Physical Realism Tests (`TestPhysicalRealism`)
- Arrhenius plot linearity across temperature range
- Detailed balance for all state pairs
- Correlation coefficient validation for temperature dependence

### Integration Tests (`TestEyring4stIntegration`)
- Model factory registration
- Multi-temperature consistency
- Parameter bounds and defaults validation
- Rate constant and population generation
- Framework integration testing

## Test Data

Tests use realistic thermodynamic parameters:
- State energies: 8-25 kJ/mol above ground state
- Activation energies: 50-120 kJ/mol (typical for protein folding)
- Temperature range: -100°C to +200°C for validation
- Entropy values: Typically set to 0 for simplicity

## Expected Behavior

### Passing Tests
All tests should pass without errors when:
- The Eyring equation is correctly implemented
- Temperature validation works properly
- Detailed balance is maintained
- Physical constants are used correctly
- Parameter generation is complete

### Common Failure Modes
Tests may fail if:
- Scipy constants are not available
- Temperature conversion is incorrect
- Activation free energy calculation has errors
- Rate constant expressions are malformed
- Model registration fails

## Performance Considerations

- Tests include caching validation to ensure performance
- Large parameter sweeps are avoided for speed
- Numerical precision is tested to 1e-10 relative tolerance
- Integration tests verify model factory efficiency

## Adding New Tests

When adding tests:
1. Follow pytest conventions with `test_` prefixes
2. Use descriptive test method names
3. Include docstrings explaining test purpose
4. Use realistic parameter values
5. Test both normal and edge cases
6. Validate physical meaning of results

## Troubleshooting

### Import Errors
If tests fail with import errors:
```bash
# Ensure the package is in PYTHONPATH
export PYTHONPATH="${PYTHONPATH}:/path/to/ChemEx/src"
```

### Missing Dependencies
Install required packages:
```bash
pip install scipy numpy pytest
```

### Numerical Precision Issues
Adjust tolerance in assertions if needed, but investigate first to ensure correctness.
