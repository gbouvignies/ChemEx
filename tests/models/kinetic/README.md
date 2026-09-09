# Eyring kinetic-model tests

This directory qualifies the shared Eyring scientific authority and all five
public Eyring model names. The focused tests cover the independent
high-precision numerical oracle, binary64 representability boundaries,
thermodynamic populations, detailed balance and four-state cycle closure,
topology integration, condition identity, deterministic uncertainty, MCMC, and
resampling.

## Running Tests

```bash
uv run pytest -q tests/models/kinetic -k eyring
```
