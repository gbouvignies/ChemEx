# Contributing to ChemEx

ChemEx is production scientific software. Contributions should remain focused,
preserve existing numerical and user-facing behaviour unless a change is
intentional, and include validation proportionate to their scientific impact.
Read [`AGENTS.md`](AGENTS.md) for the repository map, extension paths, numerical
invariants, and detailed validation matrix.

## Development setup

ChemEx requires Python 3.13 or newer. Python 3.13 is the blocking CI target;
Python 3.14 is currently tested as an informational target.

Install [uv](https://docs.astral.sh/uv/), clone your fork, and from the repository
root run:

```sh
uv sync --all-extras --dev
```

Create a branch and keep each pull request limited to one coherent change.

## Validation

Run focused tests while developing. Before submitting a Python change, run the
applicable repository checks:

```sh
uv run pytest -q
uv run ty check
uv run ruff check .
uv run ruff format --check <changed-python-files>
git diff --check
```

Use `Array` from `chemex.typing` for NumPy array annotations. Add precise type
annotations to new or modified interfaces where they clarify the contract; do
not weaken existing checks to make a change pass.

For website changes, run:

```sh
cd website
npm ci
npm run build
```

For packaging changes, also run `uv build`.

## Tests and scientific changes

- Bug fixes should include a regression test whenever practical.
- TOML or schema changes should cover valid and invalid inputs and should verify
  relevant files under `examples/`.
- Experiment and pulse-sequence changes should use realistic inputs and cover
  reference or limiting cases.
- Kinetic-model and parameter changes should verify rate/population constraints
  and relevant physical limits.
- Numerical changes should establish reference behaviour before implementation
  and use justified floating-point tolerances.
- CLI, output, or provenance changes should include the related integration tests
  and a concrete CLI run when practical.

Do not update expected values simply to accept unexplained numerical drift.
Document intentional changes to scientific behaviour, configuration, output, or
compatibility in the pull request and the user documentation or changelog when
appropriate.

## Reporting issues and proposing changes

Use the GitHub issue tracker for bugs and enhancements. Include a minimal
reproducer, input files or logs where shareable, the ChemEx and Python versions,
and the expected and actual behaviour. For a consequential scientific or
architectural proposal, describe alternatives, compatibility implications, and
a migration path before implementation.

## Pull requests and review

In the pull request description, report:

- the behaviour changed;
- compatibility implications;
- tests and checks run, with outcomes;
- numerical validation performed when relevant;
- documentation or examples updated;
- remaining risks or follow-up.

Maintainers will review the change and may request revisions before merge. Be
respectful and keep discussion focused on making ChemEx reliable for its users.

Questions are welcome in the repository's GitHub Discussions.
