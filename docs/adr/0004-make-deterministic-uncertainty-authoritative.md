---
status: accepted
---

# Make Deterministic Uncertainty authoritative

ChemEx will represent the complete scientific interpretation of uncertainty for
one accepted deterministic fit as **Deterministic Uncertainty**.  It contains
the validated root and block Evidence it interpreted, orthogonal interpretation
completeness and derivation-disposition facts, and predetermined per-parameter
reportability, unavailability classification, and boundary-warning facts.

`chemex.optimize.deterministic_uncertainty` owns the translation from narrow
authoritative accepted-fit facts into the existing production uncertainty
policy, constrained-output scope and capabilities, low-level Evidence
derivation, root-anchored block recovery, recovery precedence, completeness,
withholding, and parameter conclusions.  Its sole production interpretation
entry point is:

`AcceptedDeterministicFitFacts` → `derive_deterministic_uncertainty(...)` →
`DeterministicUncertainty`.

Direct TRF, grouped Direct TRF, and DE-seeded final Direct TRF use one continuous
uncertainty basis carrying the existing fit partition proof.  Profiled GRID has
a distinct payload-free basis and produces a complete, policy-withheld
Deterministic Uncertainty without Evidence.  An evaluation-only occurrence with
no accepted deterministic fit has no Deterministic Uncertainty.

Low-level Evidence types, covariance numerics, constraint linearization, and
root-anchored block derivation remain in `chemex.optimize.uncertainty`.
Deterministic fitting initiates interpretation after commit and then constructs
`NativeDeterministicFit`, which carries one required Deterministic Uncertainty
as its sole downstream uncertainty authority.

Publication consumes the authoritative outcome and owns only representation,
formatting, artifact lifecycle, and historical compatibility mappings.  In
particular, the internal scientific GRID state `complete + withheld` continues
to serialize as the established covariance status `incomplete` with terminal
`withheld` and the existing reason text.  This deliberate translation does not
change the scientific domain state.

## Considered options

- Keeping root Evidence, block Evidence, a printer view, and a status tuple as
  parallel downstream inputs was rejected because it leaves scientific
  authority distributed across fitting and publication.
- A resolver/interpreter object with an Evidence-injection path was rejected
  because low-level Evidence already has an independent qualification seam and
  production needs one interpretation entry point.
- Passing grouped or GRID execution aggregates was rejected because the
  uncertainty module needs only narrow accepted-fit lineage and basis facts.
- A generic `UncertaintyResult`, shared statistical-result hierarchy, combined
  status enum, or new GRID unavailability kind was rejected because the
  scientific dimensions and policies are distinct.
- New witnesses, content hashes, serialization contracts, factories, ports,
  compatibility aliases, feature flags, and dual implementations were rejected
  because existing accepted-fit and Evidence identities already establish the
  required lineage.
- Converting `uncertainty.py` into a package was rejected as an unrelated
  numerical-module migration.

## Consequences

Every accepted deterministic fit has an immutable Deterministic Uncertainty,
including complete evaluated, partially reportable, incomplete interrupted, and
complete withheld outcomes.  Root interruption retains no Evidence; block
interruption retains valid root Evidence and conclusions.  Closed scientific
Evidence failures such as rank deficiency affect parameter reportability but do
not make interpretation incomplete.

The old root/block fields, printer-side scientific classification, private
product-policy helpers, and parallel publication arguments are removed in one
cutover.  No compatibility aliases remain.  Supported CLI and Method TOML
behavior, central values, policies, numerical tolerances, exception behavior,
progress text, artifact schemas, and rendered uncertainty strings remain
unchanged.  Pre-refactor numerical qualifications and fixed representative
artifact data form the compatibility oracle; any later scientific-policy change
requires a separate decision.
