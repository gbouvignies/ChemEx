---
status: accepted
---

# Make specialized statistical results authoritative

ChemEx will distinguish **Evidence**, **Summary**, and a specialized **Analysis
Result** for resampling and MCMC. Evidence contains validated observations or
samples and their scientific state. A Summary contains scientific conclusions
derived from that Evidence under an explicit interpretation policy. A
`ResamplingAnalysisResult` or `McmcAnalysisResult` is the complete authoritative
scientific outcome consumed downstream, including its Evidence, Summary
availability, and scientific completeness.

The specialized native modules own the translation from canonical resolved
requests into immutable interpretation policy, evidence validation, retention
or burn semantics, completeness classification, summary derivation, covariance
or correlation availability, exclusions, and construction of the specialized
Analysis Result. This translation does not reinterpret facts already made
authoritative at the Method Plan and execution boundary. Existing production
and qualification policies may differ where they already do; those policies
remain explicit within the specialized module rather than being silently
unified.

Publication is a non-scientific consumer of an authoritative Analysis Result.
It owns representation, rendering, and artifact lifecycle, but does not inspect
sample arrays to derive statistics or make scientific decisions. For a
scientifically meaningful incomplete outcome with validated Evidence,
publication reproduces the existing diagnostic artifact set before outer
orchestration raises the existing compatibility exception. Operational
interruption, unexpected execution failure, filesystem failure, and publication
failure continue to propagate directly. ChemEx does not manufacture an Analysis
Result when insufficient validated Evidence exists to define a scientific
outcome.

## Considered options

- A shared `StatisticalAnalysisResult` base class or generalized statistics
  framework was rejected because resampling and MCMC do not share enough
  interpretation semantics. Burn, retention, and chain completeness remain
  MCMC-specific.
- Runtime dual implementations, feature flags, and lasting compatibility paths
  were rejected because they would preserve duplicate scientific authority.
- Immediate unification of production and qualification policies was rejected
  because this decision changes authority, not scientific policy.

## Consequences

The migration proceeds sequentially: resampling is made authoritative and
validated before MCMC work begins. For each family, pre-refactor production
behavior is captured as the oracle, production switches to the specialized
Analysis Result, equivalence is established, and superseded calculations and
obsolete helper tests are deleted.

Existing schemas, quantiles, interpolation methods, MCMC burn behavior,
numerical results, run interfaces, compatibility exceptions, and observable
artifacts remain unchanged. Stable TOML and TSV artifacts require byte-for-byte
equivalence; deterministic seeded samples require exact array equality; timing
metadata is compared structurally. Any later scientific-policy unification is
a separate explicit decision.
