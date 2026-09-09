# Chemical Exchange Analysis

ChemEx models NMR experiments that probe chemical exchange and compares their
calculated observables with measured data.

## Language

**Experiment Type**:
A named NMR acquisition and calculation protocol, such as CEST 15N or CPMG 15N
in-phase, that determines the scientific meaning of settings and observations.
_Avoid_: Experiment definition, experiment kind

**Experiment**:
One configured realization of an Experiment Type under particular conditions
and with particular measured data.
_Avoid_: Experiment Type

**Method Plan**:
An ordered, format-independent description of a fitting analysis, comprising
one or more Method Steps.
_Avoid_: Method file, Methods

**Method Step**:
One named stage in a Method Plan that defines profile selection, parameter
roles and constraints, search, and requested statistics for that stage.
_Avoid_: Method section, Method

**Generic N-state Model**:
A population-authoritative kinetic model whose public name declares a complete,
linear, or A-centered-fork structural exchange topology over three through six
states.
_Avoid_: Generic graph, N-state DSL

**Structural Exchange Edge**:
An unordered state pair present in a kinetic model's declared topology. It owns
one public `KEX_ij` and the derived directional pair `Kij`/`Kji`; setting its
`KEX_ij` to zero switches the edge off dynamically without making it absent.
_Avoid_: Enabled edge, nonzero edge

**Population Simplex**:
The closed scientific domain in which every state population is nonnegative and
the complete population vector sums to one. Generic N-state models expose the
non-A populations and derive `PA` as the validated complement.
_Avoid_: Population normalization, softmax coordinates

**Feasible Coordinates**:
Private solver coordinates compiled from model-owned scientific domains. They
keep public parameter names and continuation values unchanged while ensuring a
local optimizer proposes only representable domain states.
_Avoid_: Public reparameterization, clipping

**Evidence**:
The validated source observations or samples of a statistical analysis together
with their scientific state.
_Avoid_: Summary, Analysis Result

**Summary**:
Scientific conclusions and statistics derived from Evidence under an explicit
interpretation policy.
_Avoid_: Evidence, report

**Analysis Result**:
The complete authoritative outcome of one specific statistical analysis,
including its Evidence, Summary availability, and scientific completeness.
Resampling and MCMC have distinct Analysis Results.
_Avoid_: Evidence, Summary, output artifacts

**Deterministic Uncertainty**:
The complete interpreted uncertainty associated with one accepted deterministic
fit, including its Evidence and predetermined availability, completeness, and
reportability.
_Avoid_: Product Uncertainty, Uncertainty Result, publication artifacts

**Known ChemEx Failure**:
A termination that ChemEx deliberately classifies with trusted information
intended for concise user presentation.
_Avoid_: Any exception raised while ChemEx is running

**Unexpected Internal Failure**:
An unclassified failure for which ChemEx cannot responsibly provide a known
user or scientific explanation.
_Avoid_: Known ChemEx Failure, validation failure

**User Interruption**:
A user-requested stop whose scientifically valid committed state or Evidence may
be finalized before termination.
_Avoid_: Failed analysis, cancellation
