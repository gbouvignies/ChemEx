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
