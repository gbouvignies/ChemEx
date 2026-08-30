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
