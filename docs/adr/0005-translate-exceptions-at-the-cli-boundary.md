---
status: accepted
---

# Translate exceptions at the CLI boundary

ChemEx will keep internal Python interfaces exception-transparent while one
lightweight terminal boundary owns command-line failure translation.  The
installed `chemex` command and `python -m chemex` share that boundary.  It
covers application imports reached from the boundary as well as failures raised
during command execution.

`SystemExit` remains native so argparse and intentional command exits retain
their established status and stream behavior.  A user interruption exits 130.
Other ordinary failures exit non-zero.  Terminal diagnostics are concise,
literal, written to stderr, and contain neither traceback nor chained exception
details.  Rendering failure falls back to plain stderr without replacing the
original exit status.

Workflow layers may catch interruption temporarily to freeze scientifically
valid partial state, serialize already committed values or qualified Evidence,
write truthful incomplete status, and perform required cleanup.  Once
interruption is recognized, finalization does not start further scientific
calculation, sampling, replicates, Method Steps, or optional plotting.  The
interruption remains distinguishable from ordinary failure and propagates to
the terminal boundary after finalization.

## Considered options

- Catching failures in every command or scientific layer was rejected because
  it duplicates presentation policy and risks converting failures into success.
- Immediately re-raising every raw `KeyboardInterrupt` was rejected because it
  can discard qualified Evidence and committed results before truthful
  publication.
- A new global exception hierarchy was rejected because existing workflow
  terminals already distinguish interrupted operations from ordinary failure.
- Leaving application imports outside the boundary was rejected because common
  dependency and startup failures would still expose tracebacks.

## Consequences

Direct calls to `chemex.chemex.main()` and lower-level APIs continue to raise
their original exceptions.  Workflow-specific incomplete exceptions may carry
their existing normalized `terminal` value through orchestration.  Completed
deterministic fits, qualified uncertainty Evidence, complete MCMC transitions,
and completed resampling replicates survive a later interruption without being
promoted to complete scientific conclusions.

The boundary can only translate failures after Python has imported the ChemEx
package and reached it.  Interpreter startup failures, an unavailable package,
bootstrap syntax errors, native process crashes, and forced process termination
remain outside this guarantee.
