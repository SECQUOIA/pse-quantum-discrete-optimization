# Direct Full-QUBO Audit

This directory caches a direct audit of the original exported 36-variable
DS-MFG QUBO from `Fw_ DS mfg case qubo information.zip`.

- `direct_full_qubo_summary.csv` rescored cached direct full-QUBO QAOA/VQE
  distributions against exact auxiliary repair, the top-50/top-10/global
  repaired-flow table, encoded Gurobi-pool consistency, and Gurobi-pool
  feasibility.
- `direct_full_qubo_resource_summary.csv` records the direct full-QUBO p=2 QAOA
  circuit resource summary and gated FakeFez transpile metadata.
- `direct_full_qubo_qaoa_hardware_parameters.json` persists the direct p=2 QAOA
  depth, beta-then-gamma angle vector, source distribution, and hardware
  handoff policy used by `scripts/run_direct_full_qubo_hardware_pilot.jl`.
- `direct_full_qubo_resource_metadata.json` records the same resource audit in
  structured form. The resource cache omits wall-clock generation timestamps and
  records transpile runtime as `not_recorded_for_reproducible_artifacts` so
  rerunning the gated transpile can regenerate byte-stable artifacts when
  Qiskit output is unchanged.

The noisy-simulation feasibility note is about exact dense noise simulation,
not only the noiseless `2^36` statevector. A dense 36-qubit density matrix has
`2^72` complex entries. Aer automatic or MPS-style noisy simulation may avoid
that storage for some circuits, but its cost is entanglement-dependent; this
cache does not establish a useful-shot direct noisy run as practical.

Regenerate the cache from the case-study directory with:

```bash
julia --project=. scripts/run_direct_full_qubo_audit.jl
```

Regenerate with FakeFez transpilation metadata:

```bash
DSMFG_DIRECT_FULL_QUBO_TRANSPILE=true julia --project=. scripts/run_direct_full_qubo_audit.jl
```

Build the direct full-QUBO hardware handoff artifacts without submitting jobs
with:

```bash
QISKIT_IBM_BACKEND=ibm_fez julia --project=. scripts/run_direct_full_qubo_hardware_pilot.jl
```

No direct 36-qubit noisy samples or IBM hardware jobs are cached here. The
hardware script dry-runs unless `DSMFG_RUN_DIRECT_FULL_QUBO_HARDWARE=true` is
set with IBM Runtime credentials configured outside the repository.
