# Direct Full-QUBO Audit

This directory caches a direct audit of the original exported 36-variable
DS-MFG QUBO from `Fw_ DS mfg case qubo information.zip`.

- `direct_full_qubo_summary.csv` rescored cached direct full-QUBO QAOA/VQE
  distributions against exact auxiliary repair, the top-50/top-10/global
  repaired-flow table, encoded Gurobi-pool consistency, and Gurobi-pool
  feasibility.
- `direct_full_qubo_resource_summary.csv` records the direct full-QUBO p=2 QAOA
  circuit resource summary and gated FakeFez transpile metadata.
- `direct_full_qubo_resource_metadata.json` records the same resource audit in
  structured form.

Regenerate the cache from the case-study directory with:

```bash
julia --project=. scripts/run_direct_full_qubo_audit.jl
```

Regenerate with FakeFez transpilation metadata:

```bash
DSMFG_DIRECT_FULL_QUBO_TRANSPILE=true julia --project=. scripts/run_direct_full_qubo_audit.jl
```

No direct 36-qubit noisy samples or IBM hardware jobs are cached here.
