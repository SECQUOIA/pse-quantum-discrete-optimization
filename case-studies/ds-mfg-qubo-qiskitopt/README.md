# DS-MFG QUBO Case Study Bundle

This bundle contains a reproducible DS-MFG QUBO case study using local classical
emulations of QAOA and VQE through `QiskitOpt.jl` and Qiskit Aer, plus a small
IBM hardware pilot for the fixed-parameter QAOA handoff. It includes the
executed notebook, cached distributions, generated figures, hardware pilot
artifacts, and a manuscript-ready findings draft.

The case study is framed as a practical quantum-optimization audit. Gurobi
provides the original 19-flow-variable ground truth, local Aer emulation
provides the main QAOA/VQE evidence, and IBM hardware execution is reported as
a descriptive feasibility pilot. The bundle does not claim quantum speedup or
IBM hardware-performance superiority.

Start with:

- `README_STUDENT.md`: detailed guide to the notebook, environment, cached
  results, and rerun commands.
- `ds_mfg_qubo_qiskitopt.ipynb`: executed notebook with the full workflow.
- `MANUSCRIPT_FINDINGS.tex`: LaTeX manuscript body with introduction, method,
  results, discussion, conclusion, experiment provenance, result tables, and
  the IBM hardware pilot.
- `MANUSCRIPT_FINDINGS.md`: compact narrative summary.

The main modeling point is that the QUBO has 36 binary variables, but Gurobi's
ground truth is defined over the original 19 flow variables. The remaining 17
bits are auxiliary/slack variables. A sampled QUBO state can therefore project
to a good Gurobi flow while still carrying an inconsistent auxiliary assignment.
The notebook scores samples using exact auxiliary repair, which fixes the flow
bits and chooses the minimum-energy auxiliary completion.

Best cached results:

- Gurobi global optimum: objective `11.7095`, flow
  `1001110100111100011`.
- QAOA: top-10-targeted transferred p=5 JuliQAOA angles on the reduced 19-flow
  surrogate sampled the global optimum 1007 times in 262144 local Aer reads,
  improving on the earlier 664-hit energy-targeted endpoint.
- VQE: the final reduced-surrogate follow-up sampled the global optimum 28
  times across three 524288-read local Aer runs; seed 74018 produced 20 of
  those reads.
- Direct original-QUBO audit: cached direct 36-variable QUBO QAOA/VQE
  distributions are consolidated in `ds_mfg_direct_full_qubo_audit/`. The best
  direct full-QUBO QAOA row found 5 top-50 repaired-flow reads, 1 top-10 read,
  and no global reads in 512 final samples; the best direct full-QUBO VQE seed
  found 2 top-50 reads, 1 top-10 read, and no global reads in 8192 final
  samples. The same audit persists the direct p=2 QAOA parameter artifact used
  by the dry-run-safe full-QUBO hardware handoff script.
- Classical baselines: uniform random repaired-flow sampling found no global
  optimum hits at 262144 or 524288 samples, while the 262144-evaluation
  hill-climb restart baseline found 19 global hits. Cached outputs are in
  `ds_mfg_classical_baselines/`.
- Hit-rate uncertainty: the QAOA, reduced-surrogate VQE, and classical-baseline
  summary CSVs include top-50, top-10, and global hit rates with 95% Wilson
  intervals, feasible Gurobi-pool hit rates, and empirical 99%
  time-to-solution columns. A consolidated report is in
  `ds_mfg_hit_rate_reports/time_to_solution_report.csv`.
- IBM hardware pilot: the fixed top-10-targeted p=5 QAOA circuit was submitted
  to `ibm_fez` for 9 jobs of 4096 shots each. Across 36864 hardware reads, the
  run found 6 top-50 repaired-flow hits, 1 top-10 hit, 0 global-optimum hits,
  and a best repaired objective of `11.8105` at rank 2. Cached artifacts are in
  `ds_mfg_ibm_qaoa_pilot_fez_4096x3/`.
- Simulator-to-hardware comparison: cached tables in
  `ds_mfg_simulator_hardware_comparison/` compare the ideal Aer top-10-targeted
  p=5 transfer, the model-based FakeFez/Aer noisy simulation from
  `ds_mfg_fake_fez_qaoa_noisy_4096x3/summary.csv`, and the `ibm_fez` pilot.

The IBM hardware pilot is reported descriptively. It is not evidence of quantum
speedup or IBM hardware-performance superiority.

To refresh the direct original-QUBO audit cache without running noisy samples:

```bash
julia --project=. scripts/run_direct_full_qubo_audit.jl
```

To also record gated FakeFez transpilation metadata for the direct full-QUBO
p=2 circuit:

```bash
DSMFG_DIRECT_FULL_QUBO_TRANSPILE=true julia --project=. scripts/run_direct_full_qubo_audit.jl
```

To build the direct full-QUBO p=2 QAOA hardware handoff artifacts without
submitting jobs:

```bash
QISKIT_IBM_BACKEND=ibm_fez julia --project=. scripts/run_direct_full_qubo_hardware_pilot.jl
```

The hardware pilot consumes
`ds_mfg_direct_full_qubo_audit/direct_full_qubo_qaoa_hardware_parameters.json`
and dry-runs by default. Real IBM submission requires configured Runtime
credentials and the explicit `DSMFG_RUN_DIRECT_FULL_QUBO_HARDWARE=true` gate.

The cached resource row records 36 logical qubits, 202 logical `rzz` gates, a
FakeFez-transpiled depth of 1344 with 1157 `cz` gates, and no direct noisy
sample cache. A dense noiseless statevector would require `2^36` amplitudes
(`1099511627776` complex128 bytes), while exact dense noisy density-matrix
simulation would require `2^72` complex entries. Aer automatic or MPS-style
noisy simulation may avoid dense density-matrix storage for some circuits, but
the cost is entanglement-dependent and was not run here at a useful shot budget.

To rerun the model-based noisy simulation as an overnight job, run:

```bash
DSMFG_RUN_NOISY_SIMULATION=true julia --project=. scripts/run_noisy_qaoa_fake_backend.jl
julia --project=. scripts/update_simulator_hardware_comparison.jl
```

The default noisy run uses `qiskit_ibm_runtime.fake_provider.FakeFez`,
`qiskit_aer.AerSimulator.from_backend`, 4096 shots, three repeats, and
transpile seeds `92001:92003`, matching the hardware pilot's nominal shot
budget. Treat its output as model-based simulation, not as a calibrated
prediction of the hardware run.
