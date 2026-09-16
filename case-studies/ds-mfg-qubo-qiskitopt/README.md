# DS-MFG QUBO Case Study Bundle

This bundle contains a reproducible DS-MFG QUBO case study using local classical
emulations of QAOA and VQE through `QiskitOpt.jl` and Qiskit Aer, plus a small
IBM hardware pilot for the fixed-parameter QAOA handoff. It includes the
executed notebook, cached distributions, generated figures, hardware pilot
artifacts, the manuscript and supplement, and expanded findings.

The case study is framed as a practical quantum-optimization audit. Gurobi
provides the original 19-flow-variable ground truth, local Aer emulation
provides the main QAOA/VQE evidence, and IBM hardware execution is reported as
a descriptive feasibility pilot. The bundle does not claim quantum speedup or
IBM hardware-performance superiority.

Start with:

- [main_standalone.pdf](main_standalone.pdf): the submitted paper, including
  the reviewer robustness and cost revision. The conference submission has no
  supplement; [supplement_standalone.pdf](supplement_standalone.pdf) is
  repository documentation that the paper cites as the shared repository
  (repair component table, environment, and validation tables).
- `README_STUDENT.md`: detailed guide to the notebook, environment, cached
  results, and rerun commands.
- `ds_mfg_qubo_qiskitopt.ipynb`: executed notebook with the full workflow.
- `MANUSCRIPT_FINDINGS.tex`: expanded LaTeX findings with introduction, method,
  results, discussion, conclusion, experiment provenance, result tables, and
  the IBM hardware pilot.
- `MANUSCRIPT_FINDINGS.md`: compact narrative summary.

## Submitted Manuscript

The FOCAPO-CPC 2027 submission is titled *Gate-Based Quantum Optimization for
Discrete Process Design: A Repair-Aware Benchmark*. Its sources and supporting
files were imported from `FOCAPO_CPC_2027_YPark (1).zip` (SHA-256:
`245e18e2d24fa9e3bf6dd7e51995f92158be46fd9faf070b550643cb949adc97`).

| Submission archive | Repository file |
| --- | --- |
| `main.tex` | [main_standalone.tex](main_standalone.tex) |
| `supplement.tex` | [supplement_standalone.tex](supplement_standalone.tex) |
| All other entries | Same relative paths in this directory |

Commit `93c4f97e13936c2434c84afe5ebf8b8a4794d722` preserves all imported files
exactly, including `bibliography.bib`,
`psepress.sty`, `preprint.cls`, the figures, bundled fonts, `PSEkeywords.txt`,
and `.vscode/settings.json`. The ZIP itself is a local export excluded by the
repository's existing ignore rules. The tracked PDFs are local rebuilds; the
submission ZIP contains no PDFs. The current manuscripts extend that submission
with the [reviewer robustness and cost audit](ds_mfg_revision_audit/README.md),
within one additional main-manuscript page (six pages to seven).
`MANUSCRIPT_FINDINGS.tex` and
`MANUSCRIPT_FINDINGS.md` retain the expanded analysis and experiment provenance
that support the submitted paper.

To rebuild both PDFs, install LuaLaTeX, Biber, and latexmk, then run from this
directory:

```bash
latexmk -lualatex -interaction=nonstopmode -halt-on-error main_standalone.tex
latexmk -lualatex -interaction=nonstopmode -halt-on-error supplement_standalone.tex
```

The style loads its fonts from `fonts/`. Intermediate LaTeX files are ignored
by Git.

## Modeling and Cached Results

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
  original 512-read direct full-QUBO QAOA row found 5 top-50 repaired-flow
  reads, 1 top-10 read, and no global reads; a higher-read optimized p=2 direct
  QAOA run in `ds_mfg_direct_full_qubo_qaoa_highread/` found 4 global-optimum
  repaired-flow hits in 32768 final reads. The best direct full-QUBO VQE seed
  found 2 top-50 reads, 1 top-10 read, and no global reads in 8192 final
  samples.
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
- Direct full-QUBO hardware pilot: the optimized direct 36-variable full-QUBO
  p=2 QAOA circuit was submitted to `ibm_fez` for 1 job of 4096 shots. It found
  0 top-50, 0 top-10, and 0 global repaired-flow hits; its best repaired
  objective was `91.7295`. This was collected with the original hand-built
  hardware-circuit path and is retained as legacy descriptive evidence; rerun
  with the current `QiskitOpt.QAOA.fixed_parameter_circuit` handoff before using
  it as a validated comparison. Cached artifacts are in
  `ds_mfg_direct_full_qubo_hardware_pilot/`.
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

When present, the hardware pilot defaults to the optimized high-read artifact
`ds_mfg_direct_full_qubo_qaoa_highread/direct_full_qubo_qaoa_optimized_parameters.json`.
The older audit artifact
`ds_mfg_direct_full_qubo_audit/direct_full_qubo_qaoa_hardware_parameters.json`
remains available as the warm-start parameter record. The hardware pilot
dry-runs by default. Real IBM submission requires configured Runtime credentials
and the explicit `DSMFG_RUN_DIRECT_FULL_QUBO_HARDWARE=true` gate.

The cached resource row records 36 logical qubits, 202 logical `rzz` gates, a
FakeFez-transpiled depth of 1344 with 1157 `cz` gates, and no direct noisy
sample cache. The direct full-QUBO hardware pilot is a real hardware sample,
not a noisy simulator cache. A dense noiseless statevector would require `2^36` amplitudes
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
