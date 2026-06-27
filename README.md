# PSE Quantum Discrete Optimization

This repository collects process systems engineering case studies that audit
quantum and quantum-inspired algorithms for discrete optimization models.

The initial case study is a drug-substance manufacturing (DS-MFG) QUBO workflow
using local Qiskit Aer emulations of QAOA and VQE through `QiskitOpt.jl`.
It is framed as a practical workflow audit: Gurobi provides the 19-flow-variable
ground truth, local Aer emulation provides the current QAOA/VQE evidence,
simple classical sampling baselines provide non-quantum context, and a small
IBM hardware pilot provides a feasibility check for the fixed-parameter QAOA
handoff. The current bundle does not claim quantum speedup or IBM
hardware-performance superiority.
The repository layout is intentionally case-study based so that related PSE
instances can be added without changing the existing DS-MFG artifact paths.

## Case Studies

- `case-studies/ds-mfg-qubo-qiskitopt`: executed notebook, cached
  distributions, generated figures, scripts, and manuscript-ready findings for
  the DS-MFG QUBO case.

## Start Here

Clone the repository and enter the current DS-MFG case study:

```bash
git clone https://github.com/SECQUOIA/pse-quantum-discrete-optimization.git
cd pse-quantum-discrete-optimization
cd case-studies/ds-mfg-qubo-qiskitopt
```

Set up the Julia project and run the lightweight reproducibility check:

```bash
scripts/bootstrap_ds_mfg.sh
```

Equivalently, run the setup steps manually:

```bash
julia --project=. -e 'import Pkg; Pkg.instantiate(); Pkg.precompile()'
julia --project=. scripts/smoke_test.jl
```

The smoke test verifies that the Julia and Python dependencies import, the
case-study artifacts are present, the cached result summaries contain the
expected DS-MFG optimum, sampler hit counts including the objective-targeted
QAOA follow-up, hit-rate uncertainty columns, and time-to-solution report, and
the classical baseline distribution cache is complete. It is the fastest way for
a new user to confirm that a fresh clone is usable before rerunning the full
notebook.
See `INSTALL.md` for fresh-device prerequisites, expected generated files, and
notebook execution details.

## Reproducibility Commands

From `case-studies/ds-mfg-qubo-qiskitopt`, the primary reproducibility commands
are:

```bash
scripts/bootstrap_ds_mfg.sh
julia --project=. scripts/run_classical_baselines.jl
QISKIT_IBM_BACKEND=ibm_fez julia --project=. scripts/run_ibm_qaoa_pilot.jl
QISKIT_IBM_BACKEND=ibm_fez julia --project=. scripts/run_direct_full_qubo_hardware_pilot.jl
```

The first command is the fresh-clone smoke test. The second regenerates the
fixed-seed uniform-random and hill-climb classical baselines. The third is an
IBM pilot dry run by default: it builds and scores the fixed-parameter QAOA
handoff artifacts without submitting hardware jobs. Real IBM hardware
submission requires credentials configured outside this repository and the
explicit `DSMFG_RUN_IBM_HARDWARE=true` gate. The fourth dry-runs the direct
36-variable full-QUBO QAOA hardware handoff. It defaults to the optimized p=2
artifact that produced 4 repaired global-optimum hits in 32768 local Aer reads
when that artifact is present; real submission is separately gated by
`DSMFG_RUN_DIRECT_FULL_QUBO_HARDWARE=true`.

Then read:

1. `README_STUDENT.md` for the detailed guide to the notebook, environment,
   cached results, and rerun commands.
2. `MANUSCRIPT_FINDINGS.tex` for the manuscript-ready formulation and results.
3. `ds_mfg_qubo_qiskitopt.ipynb` for the executed analysis record.

To execute the full notebook non-interactively from the case-study directory:

```bash
scripts/nbconvert_ds_mfg.sh
```

The full notebook run uses local Qiskit Aer emulation and can take much longer
than the smoke test.

Most quantum-algorithm evidence tracked here is local classical emulation. The
bundle also includes a small IBM hardware pilot for the reduced 19-qubit
top-10-targeted p=5 QAOA circuit on `ibm_fez`: 9 jobs, 4096 shots each, 36864
total reads, 6 top-50 repaired-flow hits, 1 top-10 hit, and no global-optimum
hit. Cached simulator-to-hardware comparison tables are stored under
`ds_mfg_simulator_hardware_comparison/`. A direct 36-qubit full-QUBO p=2
hardware pilot using the optimized high-read parameters is also cached as
legacy descriptive evidence: 1 `ibm_fez` job, 4096 shots, and no
top-50/top-10/global repaired-flow hits. Rerun it with the current
`QiskitOpt.QAOA.fixed_parameter_circuit` handoff before using it as a validated
comparison.

The manuscript narrative keeps that evidence boundary explicit and does not
treat the pilot as a speedup or hardware-superiority result. The DS-MFG case
study includes `scripts/run_ibm_qaoa_pilot.jl` as a dry-run-safe IBM Runtime
runner; see `INSTALL.md` for credential handling and the explicit
`DSMFG_RUN_IBM_HARDWARE=true` submission gate. It also includes
`scripts/run_noisy_qaoa_fake_backend.jl` for an overnight model-based
FakeFez/Aer simulation; results from that script must be labeled as
model-based simulation, not as calibrated hardware predictions. The direct
full-QUBO p=2 QAOA hardware handoff is available as
`scripts/run_direct_full_qubo_hardware_pilot.jl`; it prefers
`ds_mfg_direct_full_qubo_qaoa_highread/direct_full_qubo_qaoa_optimized_parameters.json`
when present and also dry-runs unless its direct hardware gate is enabled.
