# PSE Quantum Discrete Optimization

This repository collects process systems engineering case studies that audit
quantum and quantum-inspired algorithms for discrete optimization models.

The initial case study is a demand-side manufacturing (DS-MFG) QUBO workflow
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
`ds_mfg_simulator_hardware_comparison/`.

The manuscript narrative keeps that evidence boundary explicit and does not
treat the pilot as a speedup or hardware-superiority result. The DS-MFG case
study includes `scripts/run_ibm_qaoa_pilot.jl` as a dry-run-safe IBM Runtime
runner; see `INSTALL.md` for credential handling and the explicit
`DSMFG_RUN_IBM_HARDWARE=true` submission gate. It also includes
`scripts/run_noisy_qaoa_fake_backend.jl` for an overnight model-based
FakeFez/Aer simulation; results from that script must be labeled as
model-based simulation, not as calibrated hardware predictions.
