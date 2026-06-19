# DS-MFG QUBO QiskitOpt.jl Case Study

This folder contains a worked notebook for the DS-MFG discrete IP QUBO case. It
builds the QUBO from exported CSV files, runs local Aer emulations of QAOA and
VQE through `QiskitOpt.jl`, and compares the sampled distributions to the
Gurobi solution pool.

The notebook is a practical audit of the quantum-optimization workflow rather
than a speedup study. Gurobi supplies the original 19-flow-variable ground
truth, local Aer emulation supplies the main QAOA/VQE evidence, and IBM
hardware execution is represented by a small fixed-parameter QAOA feasibility
pilot. The current artifacts do not support a quantum speedup or IBM
hardware-performance superiority claim.

## Start Here

Read these files in order:

1. `MANUSCRIPT_FINDINGS.tex` for the manuscript-ready introduction, method,
   results, discussion, conclusion, formulation, reduction, and result
   narrative.
2. `ds_mfg_qubo_qiskitopt.ipynb` for the executable analysis record.
3. `MANUSCRIPT_FINDINGS.md` for a compact summary.

Open the notebook with:

```text
ds_mfg_qubo_qiskitopt.ipynb
```

The notebook has executed outputs and explanatory text. Read the first markdown
sections before running cells; they explain the file layout, the flow/auxiliary
bit convention, and the cached result folders.

## Fresh Clone Check

From this folder, instantiate the recorded Julia environment and run the smoke
test before attempting the full notebook:

```bash
scripts/bootstrap_ds_mfg.sh
```

This is equivalent to:

```bash
julia --project=. -e 'import Pkg; Pkg.instantiate(); Pkg.precompile()'
julia --project=. scripts/smoke_test.jl
```

The smoke test imports the Julia and Python packages used by the notebook,
checks that required artifacts exist, and validates key cached result summaries.
It does not rerun QAOA or VQE, so it should complete much faster than executing
the notebook.

## Main Takeaway

The main takeaway is that meaningful QAOA/VQE samples can be obtained for this
DS-MFG instance after problem-specific reduction, auxiliary repair,
simulator-configuration, and postprocessing work. The engineering burden is part
of the result, not just notebook setup.

The original Gurobi model is the 19-variable DS-MFG flow problem. The exported
QUBO has 36 binary variables: the first 19 bits are the original flow variables,
while the remaining 17 are auxiliary/slack bits introduced by the reformulation.
Gurobi reports solutions in the original 19 flow variables. A quantum sample can
therefore project to a Gurobi-pool flow but still have the wrong auxiliary bits
and a poor raw QUBO energy.

For this reason, the notebook reports:

- projected match: whether the 19 flow bits appear in the Gurobi pool,
- encoded match: whether the full sampled 36-bit QUBO state is also encoded
  correctly,
- repaired objective: the QUBO objective after fixing the sampled flow bits and
  choosing the exact best auxiliary completion.

## Best Cached Results

Gurobi global optimum:

- objective: `11.7095`
- flow bits: `1001110100111100011`

Best sampled QAOA results in this analysis bundle:

- full QAOA with larger final sampling: repaired pool objective `14.5505`,
- reduced-surrogate QAOA with top-10-targeted JuliQAOA-learned p=5 angles:
  1007 sampled global-optimum reads out of 262144. This improves on the earlier
  664-hit energy-targeted endpoint and is the strongest cached quantum-emulation
  result in the bundle.

Best sampled VQE results in this analysis bundle:

- full VQE with 8192 final reads over 5 seeds: repaired pool objective `14.6515`,
- reduced-surrogate VQE with 8192 final reads over 5 seeds: repaired pool objective `11.8105`,
- reduced-surrogate VQE with 32768 final reads over 20 seeds: 3 sampled
  global-optimum reads,
- reduced-surrogate VQE high-read follow-up with 262144 final reads over 3
  selected seeds: 19 sampled global-optimum reads,
- final reduced-surrogate VQE follow-up with 524288 final reads over the same 3
  selected seeds: 28 sampled global-optimum reads, with seed 74018 producing 20.

Direct original-QUBO audit in this analysis bundle:

- `scripts/run_direct_full_qubo_audit.jl` reloads the original archived
  36-variable QUBO and scores cached direct full-QUBO QAOA/VQE distributions
  with exact auxiliary repair,
- the best direct full-QUBO QAOA final-sampling row has 5 top-50 repaired-flow
  hits, 1 top-10 hit, 0 global hits, and best repaired objective `14.5505` in
  512 reads,
- the best direct full-QUBO VQE seed has 2 top-50 repaired-flow hits, 1 top-10
  hit, 0 global hits, and best repaired objective `14.6515` in 8192 reads,
- the direct full-QUBO p=2 resource audit records 36 qubits, 202 logical `rzz`
  gates, FakeFez-transpiled depth 1344, 1157 `cz` gates, and no cached direct
  noisy samples.

Classical sampling baselines in this analysis bundle:

- uniform random repaired-flow sampling with 262144 samples: 24 top-50 hits, 5
  top-10 hits, 0 global-optimum hits,
- uniform random repaired-flow sampling with 524288 samples: 49 top-50 hits, 11
  top-10 hits, 0 global-optimum hits,
- steepest-descent hill climbing with random restarts and 262144
  repaired-objective evaluations: 883 top-50 hits, 177 top-10 hits, 19
  global-optimum hits.

IBM hardware pilot in this analysis bundle:

- fixed top-10-targeted p=5 QAOA circuit on `ibm_fez`: 9 submitted jobs,
  4096 shots per job, 36864 total reads, 6 top-50 hits, 1 top-10 hit, 0
  global-optimum hits, and best repaired objective `11.8105` at rank 2.
- simulator-to-hardware comparison tables in
  `ds_mfg_simulator_hardware_comparison/`: the current cache compares ideal
  Aer, model-based FakeFez/Aer noisy, and hardware hit rates using the cached
  noisy summary in `ds_mfg_fake_fez_qaoa_noisy_4096x3/summary.csv`.

The strongest cached sampler is reduced-surrogate QAOA with top-10-targeted
angles learned offline using the local JuliQAOA statevector simulator and
transferred back into `QiskitOpt.QAOA` for Aer sampling. The later VQE follow-up
also reaches the optimum, but at much lower frequency. The cached CSVs are
important for comparing distributions without rerunning the quantum emulation.

## What Was Run

Most quantum-algorithm results in this bundle are local classical emulations.
The current runnable project uses `QiskitOpt v0.6.0` with local Qiskit Aer. The
high-cost simulations use Aer's matrix-product-state backend with the
package-level Aer options exposed by `QiskitOpt.jl`; the historical cached
distributions remain tracked so the comparisons can be inspected without
rerunning expensive emulations. The included hardware pilot is a separate
fixed-parameter sampling check on IBM hardware.

The strongest QAOA run used a two-stage workflow: JuliQAOA searched for p=5
angles on the reduced 19-flow quadratic surrogate using a top-10 objective, then
those angles were transferred into `QiskitOpt.QAOA` for fixed-parameter Aer
sampling with `MaximumIterations() = 0`.

The VQE runs used `EfficientSU2`, multiple random starts, separated optimizer
reads from final sampling reads, and the same exact auxiliary-repair scoring.

The model-based noisy QAOA runner is:

```bash
DSMFG_RUN_NOISY_SIMULATION=true julia --project=. scripts/run_noisy_qaoa_fake_backend.jl
julia --project=. scripts/update_simulator_hardware_comparison.jl
```

By default, it uses `FakeFez`, 4096 shots, three repeats, transpile seeds
`92001:92003`, simulator seeds starting at `93001`, and writes
`ds_mfg_fake_fez_qaoa_noisy_4096x3/`. This is intended for an overnight run:
local probes of the depth-1557, 1141-CZ transpiled circuit took 17.5 seconds
for 64 shots and 116.6 seconds for 512 shots. Any output from this script is a
model-based noisy simulation, not a calibrated prediction of the IBM hardware
pilot.

The classical baselines are generated by:

```bash
julia --project=. scripts/run_classical_baselines.jl
```

By default, the script writes `ds_mfg_classical_baselines/`. The summary CSV
records the full sample/evaluation budgets, hit counts, 95% Wilson hit-rate
intervals, feasible Gurobi-pool hit counts, and empirical 99%
time-to-solution columns. The compressed
distribution CSV records every unique sampled flow for each baseline run. A
separate retained-flow CSV keeps the sampled top-50 flows and each run's best
sampled flow for quick inspection. A consolidated QAOA/VQE/classical report is
written to `ds_mfg_hit_rate_reports/time_to_solution_report.csv`.

## Re-running The Notebook

The simplest path is to use Jupyter and select the kernel named:

```text
QiskitOpt Julia
```

For non-interactive execution, run from this folder:

```bash
scripts/nbconvert_ds_mfg.sh
```

The helper creates a local IJulia kernel and executes the notebook with
`nbconvert`. It uses the `julia` and `jupyter` commands on your `PATH`. Julia
1.10 or newer is required by `Project.toml`. If Julia or Jupyter are not on
your `PATH`, set:

```bash
export JULIA_BIN=/path/to/julia
export JUPYTER_BIN=/path/to/jupyter
scripts/nbconvert_ds_mfg.sh
```

The helper uses `--compiled-modules=no` for the IJulia kernel because this
environment had a ZMQ precompile failure otherwise.

The baseline QAOA/VQE notebook solve uses fixed local simulator seeds. The
tracked metadata intentionally avoids volatile timestamps, absolute
machine-local paths, and solve timings so reproducibility checks do not change
cached artifacts for environment-only reasons.

The local kernel environment includes both `IJulia` and `PythonCall`, and this
packaged project also includes `IJulia` as a notebook dependency. This is
intentional: when the notebook activates this folder's project and loads
`PythonCall`, Julia may also load IJulia's `IJuliaPythonCallExt` package
extension from the active project manifest. Keeping those dependencies recorded
avoids the extension-loading error that says to run `Pkg.instantiate()`.

## Package Setup

This folder uses the registered `QiskitOpt v0.6.0` package from Julia's General
registry. The notebook project is a normal Julia environment, not a local copy
of the package source.

The final sampling control used in the notebook is
`QUBODrivers.FinalNumberOfReads()`. It applies to both QAOA and VQE and lets the
final sampler use more reads without increasing the optimizer/estimator reads at
every iteration.

The local QAOA/VQE revisit scripts use `QUBODrivers.RandomSeed()` as the
standard sampler seed. With `QiskitOpt v0.6.0`, QiskitOpt derives local Aer
simulator and transpiler seeds from that value when explicit QiskitOpt seed
attributes are not set.

## Important Files

- `Fw_ DS mfg case qubo information.zip`: original QUBO/Gurobi input archive.
- `Project.toml`, `Manifest.toml`, `CondaPkg.toml`: Julia and Python environment files for the notebook.
- `MANUSCRIPT_FINDINGS.tex`: manuscript-ready LaTeX body with introduction,
  method, results, discussion, conclusion, experiment details, cached artifact
  paths, result tables, and the hardware pilot.
- `MANUSCRIPT_FINDINGS.md`: compact companion narrative for quick reading.
- `ds_mfg_*` folders: cached result distributions and summaries used by the notebook.
- `ds_mfg_ibm_qaoa_pilot_fez_4096x3/`: sanitized IBM hardware pilot manifest,
  backend metadata, raw counts, scored counts, and summary.
- `ds_mfg_simulator_hardware_comparison/`: cached ideal Aer, model-based noisy,
  and hardware comparison tables.
- `ds_mfg_*.svg`: figures generated by the notebook.
- `scripts/nbconvert_ds_mfg.sh`: full notebook execution helper.
- `scripts/find_reduced_qaoa_angles_juliqaoa.jl`: optional angle-search helper that expects a local sibling clone of `JuliQAOA.jl`, or a path set by `JULIQAOA_PROJECT`.
- `scripts/revisit_reduced_qaoa_juliqaoa.jl`: transfers cached JuliQAOA angles into `QiskitOpt.QAOA` and scores the Aer samples.
- `scripts/run_ibm_qaoa_pilot.jl`: dry-run-safe IBM Runtime pilot for the
  reduced 19-qubit surrogate using the persisted top-10-targeted p=5 JuliQAOA
  angles; set `DSMFG_RUN_IBM_HARDWARE=true` to submit real hardware jobs.
- `scripts/run_noisy_qaoa_fake_backend.jl`: dry-run-safe overnight FakeFez/Aer
  model-based noisy simulation runner for the same fixed QAOA circuit.
- `scripts/update_simulator_hardware_comparison.jl`: refreshes the cached
  ideal/noisy/hardware comparison tables after a noisy run completes.
- `scripts/revisit_reduced_vqe.jl`: reduced-surrogate VQE rerun helper; set `DSMFG_REVISIT_OUTPUT_DIR`, `DSMFG_REVISIT_SEEDS`, and `DSMFG_REVISIT_FINAL_READS` to store new sweeps separately.
- `scripts/run_classical_baselines.jl`: fixed-seed uniform random and
  hill-climb restart baselines scored by exact auxiliary repair.
- `scripts/hit_rate_stats.jl`: shared Wilson interval and empirical
  time-to-solution helper.
- `scripts/update_cached_hit_rate_reports.jl`: updates cached summary CSVs with
  hit-rate columns and writes the consolidated time-to-solution report.

## Open Research Directions

The remaining technical questions are about model reduction and sampler bias, not notebook mechanics:

- whether the reduced exact flow objective can be exposed through a non-quadratic or higher-order interface instead of a least-squares quadratic surrogate,
- whether the surrogate fit should be weighted toward low-energy flow assignments rather than all `2^19` assignments uniformly,
- whether JuliQAOA-derived schedules can be made less instance-specific and reused as warm starts across related DS-MFG instances,
- whether VQE performance is limited primarily by the EfficientSU2 ansatz, optimizer budget, or final sampling variance,
- whether exact auxiliary repair should become a standard postprocessing hook in the QUBO/QUBODrivers workflow for reformulations with separable auxiliary structure.
