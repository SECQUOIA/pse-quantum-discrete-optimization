# DS-MFG QUBO QiskitOpt.jl Case Study

This folder contains a worked notebook for the DS-MFG discrete IP QUBO case. It
builds the QUBO from exported CSV files, runs local Aer emulations of QAOA and
VQE through `QiskitOpt.jl`, and compares the sampled distributions to the
Gurobi solution pool.

## Start Here

Read these files in order:

1. `MANUSCRIPT_FINDINGS.tex` for the manuscript-ready formulation, reduction,
   and result narrative.
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
julia --project=. -e 'import Pkg; Pkg.instantiate()'
julia --project=. scripts/smoke_test.jl
```

The smoke test imports the Julia and Python packages used by the notebook,
checks that required artifacts exist, and validates key cached result summaries.
It does not rerun QAOA or VQE, so it should complete much faster than executing
the notebook.

## Main Takeaway

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
- reduced-surrogate QAOA with JuliQAOA-learned p=5 angles: 664 sampled
  global-optimum reads out of 262144. This is the strongest cached
  quantum-emulation result in the bundle and remains the accepted QAOA endpoint.

Best sampled VQE results in this analysis bundle:

- full VQE with 8192 final reads over 5 seeds: repaired pool objective `14.6515`,
- reduced-surrogate VQE with 8192 final reads over 5 seeds: repaired pool objective `11.8105`,
- reduced-surrogate VQE with 32768 final reads over 20 seeds: 3 sampled
  global-optimum reads,
- reduced-surrogate VQE high-read follow-up with 262144 final reads over 3
  selected seeds: 19 sampled global-optimum reads,
- final reduced-surrogate VQE follow-up with 524288 final reads over the same 3
  selected seeds: 28 sampled global-optimum reads, with seed 74018 producing 20.

The strongest cached sampler is reduced-surrogate QAOA with angles learned
offline using the local JuliQAOA statevector simulator and transferred back into
`QiskitOpt.QAOA` for Aer sampling. The later VQE follow-up also reaches the
optimum, but at much lower frequency. The cached CSVs are important for
comparing distributions without rerunning the quantum emulation.

## What Was Run

All quantum-algorithm results in this bundle are local classical emulations, not
hardware jobs. QAOA and VQE were run through `QiskitOpt v0.4.2` with local
Qiskit Aer. The high-cost simulations used Aer's matrix-product-state backend
with the package-level Aer options exposed by `QiskitOpt.jl`.

The strongest QAOA run used a two-stage workflow: JuliQAOA searched for p=5
angles on the reduced 19-flow quadratic surrogate, then those angles were
transferred into `QiskitOpt.QAOA` for fixed-parameter Aer sampling with
`MaximumIterations() = 0`.

The VQE runs used `EfficientSU2`, multiple random starts, separated optimizer
reads from final sampling reads, and the same exact auxiliary-repair scoring.

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

The local kernel environment includes both `IJulia` and `PythonCall`, and this
packaged project also includes `IJulia` as a notebook dependency. This is
intentional: when the notebook activates this folder's project and loads
`PythonCall`, Julia may also load IJulia's `IJuliaPythonCallExt` package
extension from the active project manifest. Keeping those dependencies recorded
avoids the extension-loading error that says to run `Pkg.instantiate()`.

## Package Setup

This folder uses the registered `QiskitOpt v0.4.2` package from Julia's General
registry. The notebook project is a normal Julia environment, not a local copy
of the package source.

The final sampling control used in the notebook is
`QUBODrivers.FinalNumberOfReads()`. It applies to both QAOA and VQE and lets the
final sampler use more reads without increasing the optimizer/estimator reads at
every iteration.

## Important Files

- `Fw_ DS mfg case qubo information.zip`: original QUBO/Gurobi input archive.
- `Project.toml`, `Manifest.toml`, `CondaPkg.toml`: Julia and Python environment files for the notebook.
- `MANUSCRIPT_FINDINGS.tex`: manuscript-ready LaTeX section with experiment details, cached artifact paths, result tables, and hardware placeholders.
- `MANUSCRIPT_FINDINGS.md`: compact companion narrative for quick reading.
- `ds_mfg_*` folders: cached result distributions and summaries used by the notebook.
- `ds_mfg_*.svg`: figures generated by the notebook.
- `scripts/nbconvert_ds_mfg.sh`: full notebook execution helper.
- `scripts/find_reduced_qaoa_angles_juliqaoa.jl`: optional angle-search helper that expects a local sibling clone of `JuliQAOA.jl`, or a path set by `JULIQAOA_PROJECT`.
- `scripts/revisit_reduced_qaoa_juliqaoa.jl`: transfers cached JuliQAOA angles into `QiskitOpt.QAOA` and scores the Aer samples.
- `scripts/revisit_reduced_vqe.jl`: reduced-surrogate VQE rerun helper; set `DSMFG_REVISIT_OUTPUT_DIR`, `DSMFG_REVISIT_SEEDS`, and `DSMFG_REVISIT_FINAL_READS` to store new sweeps separately.

## Open Research Directions

The remaining technical questions are about model reduction and sampler bias, not notebook mechanics:

- whether the reduced exact flow objective can be exposed through a non-quadratic or higher-order interface instead of a least-squares quadratic surrogate,
- whether the surrogate fit should be weighted toward low-energy flow assignments rather than all `2^19` assignments uniformly,
- whether JuliQAOA-derived schedules can be made less instance-specific and reused as warm starts across related DS-MFG instances,
- whether VQE performance is limited primarily by the EfficientSU2 ansatz, optimizer budget, or final sampling variance,
- whether exact auxiliary repair should become a standard postprocessing hook in the QUBO/QUBODrivers workflow for reformulations with separable auxiliary structure.
