# Install On Another Device

This repository is organized as reproducible case-study artifacts. The current
installable project is the DS-MFG QUBO case study in
`case-studies/ds-mfg-qubo-qiskitopt`.

## Prerequisites

- Git.
- Julia 1.10 or newer. Julia 1.10 LTS is recommended because the checked-in
  manifest was generated with Julia 1.10.11.
- Internet access for the first install so Julia packages, binary artifacts,
  and the managed Python/Qiskit environment can be downloaded.
- Optional: a `jupyter` command on `PATH` if you want to execute the full
  notebook with `scripts/nbconvert_ds_mfg.sh`.

You do not need to install Python, Conda, Qiskit, Gurobi, or IBM credentials
manually for the smoke test. The Python/Qiskit stack is described by
`CondaPkg.toml` and created locally by `PythonCall.jl`/`CondaPkg.jl` on first
use. The tracked DS-MFG results are local classical emulations, not IBM
hardware jobs.

## Fresh Clone Check

Clone the repository and run the case-study bootstrap script:

```bash
git clone https://github.com/SECQUOIA/pse-quantum-discrete-optimization.git
cd pse-quantum-discrete-optimization
case-studies/ds-mfg-qubo-qiskitopt/scripts/bootstrap_ds_mfg.sh
```

If Julia is not named `julia` on your `PATH`, set `JULIA_BIN`:

```bash
JULIA_BIN=/path/to/julia case-studies/ds-mfg-qubo-qiskitopt/scripts/bootstrap_ds_mfg.sh
```

The script instantiates the recorded Julia environment, precompiles it, creates
the managed Python/Qiskit environment on first package import, and runs
`scripts/smoke_test.jl`. A successful setup ends with:

```text
DS-MFG smoke test passed.
```

The smoke test checks that:

- the Julia packages import at the expected versions,
- Qiskit, Qiskit Aer, and Qiskit Optimization import through PythonCall,
- required notebooks, input archives, figures, and cached result directories
  are present,
- cached result summaries contain the expected DS-MFG optimum and sampler hit
  counts.

## Manual Equivalent

From the case-study folder, the bootstrap script is equivalent to:

```bash
cd case-studies/ds-mfg-qubo-qiskitopt
julia --project=. -e 'import Pkg; Pkg.instantiate(); Pkg.precompile()'
julia --project=. scripts/smoke_test.jl
```

## Running The Notebook

The executed notebook is already tracked at:

```text
case-studies/ds-mfg-qubo-qiskitopt/ds_mfg_qubo_qiskitopt.ipynb
```

To rerun it non-interactively after the smoke test, install or expose a
`jupyter` command and run:

```bash
cd case-studies/ds-mfg-qubo-qiskitopt
scripts/nbconvert_ds_mfg.sh
```

You can override tool locations with:

```bash
JULIA_BIN=/path/to/julia JUPYTER_BIN=/path/to/jupyter scripts/nbconvert_ds_mfg.sh
```

The full notebook uses local Qiskit Aer simulation and can take much longer
than the smoke test. Several expensive experiment cells use cached CSVs when
the tracked result files are present.

The baseline QAOA/VQE notebook solve uses fixed local simulator seeds, and the
tracked metadata avoids volatile timestamps, absolute machine-local paths, and
solve timings so rerunning the notebook does not dirty cached artifacts for
environment-only reasons.

## Local Generated Files

First runs may create local directories such as `.CondaPkg/`, `.jupyter/`, and
`.julia-depot-nbconvert/` inside the case-study folder. These are machine-local
environment caches and are intentionally ignored by Git.
