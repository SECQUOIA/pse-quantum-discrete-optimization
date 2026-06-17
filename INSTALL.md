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
use. The Qiskit package versions are pinned for reproducible local emulation.
The tracked DS-MFG results are local classical emulations, not IBM hardware
jobs.

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

## IBM Hardware Pilot

The DS-MFG case study includes a dry-run-safe fixed-parameter IBM QAOA pilot
runner:

```bash
cd case-studies/ds-mfg-qubo-qiskitopt
QISKIT_IBM_BACKEND=ibm_brisbane julia --project=. scripts/run_ibm_qaoa_pilot.jl
```

Dry-run mode is the default. It validates the reduced 19-qubit surrogate, loads
the persisted top-10-targeted p=5 JuliQAOA angles, builds the fixed-parameter
QAOA circuit, and writes `job_manifest.json`, `backend_metadata.json`,
`raw_counts.csv`, `scored_counts.csv`, and `summary.csv` without contacting IBM
Runtime or submitting hardware jobs.

The runner is configured only through environment variables. The backend must
be one that the configured IBM Runtime instance can access; the tracked
hardware pilot used `ibm_fez`.

- `QISKIT_IBM_BACKEND`: required backend name, for example `ibm_fez` or another
  QPU visible to the configured instance.
- `QISKIT_IBM_CHANNEL`: optional Runtime channel, default
  `ibm_quantum_platform`.
- `QISKIT_IBM_INSTANCE`: IBM Runtime instance selector; optional only when
  Runtime can auto-resolve the account's instance.
- `DSMFG_HARDWARE_FINAL_READS`: shots per submitted sampler job, default `4096`.
- `DSMFG_HARDWARE_REPEATS`: repeat count, default `1`.
- `DSMFG_HARDWARE_TRANSPILE_SEEDS`: comma-separated seeds or ranges, default
  `92001`.
- `DSMFG_HARDWARE_OUTPUT_DIR`: output folder, default
  `ds_mfg_ibm_qaoa_pilot`.
- `DSMFG_RUN_IBM_HARDWARE`: must be `true` to submit jobs.

To run on real hardware, configure IBM credentials outside the repository,
either through IBM Runtime's normal account storage under your home directory
or through shell environment variables for the current session, for example
`QISKIT_IBM_TOKEN`. Do not place tokens, account JSON files, or backend secrets
anywhere in this repository.

After credentials are configured, the hardware gate must be explicit:

```bash
cd case-studies/ds-mfg-qubo-qiskitopt
export QISKIT_IBM_TOKEN='your-token'
QISKIT_IBM_BACKEND=ibm_fez \
QISKIT_IBM_CHANNEL=ibm_quantum_platform \
DSMFG_HARDWARE_FINAL_READS=4096 \
DSMFG_HARDWARE_REPEATS=1 \
DSMFG_HARDWARE_TRANSPILE_SEEDS=92001 \
DSMFG_RUN_IBM_HARDWARE=true \
julia --project=. scripts/run_ibm_qaoa_pilot.jl
```

The output files record backend name, job IDs, count data, repaired-flow
scoring, and summary hit rates. They intentionally do not write IBM tokens,
account-file paths, or instance values.
