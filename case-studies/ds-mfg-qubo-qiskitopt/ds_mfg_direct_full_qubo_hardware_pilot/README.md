# Direct Full-QUBO IBM Hardware Pilot

This directory caches the first real IBM Runtime run of the direct 36-variable
DS-MFG full-QUBO p=2 QAOA circuit. This run was collected with the original
hand-built direct-QAOA hardware circuit path before the hardware handoff script
was switched to `QiskitOpt.QAOA.fixed_parameter_circuit`. Retain it as
descriptive legacy hardware evidence, not as a validated comparison to the
QiskitOpt high-read noiseless reference. A comparable hardware run should be
resubmitted with the current script.

The run used the optimized high-read parameter artifact:

```text
ds_mfg_direct_full_qubo_qaoa_highread/direct_full_qubo_qaoa_optimized_parameters.json
```

Run summary:

- backend: `ibm_fez`,
- submitted jobs: 1,
- shots: 4096,
- transpile seed: 92001,
- QAOA depth: p=2,
- top-50 repaired-flow hits: 0,
- top-10 repaired-flow hits: 0,
- global-optimum repaired-flow hits: 0,
- best repaired objective: `91.7295`,
- best repaired flow: `1010110110000000000`.

The cached files are:

- `job_manifest.json`: sanitized run manifest, submitted job id, and QAOA
  parameter metadata.
- `backend_metadata.json`: resolved backend metadata for `ibm_fez`.
- `raw_counts.csv`: raw full-bitstring counts.
- `scored_counts.csv`: exact auxiliary-repair scoring of each observed
  full-bitstring.
- `summary.csv`: aggregate hit counts and hit-rate uncertainty columns.

These artifacts are descriptive legacy hardware evidence only. They are not
evidence of quantum speedup or IBM hardware-performance superiority. The run
found no top-50 repaired-flow hit in 4096 shots, whereas the corresponding
noiseless high-read Aer sample found 4 repaired global-optimum hits in 32768
reads.
