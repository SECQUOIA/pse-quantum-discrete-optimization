# Direct Full-QUBO QAOA High-Read Run

This directory caches a higher-read direct run of the original 36-variable
DS-MFG QUBO using `QiskitOpt.QAOA` with Aer matrix-product-state sampling.

- `direct_full_qubo_qaoa_highread_summary.csv` records the single retained run:
  p=2, seed 96001, 32 optimizer reads, 25 optimizer iterations, and 32768 final
  reads.
- `qaoa_direct_p2_seed96001_reads32768_distribution.csv` stores the scored final
  sample distribution with exact auxiliary repair.
- `qaoa_direct_p2_seed96001_reads32768_optimized_parameters.json` stores the
  optimizer-refined beta-then-gamma vector for this run.
- `direct_full_qubo_qaoa_optimized_parameters.json` is the canonical hardware
  handoff copy used by `scripts/run_direct_full_qubo_hardware_pilot.jl` when
  present.

The retained validation sample has 136 top-50 repaired-flow hits, 27 top-10
hits, and 4 global-optimum hits in 32768 final reads. The global flow is
`1001110100111100011` with repaired objective `11.7095`. The sampled auxiliary
bits for these global-flow hits are not encoded Gurobi-pool completions; the
global hits are established after exact auxiliary repair.

The optimized p=2 QAOA parameters are stored in Qiskit's beta-then-gamma order:

```text
0.39325958665050614;0.27183862464250103;-0.0086475034462802948;0.12401283690343518
```

Regenerate this cache from the case-study directory with:

```bash
DSMFG_DIRECT_HIGHREAD_FINAL_READS=32768 \
DSMFG_DIRECT_HIGHREAD_SEEDS=96001 \
DSMFG_DIRECT_HIGHREAD_PARAMETER_SOURCES=direct_p2 \
DSMFG_DIRECT_HIGHREAD_MAX_ITER=25 \
julia --project=. scripts/run_direct_full_qubo_qaoa_highread.jl
```

For direct hardware handoff, use fixed-parameter sampling with the optimized
artifact rather than running a hardware optimizer loop:

```bash
QISKIT_IBM_BACKEND=ibm_fez \
DSMFG_DIRECT_QAOA_PARAMETER_PATH=ds_mfg_direct_full_qubo_qaoa_highread/direct_full_qubo_qaoa_optimized_parameters.json \
julia --project=. scripts/run_direct_full_qubo_hardware_pilot.jl
```

The hardware pilot dry-runs by default. Real IBM Runtime submission requires
configured credentials and `DSMFG_RUN_DIRECT_FULL_QUBO_HARDWARE=true`.
