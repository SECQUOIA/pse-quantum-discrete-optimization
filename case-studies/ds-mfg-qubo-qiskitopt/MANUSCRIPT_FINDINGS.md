# DS-MFG QUBO Case Study: Findings Narrative

For manuscript incorporation, use `MANUSCRIPT_FINDINGS.tex`. That file is a
LaTeX-ready manuscript body organized as Introduction, Method, Results,
Discussion, and Conclusion, with experiment provenance, cached artifact paths,
result tables, and the IBM hardware pilot. This Markdown file is kept as a
compact companion summary.

## Purpose

This case study evaluates local classical emulations of QAOA and VQE, accessed
through `QiskitOpt.jl`, on a DS-MFG discrete integer-programming QUBO instance.
The analysis compares sampled quantum-optimization distributions against a
Gurobi solution pool and uses the original IP variables to interpret the sampled
QUBO bitstrings.

The central claim is a practical audit claim. Current QAOA/VQE tooling can be
made to sample meaningful DS-MFG process-design candidates, but only after
substantial formulation, reduction, parameter-search, backend-configuration,
and postprocessing work. The evidence does not support a speedup claim over
Gurobi, and it does not include IBM hardware-performance superiority evidence.

The analysis keeps three evidence tiers separate:

- Gurobi ground truth: the original 19-flow-variable integer program and
  solution pool define the optimum and scoring reference. The retained Gurobi
  export lacks solver version/status/gap fields, so a deterministic enumeration
  of all `2^19` flow assignments confirms the 36 feasible flows and the
  `11.7095` optimum. A local Gurobi 13.0.2 rerun also returns `OPTIMAL`,
  `MIPGap=0.0`, 36 solutions, and exact agreement with the retained pool.
- Local Aer emulation: `QiskitOpt.jl`/Qiskit Aer QAOA and VQE runs test whether
  simulated sampling distributions place observable mass on useful repaired
  flow assignments.
- IBM hardware feasibility: a fixed-parameter QAOA pilot on `ibm_fez` tests the
  Runtime handoff and reports returned counts descriptively. It is not evidence
  of quantum speedup or hardware-performance superiority.

## Problem Structure

The original Gurobi model is named `capex`. It is a 19-variable binary flow
model with 20 constraints and a minimization objective. Gurobi returned a
36-solution pool over those 19 original flow variables.

The exported QUBO has 36 binary variables. The first 19 variables are the
original flow decisions from the integer program, and the remaining 17 variables
are auxiliary/slack bits introduced by the QUBO reformulation. The QUBO objective
is stored as `scale * (L'x + x'Qx + offset)` in `scalars.csv`, `L_vector.csv`,
and `Q_matrix.csv`. The Gurobi solution pool is indexed by the 19 flow variables,
not by the full 36-bit QUBO encoding.

This distinction is central. A sampled QUBO bitstring can project to a
Gurobi-pool flow while still assigning the auxiliary variables incorrectly. Such
a sample is meaningful at the IP-variable level, but its raw QUBO energy can be
poor. The notebook therefore reports three related comparisons:

- projected match: whether the first 19 sampled bits appear in the Gurobi pool,
- encoded match: whether the sampled auxiliary bits also match the
  minimum-energy auxiliary completion for that flow,
- repaired objective: the QUBO objective after fixing the sampled flow bits and
  choosing the exact best auxiliary completion.

The Gurobi reference optimum has repaired objective `11.7095` and flow bitstring
`1001110100111100011`.

## Auxiliary Repair And Reduced Objective

For a fixed flow assignment `y`, the notebook defines the repaired objective
`F(y) = min_z E_QUBO(y,z)`, where `z` contains the 17 auxiliary variables. For
this instance, once the 19 flow variables are fixed, the 17 auxiliary variables
decompose into 13 independent components. Each component contains either one or
two auxiliary bits. This makes exact auxiliary elimination cheap: each sampled
flow assignment can be repaired by solving small independent auxiliary
subproblems. This repair is the appropriate bridge between the QUBO formulation
and the original IP objective.

The notebook also constructs an exact reduced 19-flow objective by enumerating
all `2^19` flow assignments and applying exact auxiliary repair. The reduced
objective recovers the Gurobi reference optimum. Because the current QAOA/VQE
interfaces accept quadratic QUBOs, the exact reduced objective is approximated
by a least-squares quadratic surrogate over all flow assignments. The surrogate
has high global fit quality (`R^2 = 0.9990819`), but its deterministic minimum is
a near-optimal Gurobi-pool flow with repaired objective `14.5505`, not the true
reference optimum. This mismatch explains why surrogate-based quantum sampling must
still be evaluated against the exact repaired objective.

## QAOA Findings

The baseline full-QUBO QAOA and fixed-angle QAOA sweeps did not reliably produce
useful encoded solutions. Fixed-angle QAOA with p=2 produced one projected
Gurobi-pool read in 128 shots, but the auxiliary bits were not encoded
correctly. Exact auxiliary repair scored that projected flow at `18.0745`.

Separating optimizer reads from final sampling improved the full 36-variable
QAOA run. A p=2 run using 128 optimizer reads and 512 final reads produced three
projected pool reads and improved the best repaired pool objective to `14.5505`.
The direct original-QUBO audit cache rescored this same run against exact repair
and the top-flow table: it contains 5 top-50 repaired-flow reads, 1 top-10 read,
0 reference reads, and 0 encoded Gurobi-pool reads.

A higher-read direct full-QUBO p=2 run starting from the same parameter source
and allowing 25 optimizer iterations found the repaired reference optimum 4 times
in 32768 final reads. The optimized beta-then-gamma vector is cached in
`ds_mfg_direct_full_qubo_qaoa_highread/direct_full_qubo_qaoa_optimized_parameters.json`.
These are repaired-flow reference hits; the sampled auxiliary bits were not the
optimal auxiliary completions.

The highest-concentration QAOA result came from moving angle search off the Qiskit/Aer optimizer path. The workflow was:

1. eliminate auxiliary variables exactly to build the reduced 19-flow objective,
2. fit a quadratic surrogate over the exact reduced objective,
3. use `JuliQAOA.jl` as a local statevector angle finder on the reduced surrogate,
4. transfer the learned p=5 QAOA angles into `QiskitOpt.QAOA`,
5. set `MaximumIterations() = 0` and use Aer only for final sampling.

The original energy-targeted transferred-angle QAOA run sampled the Gurobi
reference optimum 664 times in 262144 final reads. Objective-targeted follow-up
searches then tested top-50, top-10, and direct-reference JuliQAOA objectives. The
top-50 and top-10 p=5 targets improved the high-read Aer result to 979 and 1007
reference reads, respectively, in 262144 final reads. The bounded direct-reference
p=3 check produced 45 reference reads. The top-10-targeted p=5 transfer is now the
highest-concentration cached quantum-emulation result in the notebook. Its best repaired
objective is `11.7095`, matching the Gurobi optimum, and the sampled optimal
flow is `1001110100111100011`.

## VQE Findings

The baseline VQE run was weak: it collapsed to a small set of sampled states and
did not provide a useful Gurobi-pool hit. The VQE experiments were then
strengthened in two ways:

- optimizer reads were separated from final sampling reads through `QUBODrivers.FinalNumberOfReads()`,
- several random EfficientSU2 initializations were tested.

On the full 36-variable QUBO, a 5-seed VQE sweep with 8192 final reads found a
projected Gurobi-pool flow with repaired objective `14.6515`. On the reduced
19-flow surrogate, a 5-seed sweep improved the best repaired objective to
`11.8105`, close to the reference optimum.
The direct original-QUBO audit cache shows that the best full-QUBO VQE seed had
2 top-50 repaired-flow reads, 1 top-10 read, 0 reference reads, and 0 encoded
Gurobi-pool reads.

After package-level Aer MPS backend attributes were available, the
reduced-surrogate VQE study was repeated with more seeds and higher final
sampling. A 20-seed sweep with 32768 final reads per seed sampled the Gurobi
reference optimum 3 times. A high-read follow-up on selected seeds with 262144
final reads sampled the reference optimum 19 times. The final follow-up doubled
final sampling to 524288 reads on seeds 74018, 74007, and 74001; it sampled the
reference optimum 28 times total, with seed 74018 producing 20 of those reads.
A current metadata rerun of the same three seeds at the same final-read budget
produced 4 reference reads total (`74018: 0`, `74007: 1`, `74001: 3`). That rerun
is retained as current optimized-parameter and simulator metadata provenance,
not as a replacement for the historical cached final follow-up.

Thus, VQE can reach the reference optimum on the reduced surrogate, but its reference-hit rate is much lower than the transferred-angle QAOA run.

For direct noisy-simulator feasibility, the audit builds the direct full-QUBO
p=2 QAOA circuit from the archived `scalars.csv`, `L_vector.csv`, and
`Q_matrix.csv`. The logical circuit has 36 qubits and 202 `rzz` gates. The
gated FakeFez transpile produces depth 1344 and 1157 `cz` gates. Exact dense
noisy density-matrix simulation would require `2^72` complex entries, while Aer
automatic or MPS-style noisy simulation is entanglement-dependent. No direct
36-qubit noisy samples are cached because a useful-shot direct noisy run was
not established as practical for this circuit. The audit now persists the same
p=2 warm-start beta-then-gamma QAOA parameter vector in
`ds_mfg_direct_full_qubo_audit/direct_full_qubo_qaoa_hardware_parameters.json`.
The optimized reference-hit vector is cached in
`ds_mfg_direct_full_qubo_qaoa_highread/direct_full_qubo_qaoa_optimized_parameters.json`;
`scripts/run_direct_full_qubo_hardware_pilot.jl` prefers that artifact when
present and builds a dry-run-safe direct full-QUBO IBM Runtime handoff, with real
submission kept behind `DSMFG_RUN_DIRECT_FULL_QUBO_HARDWARE=true`.

## Classical Sampling Baselines

Classical baselines were added so the QAOA/VQE hit rates are not interpreted
only against Gurobi. The script `scripts/run_classical_baselines.jl` recomputes
exact auxiliary repair from the original 36-variable QUBO archive, samples the
19 flow variables, and writes cached outputs under
`ds_mfg_classical_baselines/`. The cache includes a summary CSV, a compressed
full aggregate distribution CSV, and a smaller retained-flow CSV for quick
inspection of top-50 and best-sampled flows.

The cached baselines are:

- uniform random sampling with seed 81001 and 262144 samples: 24 top-50 hits, 5
  top-10 hits, 0 reference hits, and best repaired objective `11.8105`,
- uniform random sampling with seed 81002 and 524288 samples: 49 top-50 hits,
  11 top-10 hits, 0 reference hits, and best repaired objective `12.2415`,
- steepest-descent hill climbing with random restarts, seed 82001, and 262144
  repaired-objective evaluations: 883 top-50 hits, 177 top-10 hits, 19 reference
  hits, and best repaired objective `11.7095`.

Against these baselines, top-10-targeted transferred p=5 QAOA remains much more
concentrated than uniform random sampling and the simple hill-climb baseline at
the same 262144 count: 62597 top-50 reads, 14326 top-10 reads, and 1007 reference
reads. The historical selected-seed VQE follow-up seed 74018 is also stronger than uniform
random sampling at 524288 reads, with 391 top-50 reads, 152 top-10 reads, and
20 reference reads.

## Hit-rate uncertainty and time to solution

The cached QAOA, reduced-surrogate VQE, and classical-baseline summary CSVs
report empirical hit rates for the top-50, top-10, reference-optimum, and
stored-Gurobi-pool events. The stored-Gurobi-pool event counts sampled flows
whose repaired projection matches either `global_optimum` or `gurobi_pool`.
Each rate is paired with a two-sided 95% Wilson score confidence interval. The
same rows also
report empirical 99% time to solution (`tts99_sec`), computed as the recorded
local solve or wall time per sample multiplied by
`ceil(log(1 - 0.99) / log(1 - hit_rate))`. Rows with zero observed hits have
infinite empirical time to solution rather than an extrapolated finite value.
For QAOA and VQE, the timing basis is the local Aer `solve_time_sec`; for the
classical baselines, it is the script `wall_time_sec`.

Representative reference-optimum and stored-Gurobi-pool-hit rows are:

| Run | Optimal hits | Optimal `tts99_sec` | Gurobi-pool feasible hits | Gurobi-pool feasible `tts99_sec` |
| --- | ---: | ---: | ---: | ---: |
| Top-10-targeted p=5 QAOA, high-read | 1007 / 262144 | 0.246 | 54820 / 262144 | 0.004 |
| Reduced-surrogate VQE, seed 74018 | 20 / 524288 | 6.926 | 389 / 524288 | 0.356 |
| Hill-climb restarts, seed 82001 | 19 / 262144 | 0.226 | 721 / 262144 | 0.006 |
| Uniform random, seed 81001 | 0 / 262144 | Inf | 21 / 262144 | 0.456 |

The consolidated report is
`ds_mfg_hit_rate_reports/time_to_solution_report.csv`; the per-run summary CSVs
carry the corresponding top-50, top-10, reference, and
`gurobi_pool_feasible_*` columns.

The comparison does not imply a runtime advantage for QAOA or VQE. It shows
that sampled-distribution quality should be interpreted relative to simple
classical sampling baselines once exact repair is available.

## Practical Workflow Burden

The Gurobi baseline is solved and interpreted directly in the original 19 flow
variables. The quantum-emulation workflow required a longer manual path:

1. import the 36-variable QUBO,
2. separate flow and auxiliary bits,
3. repair auxiliary assignments exactly,
4. enumerate all `2^19` repaired flow values,
5. fit a quadratic surrogate because the exact reduced objective is not
   generally quadratic,
6. configure Aer MPS simulator options,
7. separate optimizer reads from final sampling reads,
8. run offline QAOA angle search in `JuliQAOA.jl`,
9. transfer fixed parameters back into `QiskitOpt.jl`,
10. postprocess each sampled state against the repaired objective.

This burden is part of the result. The local QAOA and VQE experiments show that
meaningful samples can be obtained, but not through a drop-in replacement for
the classical solver.

The workflow-effort rubric used for the manuscript is:

| Workflow step | Effort class | Reproducibility status | Interpretation |
| --- | --- | --- | --- |
| QUBO reformulation | Manual interpretation of supplied artifact | Required input archive and CSV files are tracked | The 36-bit QUBO is not the same object as the 19-flow Gurobi model, so the flow/auxiliary split must be documented before sampling results are meaningful. |
| Penalty and scaling interpretation | Manual review plus scripted loading | `scalars.csv`, `L_vector.csv`, and `Q_matrix.csv` are loaded by the notebook and helper scripts | Objective values are only comparable after applying the recorded scale and offset and after distinguishing raw QUBO energy from repaired flow objective. |
| Auxiliary repair | Scripted once the split is identified | The repair components and exact repaired flow scores are cached under `ds_mfg_reduced_flow_objective/` | Exact repair is required for coherent comparison to the original Gurobi flow objective. |
| Reduced surrogate construction | Scripted, problem-specific | The notebook enumerates all `2^19` repaired flow values and writes reduced surrogate CSVs | This enables QAOA/VQE interfaces that require quadratic QUBOs, but the surrogate remains an approximation scored against exact repair. |
| Angle search | Optional scripted/manual tuning | `scripts/find_reduced_qaoa_angles_juliqaoa.jl` expects `JuliQAOA.jl` via a sibling clone or `JULIQAOA_PROJECT` | The highest-concentration QAOA result required offline statevector angle search and parameter transfer rather than a default optimizer run. |
| Backend selection | Manual choice with scripted settings | Aer MPS, FakeFez/Aer, and IBM Runtime paths are separate scripts or configuration gates | Backend choice changes the evidentiary tier: local emulation, model-based noisy simulation, or hardware feasibility. |
| Transpilation | Scripted, backend-specific | Simulator and hardware scripts set or record transpile seeds | Hardware-comparable circuits require explicit transpilation provenance. |
| IBM credential handling | Manual, hardware-specific, optional | `scripts/run_ibm_qaoa_pilot.jl` dry-runs by default and submits only when `DSMFG_RUN_IBM_HARDWARE=true` | IBM tokens and account files must remain outside the repository; dry-run artifacts should not expose secrets. |
| Postprocessing and exact repair | Scripted and required for reported metrics | Scoring scripts and cached summaries report top-50, top-10, reference, feasible-pool, Wilson interval, and time-to-solution fields | The reported result is the repaired-flow distribution, not merely raw sampled QUBO bitstrings. |

## Hardware Pilot Results

The QAOA hardware candidate was executed as a fixed-parameter sampling pilot.
The run used the reduced 19-flow surrogate, p=5, and the beta-then-gamma
top-10-targeted parameter vector learned by `JuliQAOA.jl`. It submitted the
fixed circuit to `ibm_fez` through IBM Runtime, using three repeats and
transpile seeds `92001`, `92002`, and `92003`; each repeat/seed pair used 4096
shots, for 9 submitted jobs and 36864 hardware reads.

The hardware run returned 34573 unique flow assignments. Exact auxiliary repair
scored 6 top-50 repaired-flow hits, 1 top-10 hit, 0 reference-optimum hits, and 4
stored-Gurobi-pool hits. The best repaired flow was rank 2 with objective
`11.8105` and flow bitstring `1001110100111010011`. The hardware artifacts are
stored in `ds_mfg_ibm_qaoa_pilot_fez_4096x3/`.

The cached simulator-to-hardware comparison is stored in
`ds_mfg_simulator_hardware_comparison/`.

| Evidence | Budget or status | Top-50 hits | Top-10 hits | Global hits |
| --- | ---: | ---: | ---: | ---: |
| Ideal Aer transferred-angle QAOA | 262144 reads | 62597 (0.238789) | 14326 (0.0546494) | 1007 (0.0038414) |
| FakeFez/Aer model-based noisy simulation | 36864 reads | 516 (0.0139974) | 108 (0.00292969) | 12 (0.000325521) |
| IBM `ibm_fez` hardware pilot | 36864 reads | 6 (0.000162760) | 1 (0.0000271267) | 0 (0) |

A direct full-QUBO p=2 `ibm_fez` hardware pilot using the optimized high-read
parameters was also run as 1 job with 4096 shots. It found 0 top-50, 0 top-10,
and 0 reference-optimum repaired-flow hits; the best repaired objective was `91.7295`.
The artifacts are in `ds_mfg_direct_full_qubo_hardware_pilot/`. This cached run
used the original hand-built hardware-circuit path and is retained as legacy
descriptive evidence; rerun with the current `QiskitOpt.QAOA.fixed_parameter_circuit`
handoff before treating it as a validated comparison to the high-read noiseless
reference.

Relative to the ideal Aer rates, the model-based noisy simulation retained about
5.9% of the top-50 hit rate, 5.4% of the top-10 hit rate, and 8.5% of the
reference-optimum hit rate. The hardware pilot retained about 0.068% of the top-50
hit rate and about 0.050% of the top-10 hit rate, with no observed
reference-optimum reads. The model-based noisy simulation was generated by
`scripts/run_noisy_qaoa_fake_backend.jl`, using
`qiskit_ibm_runtime.fake_provider.FakeFez` and
`qiskit_aer.AerSimulator.from_backend` with hardware-comparable default shots.
The cached 9-job run took 9949.7 seconds locally on transpiled circuits with
depths 1442--1557 and 1136--1156 `cz` gates. Any FakeFez/Aer result should be
labeled as model-based simulation, not as a calibrated prediction of the
hardware run.

For VQE, the current cache identifies the best configuration rather than a persisted optimized parameter vector. The best final follow-up seed was 74018, using EfficientSU2, 128 optimizer reads, 25 COBYLA iterations, and high final sampling. A VQE hardware follow-up should persist the optimized `result.x` ansatz vector before submission so optimization quality can be separated from hardware sampling noise.

The QAOA hardware pilot is a feasibility measurement of the fixed-circuit
handoff and returned count scoring. It is not a hardware-performance or speedup
claim.

## Final Artifact Audit

The manuscript claims are supported by tracked artifacts rather than by
uncommitted reruns. The principal evidence paths are:

- original QUBO/Gurobi archive: `Fw_ DS mfg case qubo information.zip`,
- Gurobi-pool provenance, exact 19-flow IP enumeration, and local Gurobi rerun:
  `ds_mfg_gurobi_provenance/`,
- exact repair and reduced objective: `ds_mfg_reduced_flow_objective/`,
- QAOA and VQE sampled distributions: `ds_mfg_*summary.csv` files listed in
  the workflow table and checked by `scripts/smoke_test.jl`,
- direct full-QUBO audit and high-read QAOA evidence:
  `ds_mfg_direct_full_qubo_audit/` and
  `ds_mfg_direct_full_qubo_qaoa_highread/`,
- selected-seed VQE optimized-parameter metadata:
  `ds_mfg_vqe_reduced_flow_objective_seed74018_metadata/` and
  `ds_mfg_vqe_reduced_flow_objective_selected_metadata/`,
- classical context: `ds_mfg_classical_baselines/` and
  `ds_mfg_hit_rate_reports/time_to_solution_report.csv`,
- hardware and simulator comparison:
  `ds_mfg_ibm_qaoa_pilot_fez_4096x3/`,
  `ds_mfg_direct_full_qubo_hardware_pilot/`,
  `ds_mfg_fake_fez_qaoa_noisy_4096x3/`, and
  `ds_mfg_simulator_hardware_comparison/`.

The committed hardware and simulator manifests are sanitized provenance
records. They retain public backend names, job IDs, shot counts, transpile
seeds, parameter-artifact paths, and scoring paths, but they do not include IBM
tokens, account files, Runtime instance CRNs, or absolute local output
directories. Real submissions remain gated by explicit environment variables.

The remaining low-energy-weighted surrogate study is future work, not a
blocking result for this manuscript state. It would test a different surrogate
design question; the current manuscript instead reports the uniform quadratic
surrogate, exact repair scoring, direct full-QUBO evidence, classical
baselines, noisy-model comparison, and small hardware pilots.

The executed notebook still contains historical machine-local paths in captured
outputs. Those paths are descriptive notebook provenance, not secrets or
submission artifacts; sanitizing them is deferred to a future notebook rerun.
For FOCAPO-CPC/PSE Press submission, `MANUSCRIPT_FINDINGS.tex` is body content
intended to be inserted after `\psemaketitle` in the journal template wrapper.

## Interpretation

The main methodological finding is that model interpretation matters as much as sampler behavior. Raw QUBO energies alone can be misleading because auxiliary-variable assignments may be inconsistent even when the projected flow is meaningful. Exact auxiliary repair makes the comparison to the original IP objective coherent.

The main algorithmic finding is that reduced-surrogate QAOA with offline
statevector angle search is the highest-concentration approach tested here, and the
objective-targeted top-10 p=5 search improves the earlier energy-targeted
endpoint. VQE benefits from final-sampling separation and multiple
initializations, but it remains less concentrated on the reference optimum. The
best VQE results are useful as evidence that the reduced surrogate contains the
optimum in its sampled distribution, whereas the QAOA transfer result provided
the fixed-parameter circuit used for the hardware pilot.

The classical baselines add sampling context. Uniform random repaired-flow
sampling almost never reaches the reference optimum at these budgets, so both the
QAOA transfer result and the historical selected-seed VQE follow-up are meaningfully more
concentrated than random sampling. A simple hill-climb baseline can also reach
the reference optimum once it has direct access to the repaired objective, which
underscores that exact repair is a strong classical tool as well as a scoring
method.

The practical conclusion is more conservative: the DS-MFG instance can be
reformulated and sampled with current QAOA/VQE software, and local simulations
can recover the Gurobi optimum after reduction and repair, but the workflow
burden is substantial. The current evidence supports feasibility,
interpretability, and engineering-effort claims; it does not support replacing
Gurobi for this process-design instance.

## Limitations And Follow-Up Work

Most quantum-algorithm results here are Aer simulations. The included IBM
hardware pilot is small and should not be interpreted as quantum speedup or IBM
hardware-performance superiority. The reduced surrogate is quadratic and
approximate, even though the exact repaired flow objective is known by
enumeration. Future work should test interfaces that can optimize the exact
reduced objective directly or fit a surrogate weighted toward low-energy flow
assignments. It would also be valuable to persist optimized VQE parameters, not
just initial seeds and sampled distributions, so hardware follow-up can
separate optimization quality from hardware sampling noise.
