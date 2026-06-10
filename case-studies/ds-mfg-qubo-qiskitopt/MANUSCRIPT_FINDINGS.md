# DS-MFG QUBO Case Study: Findings Narrative

For manuscript incorporation, use `MANUSCRIPT_FINDINGS.tex`. That file is a
LaTeX-ready section with experiment provenance, cached artifact paths, result
tables, and hardware-run placeholders. This Markdown file is kept as a compact
companion summary.

## Purpose

This case study evaluates local classical emulations of QAOA and VQE, accessed
through `QiskitOpt.jl`, on a DS-MFG discrete integer-programming QUBO instance.
The analysis compares sampled quantum-optimization distributions against a
Gurobi solution pool and uses the original IP variables to interpret the sampled
QUBO bitstrings.

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

The Gurobi global optimum has repaired objective `11.7095` and flow bitstring
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
objective recovers the Gurobi global optimum. Because the current QAOA/VQE
interfaces accept quadratic QUBOs, the exact reduced objective is approximated
by a least-squares quadratic surrogate over all flow assignments. The surrogate
has high global fit quality (`R^2 = 0.9990819`), but its deterministic minimum is
a near-optimal Gurobi-pool flow with repaired objective `14.5505`, not the true
global optimum. This mismatch explains why surrogate-based quantum sampling must
still be evaluated against the exact repaired objective.

## QAOA Findings

The baseline full-QUBO QAOA and fixed-angle QAOA sweeps did not reliably produce
useful encoded solutions. Fixed-angle QAOA with p=2 produced one projected
Gurobi-pool read in 128 shots, but the auxiliary bits were not encoded
correctly. Exact auxiliary repair scored that projected flow at `18.0745`.

Separating optimizer reads from final sampling improved the full 36-variable
QAOA run. A p=2 run using 128 optimizer reads and 512 final reads produced three
projected pool reads and improved the best repaired pool objective to `14.5505`.

The strongest QAOA result came from moving angle search off the Qiskit/Aer optimizer path. The workflow was:

1. eliminate auxiliary variables exactly to build the reduced 19-flow objective,
2. fit a quadratic surrogate over the exact reduced objective,
3. use `JuliQAOA.jl` as a local statevector angle finder on the reduced surrogate,
4. transfer the learned p=5 QAOA angles into `QiskitOpt.QAOA`,
5. set `MaximumIterations() = 0` and use Aer only for final sampling.

This transferred-angle QAOA run sampled the Gurobi global optimum 664 times in
262144 final reads. It is the strongest cached quantum-emulation result in the
notebook. Its best repaired objective is `11.7095`, matching the Gurobi optimum,
and the sampled optimal flow is `1001110100111100011`.

## VQE Findings

The baseline VQE run was weak: it collapsed to a small set of sampled states and
did not provide a useful Gurobi-pool hit. The VQE experiments were then
strengthened in two ways:

- optimizer reads were separated from final sampling reads through `QUBODrivers.FinalNumberOfReads()`,
- several random EfficientSU2 initializations were tested.

On the full 36-variable QUBO, a 5-seed VQE sweep with 8192 final reads found a
projected Gurobi-pool flow with repaired objective `14.6515`. On the reduced
19-flow surrogate, a 5-seed sweep improved the best repaired objective to
`11.8105`, close to the global optimum.

After package-level Aer MPS backend attributes were available, the
reduced-surrogate VQE study was repeated with more seeds and higher final
sampling. A 20-seed sweep with 32768 final reads per seed sampled the Gurobi
global optimum 3 times. A high-read follow-up on selected seeds with 262144
final reads sampled the global optimum 19 times. The final follow-up doubled
final sampling to 524288 reads on seeds 74018, 74007, and 74001; it sampled the
global optimum 28 times total, with seed 74018 producing 20 of those reads.

Thus, VQE can reach the global optimum on the reduced surrogate, but its global-hit rate is much lower than the transferred-angle QAOA run.

## Hardware Run Plan

The notebook includes placeholders for hardware execution using the best local findings.

For QAOA, the hardware candidate is well-defined: use the reduced 19-flow surrogate, p=5, and the beta-then-gamma parameter vector learned by `JuliQAOA.jl`. This is a fixed-angle sampling run with `MaximumIterations() = 0`, so the hardware execution focuses on sampling the transferred p=5 circuit.

For VQE, the current cache identifies the best configuration rather than a persisted optimized parameter vector. The best final follow-up seed was 74018, using EfficientSU2, 128 optimizer reads, 25 COBYLA iterations, and high final sampling. The hardware placeholder therefore reruns VQE from the best initial seed unless the scripts are extended to persist the optimized `result.x` ansatz vector.

## Interpretation

The main methodological finding is that model interpretation matters as much as sampler behavior. Raw QUBO energies alone can be misleading because auxiliary-variable assignments may be inconsistent even when the projected flow is meaningful. Exact auxiliary repair makes the comparison to the original IP objective coherent.

The main algorithmic finding is that reduced-surrogate QAOA with offline statevector angle search is the strongest approach tested here. VQE benefits from final-sampling separation and multiple initializations, but it remains less concentrated on the global optimum. The best VQE results are useful as evidence that the reduced surrogate contains the optimum in its sampled distribution, whereas the QAOA transfer result provides the clearest candidate for hardware sampling.

## Limitations And Follow-Up Work

The local results are Aer simulations, not hardware executions. The reduced surrogate is quadratic and approximate, even though the exact repaired flow objective is known by enumeration. Future work should test interfaces that can optimize the exact reduced objective directly or fit a surrogate weighted toward low-energy flow assignments. It would also be valuable to persist optimized VQE parameters, not just initial seeds and sampled distributions, so hardware follow-up can separate optimization quality from hardware sampling noise.
