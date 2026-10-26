# Reviewer robustness and computational-cost audit

This revision extends submission commit
`93c4f97e13936c2434c84afe5ebf8b8a4794d722`. The main paper is limited to one
additional page. Detailed configurations and results live here; the original
experiment caches remain unchanged. New measurements use one Julia thread and
one BLAS thread on the CPU recorded in the metadata files.

| Reviewer request | Evidence and scope |
| --- | --- |
| Penalty sensitivity | `penalty_sensitivity.csv`: eight common multipliers of all penalty terms, with capital costs fixed. Every flow is audited exactly. This is not independent tuning of each constraint family. |
| Surrogate construction | `surrogate_sensitivity.csv`: the full uniform fit and uniform subsets of 1024, 8192, and 65536 flows, with three seeds per size. No fit is selected for favorable performance. |
| QAOA depth, optimizer, initialization | `qaoa_sensitivity.csv`: every depth 1–5 for seeds 91001–91003 at one basin-hopping iteration, plus seed 91001 at five iterations. `cached_qaoa_sweeps.csv` inventories the historical native-QAOA depth, iteration-cap, and angle-perturbation runs. |
| VQE selection bias | `all_vqe_seed_runs.csv` inventories 43 runs across ten campaigns, including the original 128-read baseline. `vqe_cohorts.csv` aggregates only within each campaign. The main comparison uses the complete 20-seed screen. |
| Full computational cost | `historical_cost_ledger.csv` explicitly represents unrecorded stages as unknown. New preprocessing, angle-optimization, environment/warmup, and offline compilation measurements are in separate files; they are not added to historical sampling times. |
| Additional classical methods | `classical_runs.csv`: uniform sampling, hill-climb restarts, and simulated annealing, each with five seeds and 262144 evaluations per run. |
| Scalability and references | The main text states the exponential enumeration/statevector and component-repair costs. No larger-instance scalability claim is made. The classical annealing reference is Kirkpatrick, Gelatt, and Vecchi (1983), DOI `10.1126/science.220.4598.671`. |

## Definitions and experimental protocol

Penalty scaling uses `E_alpha(y,z) = c(y) + alpha * P(y,z)`, with the capital
costs from `scripts/enumerate_ip_provenance.py`. Since `alpha > 0` and `c` does
not depend on the auxiliaries, the exact repaired score is
`F_alpha(y) = c(y) + alpha * (F_1(y) - c(y))`. The feasible set is the 36 designs
independently certified in `ds_mfg_gurobi_provenance/`. At the submitted weight,
all feasible designs precede all infeasible designs; the top-50 set contains
infeasible flows only at ranks 37–50. The revised paper corrects the earlier
statement that these flows interleave with feasible ones.

Surrogates contain an intercept, 19 linear terms, and 171 pair terms. The full
uniform least-squares fit is evaluated through an orthonormal Walsh basis:
transform all repaired values and retain degrees zero, one, and two. This
recovers the submitted polynomial and its optimum/rank failure without a dense
524288-by-191 design matrix. The subset fits use ordinary least squares on
random samples without replacement. Full enumeration is still used for the
audit, so the fast Walsh transform is not evidence of scalable optimization.
`fit_seconds` includes training-index selection, feature construction, and the
linear solve for subset fits; full-fit time covers the transform and coefficient
extraction. It excludes obtaining labels and evaluating the resulting model.

QAOA uses the submitted reduced Hamiltonian, its recorded normalization, the
existing JuliQAOA `find_angles_bh` path, and energy minimization. Each optimizer
timer covers the entire incremental depth-1-to-5 sequence, not just depth 5.
The seeds share the algorithm's deterministic initial depth guesses; randomness
enters its basin-hopping proposals and affects subsequent depth initialization.
One versus five iterations is paired at seed 91001. The three-seed comparison
is at one iteration only, so this is a limited sensitivity study rather than a
full factorial test or a population-level confidence statement. The native
historical COBYLA iteration sweep also changes seeds; it cannot isolate an
iteration-cap effect. Its complete exploratory results are retained as such.

For the classical comparison, each run has its own empty repair cache and
counts every queried candidate against the evaluation budget, including
revisits and rejected moves. Each method uses the same exact scorer and
distribution bookkeeping. Search is run before exhaustive enumeration and has
no access to the enumerated objective, surrogate, optimum, or feasible pool.
The optimum is used only for reporting. Compilation is warmed with discarded
256-evaluation runs. Hill climbing takes steepest improving single-bit moves,
with at most 100 steps per restart. Annealing proposes a random single-bit flip
and accepts uphill changes with probability `exp(-delta/T)`. It cools
geometrically from 100 to 0.01 over each 4096-evaluation restart. Settings are
fixed across seeds; no tuning or best-seed selection is applied. Reported
`search_seconds` includes cache lookup, exact repair on a cache miss, and
candidate recording. It excludes process startup and output writing.

The original baseline's shared cache and historical machine/timing scope differ
from this fresh-cache comparison. Compare the new methods with each other;
do not splice their times into the historical Table 1. Heuristic evaluations
are correlated. We report run successes and full-budget elapsed times, not a
binomial confidence interval or an independent-trial TTS guarantee.

## VQE selection and complete reporting

The full screen uses seeds 74001–74020, 128 optimizer reads, 25 COBYLA
iterations, and 32768 final reads per seed. There are three optimum hits:
seed 74001 has one, seed 74007 has two, and the other 18 seeds have none.
Seed 74018 has the highest top-50 count (23) but no optimum hit in that screen.
Selection of 74001 and 74007 is consistent with their optimum hits; selection
of 74018 is consistent with its top-50 count. There is no retained preregistered
selection rule, so all three follow-ups are described as exploratory.

At 524288 reads per seed, the historical follow-ups have 5, 3, and 20 optimum
hits for seeds 74001, 74007, and 74018, respectively. Later metadata reruns at
that budget have 3, 1, and 0. The metadata reruns preserve optimized parameters
for their own samples, not the historical selected-seed result. The main table
now uses the whole screen (655360 reads), and the supplement links to every
seed/campaign instead of representing the best follow-up as typical.

Pooled Wilson intervals and TTS values in `vqe_cohorts.csv` are descriptive
calculations for the recorded read mixture. They do not capture between-seed
uncertainty or establish a guaranteed repeated-run success probability. The
original baseline's unrecorded solve time remains blank, including its TTS.

## Cost accounting

A complete run must account for environment/setup, IP-to-QUBO reformulation,
exhaustive flow enumeration with auxiliary repair, surrogate construction,
offline angle optimization, logical circuit construction, compilation, sampling
or execution, hardware queueing and retrieval, and final projection/repair and
scoring. Some branches reuse setup; report that reuse and its amortization
explicitly rather than dropping the stages.

The historical local numbers time solver calls, not the complete workflow.
The hardware archive provides a 73.488455-second runner interval and no
per-stage compilation, queue, or execution timestamps. It cannot establish
which component dominates. Original reformulation, enumeration/fit, and
offline optimization timings are absent. Consequently, a complete historical
total cannot be reconstructed and missing costs are not zero. No IBM account
was configured for retrospective job-metric retrieval; no new jobs were sent.

The new measurements are reference reruns, not recovered historical times:

- `preprocessing_times.csv` times input loading/validation and exhaustive
  repair; `surrogate_sensitivity.csv` separately times the fit.
- `qaoa_timing.csv` records four complete optimizer runs and probability
  scoring. `qaoa_metadata.toml` also records environment load, a small optimizer
  warmup, Hamiltonian enumeration/setup, and total script wall time.
- `compilation_proxy.csv` times three transpiler seeds at optimization level 3
  on the bundled FakeFez model, using the actual pilot circuit builder. The
  logical circuit is checked against the archived depth and operation counts.
  `compilation_metadata.toml` records the other setup/build/warmup costs. These
  offline model timings cannot replace historical hardware compilation or
  queue times.

The snapshot's CPU/software versions and timer boundaries matter. Startup,
compilation, and cache effects are stated separately; no timing from a new
implementation is added to the old Aer or hardware interval to manufacture an
end-to-end total.

## Reproduce and check

From `case-studies/ds-mfg-qubo-qiskitopt`, use Julia 1.10 and the project's
documented Python environment. The first script uses only Julia standard
libraries; the evidence summarizer uses only the Python standard library.

```bash
julia --startup-file=no --threads=1 --project=. scripts/run_revision_sensitivity.jl
python scripts/summarize_revision_evidence.py
julia --startup-file=no --threads=1 --project=. scripts/test_revision_sensitivity.jl
python scripts/summarize_revision_evidence.py --check
```

The angle audit uses `lanl/JuliQAOA.jl` at commit
`6b28d23febe092954cae55b95305f9e5060ef374`. Its `Project.toml` and `Manifest.toml`
are archived as `JuliQAOA_Project.toml` and `JuliQAOA_Manifest.toml`; install those
in a checkout of that commit and instantiate it before running:

```bash
julia --startup-file=no --threads=1 --project=/path/to/JuliQAOA.jl scripts/run_revision_qaoa.jl
julia --startup-file=no --threads=1 --project=. scripts/run_revision_compilation.jl
```

Run timing commands sequentially. They overwrite only this revision's result
files; historical caches are unchanged. Rerun times can vary. Tests compare the
repair decomposition and scaled penalties with brute-force enumeration over
all 131072 auxiliary states for feasible and random flows, verify Walsh/OLS
equivalence and agreement with the submitted surrogate, and check search
budgets and deterministic seed reuse. CI runs these tests and regenerates the
cached evidence in a temporary directory to check exact table parity.

`source_checksums.json` pins the source records for the cached-evidence tables.
The measurement metadata records the scripts and input/environment hashes used
for the new runs. The manuscript PDFs are rebuilt with the LuaLaTeX commands in
the parent README; the main PDF must remain at most seven pages.
