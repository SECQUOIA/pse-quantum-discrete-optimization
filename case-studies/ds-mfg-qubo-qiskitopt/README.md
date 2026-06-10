# DS-MFG QUBO Case Study Bundle

This bundle contains a reproducible DS-MFG QUBO case study using local classical
emulations of QAOA and VQE through `QiskitOpt.jl` and Qiskit Aer. It includes
the executed notebook, cached distributions, generated figures, and a
manuscript-ready findings draft.

Start with:

- `README_STUDENT.md`: detailed guide to the notebook, environment, cached
  results, and rerun commands.
- `ds_mfg_qubo_qiskitopt.ipynb`: executed notebook with the full workflow.
- `MANUSCRIPT_FINDINGS.tex`: LaTeX manuscript section with problem formulation,
  reduction, experiment provenance, result tables, and hardware placeholders.
- `MANUSCRIPT_FINDINGS.md`: compact narrative summary.

The main modeling point is that the QUBO has 36 binary variables, but Gurobi's
ground truth is defined over the original 19 flow variables. The remaining 17
bits are auxiliary/slack variables. A sampled QUBO state can therefore project
to a good Gurobi flow while still carrying an inconsistent auxiliary assignment.
The notebook scores samples using exact auxiliary repair, which fixes the flow
bits and chooses the minimum-energy auxiliary completion.

Best cached results:

- Gurobi global optimum: objective `11.7095`, flow
  `1001110100111100011`.
- QAOA: transferred p=5 JuliQAOA angles on the reduced 19-flow surrogate sampled
  the global optimum 664 times in 262144 local Aer reads.
- VQE: the final reduced-surrogate follow-up sampled the global optimum 28
  times across three 524288-read local Aer runs; seed 74018 produced 20 of
  those reads.

No IBM hardware results are included. The notebook contains disabled hardware
placeholder cells for the best QAOA and VQE local configurations.
