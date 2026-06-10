# PSE Quantum Discrete Optimization

This repository collects process systems engineering case studies that use
quantum and quantum-inspired algorithms for discrete optimization models.

The initial case study is a demand-side manufacturing (DS-MFG) QUBO workflow
using local Qiskit Aer emulations of QAOA and VQE through `QiskitOpt.jl`.
The repository layout is intentionally case-study based so that related PSE
instances can be added without changing the existing DS-MFG artifact paths.

## Case Studies

- `case-studies/ds-mfg-qubo-qiskitopt`: executed notebook, cached
  distributions, generated figures, scripts, and manuscript-ready findings for
  the DS-MFG QUBO case.

## Start Here

For the current DS-MFG case:

```bash
cd case-studies/ds-mfg-qubo-qiskitopt
```

Then read:

1. `README_STUDENT.md` for the detailed guide to the notebook, environment,
   cached results, and rerun commands.
2. `MANUSCRIPT_FINDINGS.tex` for the manuscript-ready formulation and results.
3. `ds_mfg_qubo_qiskitopt.ipynb` for the executed analysis record.

All quantum-algorithm results currently tracked here are local classical
emulations. No IBM hardware results are included in the initial DS-MFG bundle.
