#!/usr/bin/env julia

# Lightweight reproducibility check for a fresh clone. This intentionally
# validates imports and cached artifacts without rerunning the expensive
# notebook simulations.

using PythonCall
using QUBO
using QUBODrivers
using QUBOTools
using QiskitOpt

const ROOT = normpath(joinpath(@__DIR__, ".."))
const GLOBAL_FLOW_BITS = "1001110100111100011"
const GLOBAL_OBJECTIVE = 11.7095

cd(ROOT)

function smoke_error(message)
    error("DS-MFG smoke test failed: $message")
end

function require_file(relpath)
    path = joinpath(ROOT, relpath)
    isfile(path) || smoke_error("missing required file: $relpath")
    filesize(path) > 0 || smoke_error("required file is empty: $relpath")
    return path
end

function require_dir(relpath)
    path = joinpath(ROOT, relpath)
    isdir(path) || smoke_error("missing required directory: $relpath")
    return path
end

function parse_csv_rows(relpath)
    path = require_file(relpath)
    lines = filter(line -> !isempty(strip(line)), readlines(path))
    length(lines) >= 2 || smoke_error("expected a header and at least one row in $relpath")

    header = split(lines[1], ","; keepempty = true)
    rows = Dict{String,String}[]
    for (line_number, line) in enumerate(lines[2:end])
        fields = split(line, ","; keepempty = true)
        length(fields) == length(header) || smoke_error(
            "row $(line_number + 1) in $relpath has $(length(fields)) fields; expected $(length(header))",
        )
        push!(rows, Dict(zip(header, fields)))
    end
    return rows
end

function parse_key_value_csv(relpath)
    rows = parse_csv_rows(relpath)
    all(haskey(row, "key") && haskey(row, "value") for row in rows) ||
        smoke_error("$relpath must contain key,value columns")
    return Dict(row["key"] => row["value"] for row in rows)
end

function require_value(row, key, expected)
    actual = get(row, key, nothing)
    actual == expected || smoke_error("expected $key=$expected, got $actual")
end

function require_int(row, key, expected)
    actual = parse(Int, row[key])
    actual == expected || smoke_error("expected $key=$expected, got $actual")
end

function require_float(row, key, expected; atol = 1e-8)
    actual = parse(Float64, row[key])
    isapprox(actual, expected; atol = atol, rtol = 0.0) ||
        smoke_error("expected $key to be approximately $expected, got $actual")
end

println("Checking Julia and Python package imports...")
Base.pkgversion(QiskitOpt) == v"0.4.2" ||
    smoke_error("expected QiskitOpt v0.4.2, got $(Base.pkgversion(QiskitOpt))")
PythonCall.pyimport("qiskit")
PythonCall.pyimport("qiskit_aer")
PythonCall.pyimport("qiskit_optimization")

println("Checking required case-study files...")
for relpath in (
    "Project.toml",
    "Manifest.toml",
    "CondaPkg.toml",
    "README_STUDENT.md",
    "MANUSCRIPT_FINDINGS.tex",
    "MANUSCRIPT_FINDINGS.md",
    "ds_mfg_qubo_qiskitopt.ipynb",
    "Fw_ DS mfg case qubo information.zip",
    "ds_mfg_qaoa_vqe_distribution.svg",
)
    require_file(relpath)
end

for relpath in (
    "ds_mfg_saved_distributions",
    "ds_mfg_qaoa_juliqaoa_transfer_highread",
    "ds_mfg_vqe_reduced_flow_objective_final",
    "ds_mfg_reduced_flow_objective",
    "scripts",
)
    require_dir(relpath)
end

println("Validating cached result summaries...")
metadata = parse_key_value_csv("ds_mfg_saved_distributions/run_metadata.csv")
isapprox(parse(Float64, metadata["gurobi_best_objective"]), GLOBAL_OBJECTIVE; atol = 1e-8, rtol = 0.0) ||
    smoke_error("unexpected Gurobi best objective in run metadata")
metadata["qubo_dimension"] == "36" || smoke_error("unexpected QUBO dimension in run metadata")

reduced = only(parse_csv_rows("ds_mfg_reduced_flow_objective/reduced_flow_summary.csv"))
require_value(reduced, "n_flow_variables", "19")
require_value(reduced, "n_auxiliary_variables_eliminated", "17")
require_value(reduced, "exact_global_flow_bits", GLOBAL_FLOW_BITS)
require_float(reduced, "exact_global_reduced_energy", GLOBAL_OBJECTIVE)

qaoa = only(parse_csv_rows("ds_mfg_qaoa_juliqaoa_transfer_highread/qaoa_juliqaoa_transfer_summary.csv"))
require_value(qaoa, "algorithm", "QAOA_reduced_surrogate_JuliQAOA_transfer")
require_int(qaoa, "p", 5)
require_int(qaoa, "total_reads", 262144)
require_int(qaoa, "global_hits", 664)
require_value(qaoa, "best_top50_match", "global_optimum")
require_value(qaoa, "best_top50_flow_bits", GLOBAL_FLOW_BITS)

vqe_rows = parse_csv_rows("ds_mfg_vqe_reduced_flow_objective_final/vqe_reduced_top50_sampling_summary.csv")
length(vqe_rows) == 3 || smoke_error("expected three final VQE follow-up rows")
sum(parse(Int, row["global_hits"]) for row in vqe_rows) == 28 ||
    smoke_error("unexpected final VQE global-hit total")
all(row["best_top50_match"] == "global_optimum" for row in vqe_rows) ||
    smoke_error("not every final VQE row reaches the global optimum")
all(row["best_top50_flow_bits"] == GLOBAL_FLOW_BITS for row in vqe_rows) ||
    smoke_error("unexpected VQE optimum flow bits")

println("DS-MFG smoke test passed.")
