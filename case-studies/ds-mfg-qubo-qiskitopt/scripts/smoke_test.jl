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

include(joinpath(@__DIR__, "hit_rate_stats.jl"))

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

function validate_classical_distribution(relpath, summary_rows)
    path = require_file(relpath)
    expected = Dict(
        (row["algorithm"], row["seed"], row["sample_budget"]) => row for row in summary_rows
    )
    counts = Dict{Tuple{String,String,String},Int}()
    visits = Dict{Tuple{String,String,String},Int}()
    probabilities = Dict{Tuple{String,String,String},Float64}()

    open(`gzip -cd $path`, "r") do io
        header = split(readline(io), ","; keepempty = true)
        for (line_number, line) in enumerate(eachline(io))
            isempty(strip(line)) && continue
            fields = split(line, ","; keepempty = true)
            length(fields) == length(header) || smoke_error(
                "row $(line_number + 1) in $relpath has $(length(fields)) fields; expected $(length(header))",
            )
            row = Dict(zip(header, fields))
            key = (row["algorithm"], row["seed"], row["sample_budget"])
            haskey(expected, key) || smoke_error("unexpected classical distribution key: $key")
            counts[key] = get(counts, key, 0) + 1
            visits[key] = get(visits, key, 0) + parse(Int, row["visits"])
            probabilities[key] = get(probabilities, key, 0.0) + parse(Float64, row["probability"])
        end
    end

    for (key, row) in expected
        get(counts, key, 0) == parse(Int, row["unique_flows"]) ||
            smoke_error("classical distribution unique-flow count mismatch for $key")
        get(visits, key, 0) == parse(Int, row["total_samples"]) ||
            smoke_error("classical distribution visit total mismatch for $key")
        isapprox(get(probabilities, key, 0.0), 1.0; atol = 1e-8, rtol = 0.0) ||
            smoke_error("classical distribution probability mass mismatch for $key")
    end
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

function require_hit_rate_stats(row; total_key = "total_reads", time_key = "solve_time_sec")
    total_reads = parse(Int, row[total_key])
    stat_values = hit_rate_stat_values(
        total_reads,
        parse(Int, row["top50_hits"]),
        parse(Int, row["top10_hits"]),
        parse(Int, row["global_hits"]),
        parse(Int, row["gurobi_pool_feasible_hits"]),
        parse(Float64, row[time_key]),
    )
    for (column, expected) in zip(hit_rate_stat_headers(), stat_values)
        require_value(row, column, expected)
    end
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
    "ds_mfg_classical_baselines",
    "ds_mfg_hit_rate_reports",
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
require_int(qaoa, "gurobi_pool_feasible_hits", 32502)
require_value(qaoa, "best_top50_match", "global_optimum")
require_value(qaoa, "best_top50_flow_bits", GLOBAL_FLOW_BITS)
require_hit_rate_stats(qaoa)

vqe_rows = parse_csv_rows("ds_mfg_vqe_reduced_flow_objective_final/vqe_reduced_top50_sampling_summary.csv")
length(vqe_rows) == 3 || smoke_error("expected three final VQE follow-up rows")
sum(parse(Int, row["global_hits"]) for row in vqe_rows) == 28 ||
    smoke_error("unexpected final VQE global-hit total")
sum(parse(Int, row["gurobi_pool_feasible_hits"]) for row in vqe_rows) == 601 ||
    smoke_error("unexpected final VQE Gurobi-pool feasible-hit total")
all(row["best_top50_match"] == "global_optimum" for row in vqe_rows) ||
    smoke_error("not every final VQE row reaches the global optimum")
all(row["best_top50_flow_bits"] == GLOBAL_FLOW_BITS for row in vqe_rows) ||
    smoke_error("unexpected VQE optimum flow bits")
for row in vqe_rows
    require_hit_rate_stats(row)
end

classical_rows = parse_csv_rows("ds_mfg_classical_baselines/classical_baseline_summary.csv")
length(classical_rows) == 3 || smoke_error("expected three classical baseline rows")
uniform_262 = only(filter(
    row -> row["algorithm"] == "uniform_random_repaired_flow" && row["sample_budget"] == "262144",
    classical_rows,
))
require_int(uniform_262, "top50_hits", 24)
require_int(uniform_262, "top10_hits", 5)
require_int(uniform_262, "global_hits", 0)
require_int(uniform_262, "gurobi_pool_feasible_hits", 21)
require_value(uniform_262, "best_top50_rank", "2")
require_hit_rate_stats(uniform_262; total_key = "total_samples", time_key = "wall_time_sec")

uniform_524 = only(filter(
    row -> row["algorithm"] == "uniform_random_repaired_flow" && row["sample_budget"] == "524288",
    classical_rows,
))
require_int(uniform_524, "top50_hits", 49)
require_int(uniform_524, "top10_hits", 11)
require_int(uniform_524, "global_hits", 0)
require_int(uniform_524, "gurobi_pool_feasible_hits", 33)
require_value(uniform_524, "best_top50_rank", "3")
require_hit_rate_stats(uniform_524; total_key = "total_samples", time_key = "wall_time_sec")

hill = only(filter(row -> row["algorithm"] == "hill_climb_restarts_repaired_flow", classical_rows))
require_int(hill, "top50_hits", 883)
require_int(hill, "top10_hits", 177)
require_int(hill, "global_hits", 19)
require_int(hill, "gurobi_pool_feasible_hits", 721)
require_value(hill, "best_match", "global_optimum")
require_value(hill, "best_flow_bits", GLOBAL_FLOW_BITS)
require_hit_rate_stats(hill; total_key = "total_samples", time_key = "wall_time_sec")

tts_rows = parse_csv_rows("ds_mfg_hit_rate_reports/time_to_solution_report.csv")
all(haskey(row, "total_trials") for row in tts_rows) ||
    smoke_error("time-to-solution report must use total_trials column")
qaoa_global_tts = only(filter(
    row -> row["source_summary"] == "ds_mfg_qaoa_juliqaoa_transfer_highread/qaoa_juliqaoa_transfer_summary.csv" &&
        row["event"] == "global",
    tts_rows,
))
require_int(qaoa_global_tts, "hits", 664)
require_value(qaoa_global_tts, "hit_rate", qaoa["global_hit_rate"])
require_value(qaoa_global_tts, "tts99_sec", qaoa["global_tts99_sec"])

qaoa_feasible_tts = only(filter(
    row -> row["source_summary"] == "ds_mfg_qaoa_juliqaoa_transfer_highread/qaoa_juliqaoa_transfer_summary.csv" &&
        row["event"] == "gurobi_pool_feasible",
    tts_rows,
))
require_int(qaoa_feasible_tts, "hits", 32502)
require_value(qaoa_feasible_tts, "hit_rate", qaoa["gurobi_pool_feasible_hit_rate"])
require_value(qaoa_feasible_tts, "tts99_sec", qaoa["gurobi_pool_feasible_tts99_sec"])

uniform_global_tts = only(filter(
    row -> row["source_summary"] == "ds_mfg_classical_baselines/classical_baseline_summary.csv" &&
        row["algorithm"] == "uniform_random_repaired_flow" &&
        occursin("sample_budget=262144", row["run_label"]) &&
        row["event"] == "global",
    tts_rows,
))
require_value(uniform_global_tts, "tts99_sec", "Inf")

validate_classical_distribution(
    "ds_mfg_classical_baselines/classical_baseline_distribution.csv.gz",
    classical_rows,
)
retained_rows = parse_csv_rows("ds_mfg_classical_baselines/classical_baseline_retained_flows.csv")
all(!isempty(get(row, "retained_reason", "")) for row in retained_rows) ||
    smoke_error("every retained classical baseline row must include a retained_reason")

println("DS-MFG smoke test passed.")
