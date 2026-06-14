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
    "ds_mfg_qaoa_juliqaoa_objective_angle_search",
    "ds_mfg_qaoa_juliqaoa_objective_transfer_highread",
    "ds_mfg_qaoa_juliqaoa_global_angle_search_p3",
    "ds_mfg_qaoa_juliqaoa_global_transfer_p3",
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

objective_angle_rows = parse_csv_rows("ds_mfg_qaoa_juliqaoa_objective_angle_search/juliqaoa_angle_summary.csv")
length(objective_angle_rows) == 10 || smoke_error("expected ten objective-targeted QAOA angle rows")
top50_angle = only(filter(row -> row["target"] == "top50" && row["p"] == "5", objective_angle_rows))
require_int(top50_angle, "basinhopping_niter", 5)
require_float(top50_angle, "top50_probability", 0.23122699481130704; atol = 1e-12)
require_float(top50_angle, "top10_probability", 0.0522370454529187; atol = 1e-12)
require_float(top50_angle, "global_probability", 0.003392140410350963; atol = 1e-12)
top10_angle = only(filter(row -> row["target"] == "top10" && row["p"] == "5", objective_angle_rows))
require_int(top10_angle, "basinhopping_niter", 5)
require_float(top10_angle, "top50_probability", 0.22967538871326482; atol = 1e-12)
require_float(top10_angle, "top10_probability", 0.05261172053991265; atol = 1e-12)
require_float(top10_angle, "global_probability", 0.003456748192282216; atol = 1e-12)

global_angle_rows = parse_csv_rows("ds_mfg_qaoa_juliqaoa_global_angle_search_p3/juliqaoa_angle_summary.csv")
length(global_angle_rows) == 3 || smoke_error("expected three direct-global QAOA angle rows")
global_angle = only(filter(row -> row["target"] == "global" && row["p"] == "3", global_angle_rows))
require_int(global_angle, "basinhopping_niter", 1)
require_float(global_angle, "top50_probability", 0.003936791098224448; atol = 1e-12)
require_float(global_angle, "top10_probability", 0.0010251847647793182; atol = 1e-12)
require_float(global_angle, "global_probability", 9.831625827252306e-5; atol = 1e-12)

objective_qaoa_rows = parse_csv_rows(
    "ds_mfg_qaoa_juliqaoa_objective_transfer_highread/qaoa_juliqaoa_transfer_summary.csv",
)
length(objective_qaoa_rows) == 2 || smoke_error("expected two objective-targeted QAOA transfer rows")
top50_qaoa = only(filter(row -> row["angle_target"] == "top50" && row["p"] == "5", objective_qaoa_rows))
require_int(top50_qaoa, "total_reads", 262144)
require_int(top50_qaoa, "top50_hits", 62964)
require_int(top50_qaoa, "top10_hits", 14221)
require_int(top50_qaoa, "global_hits", 979)
require_int(top50_qaoa, "gurobi_pool_feasible_hits", 55284)
require_value(top50_qaoa, "best_top50_match", "global_optimum")
require_value(top50_qaoa, "best_top50_flow_bits", GLOBAL_FLOW_BITS)
require_hit_rate_stats(top50_qaoa)
top10_qaoa = only(filter(row -> row["angle_target"] == "top10" && row["p"] == "5", objective_qaoa_rows))
require_int(top10_qaoa, "total_reads", 262144)
require_int(top10_qaoa, "top50_hits", 62597)
require_int(top10_qaoa, "top10_hits", 14326)
require_int(top10_qaoa, "global_hits", 1007)
require_int(top10_qaoa, "gurobi_pool_feasible_hits", 54820)
require_value(top10_qaoa, "best_top50_match", "global_optimum")
require_value(top10_qaoa, "best_top50_flow_bits", GLOBAL_FLOW_BITS)
require_hit_rate_stats(top10_qaoa)

global_qaoa = only(parse_csv_rows("ds_mfg_qaoa_juliqaoa_global_transfer_p3/qaoa_juliqaoa_transfer_summary.csv"))
require_value(global_qaoa, "angle_target", "global")
require_int(global_qaoa, "p", 3)
require_int(global_qaoa, "total_reads", 262144)
require_int(global_qaoa, "top50_hits", 1015)
require_int(global_qaoa, "top10_hits", 312)
require_int(global_qaoa, "global_hits", 45)
require_int(global_qaoa, "gurobi_pool_feasible_hits", 924)
require_value(global_qaoa, "best_top50_match", "global_optimum")
require_value(global_qaoa, "best_top50_flow_bits", GLOBAL_FLOW_BITS)
require_hit_rate_stats(global_qaoa)

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

println("Checking IBM QAOA pilot dry-run schema...")
mktempdir() do pilot_output_dir
    withenv(
        "QISKIT_IBM_BACKEND" => "ibm_brisbane",
        "QISKIT_IBM_INSTANCE" => "",
        "DSMFG_HARDWARE_FINAL_READS" => "64",
        "DSMFG_HARDWARE_REPEATS" => "1",
        "DSMFG_HARDWARE_TRANSPILE_SEEDS" => "123",
        "DSMFG_HARDWARE_OUTPUT_DIR" => pilot_output_dir,
        "DSMFG_RUN_IBM_HARDWARE" => "false",
    ) do
        run(`$(Base.julia_cmd()) --project=$(ROOT) scripts/run_ibm_qaoa_pilot.jl`)
    end

    for filename in ("job_manifest.json", "backend_metadata.json", "raw_counts.csv", "scored_counts.csv", "summary.csv")
        path = joinpath(pilot_output_dir, filename)
        isfile(path) || smoke_error("IBM pilot dry run did not write $(filename)")
        filesize(path) > 0 || smoke_error("IBM pilot dry-run file is empty: $(filename)")
    end

    manifest_text = read(joinpath(pilot_output_dir, "job_manifest.json"), String)
    occursin("\"mode\":\"dry_run\"", manifest_text) ||
        smoke_error("IBM pilot manifest must record dry_run mode")
    occursin("\"channel\":\"ibm_quantum_platform\"", manifest_text) ||
        smoke_error("IBM pilot manifest must record the default Runtime channel")
    occursin("\"submitted\":false", manifest_text) ||
        smoke_error("IBM pilot dry-run manifest must not mark jobs submitted")
    !occursin("QISKIT_IBM_TOKEN", manifest_text) ||
        smoke_error("IBM pilot manifest must not include token environment names")
    !occursin("qiskit-ibm.json", manifest_text) ||
        smoke_error("IBM pilot manifest must not include account file paths")

    summary_lines = filter(line -> !isempty(strip(line)), readlines(joinpath(pilot_output_dir, "summary.csv")))
    length(summary_lines) == 2 || smoke_error("IBM pilot dry-run summary must contain one data row")
    occursin("QAOA_reduced_surrogate_JuliQAOA_IBM_pilot,dry_run,ibm_brisbane,top10,5,64,1,123,0,", summary_lines[2]) ||
        smoke_error("IBM pilot dry-run summary row has unexpected configuration values")
end

println("Checking IBM QAOA pilot manifest durability helpers...")
include(joinpath(@__DIR__, "run_ibm_qaoa_pilot.jl"))
mktempdir() do pilot_output_dir
    config = PilotConfig(
        DEFAULT_IBM_RUNTIME_CHANNEL,
        "ibm_brisbane",
        nothing,
        64,
        1,
        [123],
        pilot_output_dir,
        true,
    )
    jobs = planned_jobs(config)
    data = (
        scalars = (n = 19, scale = 1.0, offset = 0.0),
        angle_record = Dict(
            "seed" => "91001",
            "basinhopping_niter" => "5",
            "top50_probability" => "0.22967538871326482",
            "top10_probability" => "0.05261172053991265",
            "global_probability" => "0.003456748192282216",
        ),
        angle_path = joinpath(ANGLE_DIR, "juliqaoa_angle_summary.csv"),
        angles = fill(0.0, 2 * ANGLE_P),
        top_flows = Dict{String,NamedTuple}(),
    )
    circuit_info = Dict{String,Any}("num_qubits" => 19, "depth" => 0)
    paths = output_paths(config)

    jobs[1]["submitted"] = true
    jobs[1]["job_id"] = "synthetic-runtime-job"
    jobs[1]["status"] = "RUNNING"
    persist_job_manifest!(config, data, circuit_info, jobs, paths)
    manifest_text = read(paths["job_manifest"], String)
    occursin("\"submitted\":true", manifest_text) ||
        smoke_error("IBM pilot manifest must persist submitted job state")
    occursin("\"job_id\":\"synthetic-runtime-job\"", manifest_text) ||
        smoke_error("IBM pilot manifest must persist Runtime job IDs")

    try
        error("synthetic after-submission failure")
    catch err
        record_job_failure!(jobs[1], nothing, true, err)
    end
    persist_job_manifest!(config, data, circuit_info, jobs, paths)
    manifest_text = read(paths["job_manifest"], String)
    occursin("\"error_type\":\"ErrorException\"", manifest_text) ||
        smoke_error("IBM pilot manifest must persist after-submission failure type")
    occursin("Runtime job result retrieval or scoring failed after submission.", manifest_text) ||
        smoke_error("IBM pilot manifest must persist after-submission failure state")
end

println("DS-MFG smoke test passed.")
