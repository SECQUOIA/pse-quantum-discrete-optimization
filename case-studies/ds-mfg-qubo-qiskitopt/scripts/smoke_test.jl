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

function count_csv_data_rows(relpath)
    path = require_file(relpath)
    lines = filter(line -> !isempty(strip(line)), readlines(path))
    length(lines) >= 1 || smoke_error("expected a header row in $relpath")
    return length(lines) - 1
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

function require_text(relpath, required)
    text = read(require_file(relpath), String)
    for value in required
        occursin(value, text) || smoke_error("expected $(relpath) to contain $(value)")
    end
    return text
end

println("Checking Julia and Python package imports...")
Base.pkgversion(QiskitOpt) == v"0.7.0" ||
    smoke_error("expected QiskitOpt v0.7.0, got $(Base.pkgversion(QiskitOpt))")
Base.pkgversion(QUBODrivers) >= v"0.6.1" ||
    smoke_error("expected QUBODrivers v0.6.1 or newer, got $(Base.pkgversion(QUBODrivers))")
isdefined(QUBODrivers, :RandomSeed) ||
    smoke_error("QUBODrivers.RandomSeed() is required for QiskitOpt sampler seeding")
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
    "scripts/artifact_paths.jl",
    "scripts/run_noisy_qaoa_fake_backend.jl",
    "scripts/run_direct_full_qubo_audit.jl",
    "scripts/run_direct_full_qubo_hardware_pilot.jl",
    "scripts/run_direct_full_qubo_qaoa_highread.jl",
    "scripts/update_simulator_hardware_comparison.jl",
    "scripts/enumerate_ip_provenance.py",
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
    "ds_mfg_vqe_reduced_flow_objective_seed74018_metadata",
    "ds_mfg_vqe_reduced_flow_objective_selected_metadata",
    "ds_mfg_reduced_flow_objective",
    "ds_mfg_ibm_qaoa_pilot_fez_4096x3",
    "ds_mfg_simulator_hardware_comparison",
    "ds_mfg_direct_full_qubo_audit",
    "ds_mfg_direct_full_qubo_qaoa_highread",
    "ds_mfg_direct_full_qubo_hardware_pilot",
    "ds_mfg_gurobi_provenance",
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

println("Checking Gurobi provenance artifact...")
gurobi_metadata = require_text(
    "ds_mfg_gurobi_provenance/gurobi_pool_metadata.json",
    (
        "\"artifact_role\": \"ds_mfg_gurobi_pool_provenance\"",
        "\"pool_search_mode\": 2",
        "\"pool_solutions\": 36",
        "\"solution_count\": 36",
        "\"objective_value\": 11.7095",
        "\"mip_gap\": null",
        "\"gurobi_version\": null",
        "\"pool_completeness_certificate\"",
    ),
)
!occursin("/home/", gurobi_metadata) ||
    smoke_error("Gurobi provenance artifact must not include absolute home-directory paths")
ip_enumeration = require_text(
    "ds_mfg_gurobi_provenance/ip_exact_enumeration_summary.json",
    (
        "\"artifact_role\": \"ds_mfg_exact_ip_enumeration_provenance\"",
        "\"total_binary_assignments\": 524288",
        "\"feasible_assignment_count\": 36",
        "\"optimum_objective\": 11.7095",
        "\"optimum_flow_bits\": \"1001110100111100011\"",
        "\"feasible_set_matches_retained_pool\": true",
    ),
)
!occursin("/home/", ip_enumeration) ||
    smoke_error("exact IP enumeration artifact must not include absolute home-directory paths")
count_csv_data_rows("ds_mfg_gurobi_provenance/ip_exact_feasible_flows.csv") == 36 ||
    smoke_error("exact IP enumeration must list 36 feasible flows")
ip_flows = parse_csv_rows("ds_mfg_gurobi_provenance/ip_exact_feasible_flows.csv")
first_ip_flow = first(ip_flows)
require_value(first_ip_flow, "rank", "1")
require_value(first_ip_flow, "ip_obj_value", "11.7095")
require_value(first_ip_flow, "flow_bits", GLOBAL_FLOW_BITS)
last_ip_flow = last(ip_flows)
require_value(last_ip_flow, "rank", "36")
require_value(last_ip_flow, "ip_obj_value", "22.9925")

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

println("Checking direct original-QUBO audit artifacts...")
for filename in (
    "README.md",
    "direct_full_qubo_summary.csv",
    "direct_full_qubo_qaoa_hardware_parameters.json",
    "direct_full_qubo_resource_summary.csv",
    "direct_full_qubo_resource_metadata.json",
)
    require_file(joinpath("ds_mfg_direct_full_qubo_audit", filename))
end

parameter_text = read(
    require_file("ds_mfg_direct_full_qubo_audit/direct_full_qubo_qaoa_hardware_parameters.json"),
    String,
)
occursin("\"artifact_role\":\"direct_full_qubo_qaoa_hardware_parameter_handoff\"", parameter_text) ||
    smoke_error("direct full-QUBO parameter artifact must declare its hardware handoff role")
occursin("\"p\":2", parameter_text) ||
    smoke_error("direct full-QUBO parameter artifact must record p=2")
occursin("\"qiskit_angles_beta_then_gamma\":[0.55500000000000005,0.29299999999999998,-0.48799999999999999,-0.89800000000000002]", parameter_text) ||
    smoke_error("direct full-QUBO parameter artifact must record the fixed QAOA angles")
occursin("fixed-parameter sampling only", parameter_text) ||
    smoke_error("direct full-QUBO parameter artifact must record hardware sampling policy")

direct_rows = parse_csv_rows("ds_mfg_direct_full_qubo_audit/direct_full_qubo_summary.csv")
length(direct_rows) == 20 || smoke_error("expected 20 direct full-QUBO audit rows")
direct_qaoa = only(filter(
    row -> row["source_artifact"] == "ds_mfg_final_sampling_sweep_v2/qaoa_p2_optimizer_reads128_final_reads512_iter25_distribution.csv",
    direct_rows,
))
require_value(direct_qaoa, "problem", "direct_original_36_variable_qubo")
require_value(direct_qaoa, "mode", "cached_local_aer_emulation")
require_int(direct_qaoa, "total_reads", 512)
require_int(direct_qaoa, "projected_pool_hits", 3)
require_int(direct_qaoa, "encoded_pool_hits", 0)
require_int(direct_qaoa, "top50_hits", 5)
require_int(direct_qaoa, "top10_hits", 1)
require_int(direct_qaoa, "global_hits", 0)
require_int(direct_qaoa, "gurobi_pool_feasible_hits", 3)
require_value(direct_qaoa, "best_repaired_match", "gurobi_pool")
require_value(direct_qaoa, "best_repaired_flow_bits", "1001110110011100011")
require_float(direct_qaoa, "best_repaired_qubo_energy", 14.5505)
require_hit_rate_stats(direct_qaoa)

direct_vqe = only(filter(
    row -> row["source_artifact"] == "ds_mfg_vqe_final_sampling_sweep_v2/vqe_full_random_seed73005_optimizer_reads128_final_reads8192_iter25_distribution.csv",
    direct_rows,
))
require_int(direct_vqe, "total_reads", 8192)
require_int(direct_vqe, "projected_pool_hits", 1)
require_int(direct_vqe, "encoded_pool_hits", 0)
require_int(direct_vqe, "top50_hits", 2)
require_int(direct_vqe, "top10_hits", 1)
require_int(direct_vqe, "global_hits", 0)
require_value(direct_vqe, "best_repaired_match", "gurobi_pool")
require_value(direct_vqe, "best_repaired_flow_bits", "1001110110011010011")
require_float(direct_vqe, "best_repaired_qubo_energy", 14.6515)
require_hit_rate_stats(direct_vqe)

println("Checking direct full-QUBO QAOA high-read artifacts...")
for filename in (
    "direct_full_qubo_qaoa_highread_summary.csv",
    "direct_full_qubo_qaoa_optimized_parameters.json",
    "qaoa_direct_p2_seed96001_reads32768_distribution.csv",
    "qaoa_direct_p2_seed96001_reads32768_optimized_parameters.json",
)
    require_file(joinpath("ds_mfg_direct_full_qubo_qaoa_highread", filename))
end
count_csv_data_rows(
    "ds_mfg_direct_full_qubo_qaoa_highread/qaoa_direct_p2_seed96001_reads32768_distribution.csv",
) == 32696 || smoke_error("unexpected direct full-QUBO high-read distribution row count")

highread = only(parse_csv_rows("ds_mfg_direct_full_qubo_qaoa_highread/direct_full_qubo_qaoa_highread_summary.csv"))
require_value(highread, "algorithm", "QAOA_direct_full_qubo_fixed_parameter_highread")
require_value(highread, "parameter_source", "direct_p2")
require_value(
    highread,
    "optimized_parameter_artifact",
    "ds_mfg_direct_full_qubo_qaoa_highread/qaoa_direct_p2_seed96001_reads32768_optimized_parameters.json",
)
require_int(highread, "p", 2)
require_int(highread, "seed", 96001)
require_int(highread, "optimizer_reads", 32)
require_int(highread, "final_reads", 32768)
require_int(highread, "maximum_iterations", 25)
require_int(highread, "total_reads", 32768)
require_int(highread, "unique_states", 32696)
require_int(highread, "top50_hits", 136)
require_int(highread, "top10_hits", 27)
require_int(highread, "global_hits", 4)
require_int(highread, "gurobi_pool_feasible_hits", 99)
require_value(highread, "best_top50_rank", "1")
require_float(highread, "best_repaired_qubo_energy", GLOBAL_OBJECTIVE)
require_value(highread, "best_repaired_match", "global_optimum")
require_value(highread, "best_repaired_flow_bits", GLOBAL_FLOW_BITS)
require_value(highread, "parameter_order", "beta_then_gamma")
require_value(
    highread,
    "optimized_qiskit_angles_beta_then_gamma",
    "0.39325958665050614;0.271838624642501;-0.008647503446280295;0.12401283690343518",
)
require_hit_rate_stats(highread)

optimized_parameter_text = read(
    require_file("ds_mfg_direct_full_qubo_qaoa_highread/direct_full_qubo_qaoa_optimized_parameters.json"),
    String,
)
occursin("\"artifact_role\":\"direct_full_qubo_qaoa_hardware_parameter_handoff\"", optimized_parameter_text) ||
    smoke_error("direct high-read optimized parameter artifact must declare its hardware handoff role")
occursin("\"global_hits\":4", optimized_parameter_text) ||
    smoke_error("direct high-read optimized parameter artifact must record the global-hit validation sample")
occursin("\"qiskit_angles_beta_then_gamma\":[0.39325958665050614,0.27183862464250103,-0.0086475034462802948,0.12401283690343518]", optimized_parameter_text) ||
    smoke_error("direct high-read optimized parameter artifact must record optimized QAOA angles")
occursin("fixed-parameter sampling only", optimized_parameter_text) ||
    smoke_error("direct high-read optimized parameter artifact must record hardware sampling policy")

direct_highread = only(filter(
    row -> row["source_artifact"] == "ds_mfg_direct_full_qubo_qaoa_highread/qaoa_direct_p2_seed96001_reads32768_distribution.csv",
    direct_rows,
))
require_value(direct_highread, "algorithm", "QAOA_direct_full_qubo_fixed_parameter_highread")
require_value(direct_highread, "problem", "direct_original_36_variable_qubo")
require_value(direct_highread, "mode", "cached_local_aer_emulation")
require_int(direct_highread, "total_reads", 32768)
require_int(direct_highread, "projected_pool_hits", 99)
require_int(direct_highread, "encoded_pool_hits", 0)
require_int(direct_highread, "encoded_global_hits", 0)
require_int(direct_highread, "top50_hits", 136)
require_int(direct_highread, "top10_hits", 27)
require_int(direct_highread, "global_hits", 4)
require_int(direct_highread, "gurobi_pool_feasible_hits", 99)
require_value(direct_highread, "best_repaired_rank", "1")
require_float(direct_highread, "best_repaired_qubo_energy", GLOBAL_OBJECTIVE)
require_value(direct_highread, "best_repaired_match", "global_optimum")
require_value(direct_highread, "best_repaired_flow_bits", GLOBAL_FLOW_BITS)
require_hit_rate_stats(direct_highread)

println("Checking direct full-QUBO IBM hardware pilot artifacts...")
for filename in ("README.md", "LEGACY_ARTIFACT.md", "job_manifest.json", "backend_metadata.json", "raw_counts.csv", "scored_counts.csv", "summary.csv")
    require_file(joinpath("ds_mfg_direct_full_qubo_hardware_pilot", filename))
end

direct_hardware_legacy_notice = read(
    require_file("ds_mfg_direct_full_qubo_hardware_pilot/LEGACY_ARTIFACT.md"),
    String,
)
occursin("quarantined legacy hardware artifact", direct_hardware_legacy_notice) ||
    smoke_error("direct full-QUBO hardware directory must include a legacy quarantine notice")
occursin("Do not treat these cached counts as reproducible from the current script", direct_hardware_legacy_notice) ||
    smoke_error("direct full-QUBO hardware legacy notice must mark cached counts non-reproducible from current script")

direct_hardware_summary = only(parse_csv_rows("ds_mfg_direct_full_qubo_hardware_pilot/summary.csv"))
require_value(direct_hardware_summary, "algorithm", "QAOA_direct_full_qubo_IBM_handoff")
require_value(direct_hardware_summary, "mode", "hardware")
require_value(direct_hardware_summary, "backend", "ibm_fez")
require_value(
    direct_hardware_summary,
    "parameter_artifact",
    "ds_mfg_direct_full_qubo_qaoa_highread/direct_full_qubo_qaoa_optimized_parameters.json",
)
require_int(direct_hardware_summary, "p", 2)
require_int(direct_hardware_summary, "final_reads", 4096)
require_int(direct_hardware_summary, "repeats", 1)
require_value(direct_hardware_summary, "transpile_seeds", "92001")
require_int(direct_hardware_summary, "submitted_jobs", 1)
require_int(direct_hardware_summary, "unique_full_states", 4096)
require_int(direct_hardware_summary, "unique_flows", 4057)
require_int(direct_hardware_summary, "total_reads", 4096)
require_int(direct_hardware_summary, "projected_pool_hits", 0)
require_int(direct_hardware_summary, "encoded_pool_hits", 0)
require_int(direct_hardware_summary, "encoded_global_hits", 0)
require_int(direct_hardware_summary, "top50_hits", 0)
require_int(direct_hardware_summary, "top10_hits", 0)
require_int(direct_hardware_summary, "global_hits", 0)
require_int(direct_hardware_summary, "gurobi_pool_feasible_hits", 0)
require_value(direct_hardware_summary, "best_repaired_rank", "")
require_float(direct_hardware_summary, "best_repaired_qubo_energy", 91.7295)
require_value(direct_hardware_summary, "best_repaired_match", "not_top50")
require_value(direct_hardware_summary, "best_repaired_flow_bits", "1010110110000000000")
require_hit_rate_stats(direct_hardware_summary; time_key = "elapsed_sec")
count_csv_data_rows("ds_mfg_direct_full_qubo_hardware_pilot/raw_counts.csv") == 4096 ||
    smoke_error("unexpected direct full-QUBO hardware raw-count row count")
count_csv_data_rows("ds_mfg_direct_full_qubo_hardware_pilot/scored_counts.csv") == 4096 ||
    smoke_error("unexpected direct full-QUBO hardware scored-count row count")

direct_hardware_manifest = read(require_file("ds_mfg_direct_full_qubo_hardware_pilot/job_manifest.json"), String)
occursin("\"mode\":\"hardware\"", direct_hardware_manifest) ||
    smoke_error("direct full-QUBO hardware manifest must record hardware mode")
occursin("\"backend_name\":\"ibm_fez\"", direct_hardware_manifest) ||
    smoke_error("direct full-QUBO hardware manifest must record ibm_fez backend")
occursin("\"output_dir\":\"ds_mfg_direct_full_qubo_hardware_pilot\"", direct_hardware_manifest) ||
    smoke_error("direct full-QUBO hardware manifest must store a repo-relative output directory")
!occursin("/home/", direct_hardware_manifest) ||
    smoke_error("direct full-QUBO hardware manifest must not include local home paths")
occursin("\"run_hardware\":true", direct_hardware_manifest) ||
    smoke_error("direct full-QUBO hardware manifest must record hardware execution")
occursin("\"submitted\":true", direct_hardware_manifest) ||
    smoke_error("direct full-QUBO hardware manifest must record submitted job")
occursin("\"status\":\"DONE\"", direct_hardware_manifest) ||
    smoke_error("direct full-QUBO hardware manifest must record completed job")
occursin("\"qiskit_angles_beta_then_gamma\":[0.39325958665050614,0.27183862464250103,-0.0086475034462802948,0.12401283690343518]", direct_hardware_manifest) ||
    smoke_error("direct full-QUBO hardware manifest must record optimized QAOA angles")

direct_hardware_backend = read(require_file("ds_mfg_direct_full_qubo_hardware_pilot/backend_metadata.json"), String)
occursin("\"backend_name_resolved\":\"ibm_fez\"", direct_hardware_backend) ||
    smoke_error("direct full-QUBO backend metadata must record resolved backend")
occursin("\"num_qubits\":156", direct_hardware_backend) ||
    smoke_error("direct full-QUBO backend metadata must record ibm_fez qubit count")
occursin("\"simulator\":false", direct_hardware_backend) ||
    smoke_error("direct full-QUBO backend metadata must record non-simulator backend")
for text in (direct_hardware_manifest, direct_hardware_backend)
    !occursin("QISKIT_IBM_TOKEN", text) ||
        smoke_error("direct full-QUBO hardware artifacts must not include token environment names")
    !occursin("qiskit-ibm.json", text) ||
        smoke_error("direct full-QUBO hardware artifacts must not include account file paths")
    !occursin("crn:v1:", lowercase(text)) ||
        smoke_error("direct full-QUBO hardware artifacts must not include Runtime instance CRNs")
end

resource_row = only(parse_csv_rows("ds_mfg_direct_full_qubo_audit/direct_full_qubo_resource_summary.csv"))
require_value(resource_row, "problem", "direct_original_36_variable_qubo")
require_value(resource_row, "mode", "resource_audit")
require_value(resource_row, "fake_backend", "FakeFez")
require_int(resource_row, "logical_qubits", 36)
require_int(resource_row, "logical_rzz_count", 202)
require_int(resource_row, "logical_two_qubit_ops", 202)
require_value(resource_row, "statevector_amplitudes", "68719476736")
require_value(resource_row, "statevector_complex128_bytes", "1099511627776")
require_value(resource_row, "density_matrix_elements", "4722366482869645213696")
require_value(resource_row, "density_matrix_complex128_bytes", "75557863725914323419136")
require_value(resource_row, "noisy_method_considered", "qiskit_aer.AerSimulator.from_backend(FakeFez)")
require_value(resource_row, "transpile_requested", "true")
require_value(resource_row, "transpile_status", "DONE")
require_int(resource_row, "transpiled_depth", 1344)
require_int(resource_row, "transpiled_cz_count", 1157)
require_value(resource_row, "transpile_elapsed_sec", "not_recorded_for_reproducible_artifacts")
require_value(resource_row, "simulation_status", "not_run")
resource_metadata = read(require_file("ds_mfg_direct_full_qubo_audit/direct_full_qubo_resource_metadata.json"), String)
occursin("\"n_qubits\":36", resource_metadata) ||
    smoke_error("direct full-QUBO resource metadata must record 36 qubits")
occursin("\"cz_count\":1157", resource_metadata) ||
    smoke_error("direct full-QUBO resource metadata must record FakeFez CZ count")
occursin("\"dense_density_matrix_elements\":\"4722366482869645213696\"", resource_metadata) ||
    smoke_error("direct full-QUBO metadata must record density-matrix size")
occursin("not_recorded_for_reproducible_artifacts", resource_metadata) ||
    smoke_error("direct full-QUBO metadata must avoid volatile transpile runtime")
!occursin("created_at_utc", resource_metadata) ||
    smoke_error("direct full-QUBO metadata must avoid volatile timestamps")
occursin("Aer automatic or MPS simulation is entanglement-dependent", resource_metadata) ||
    smoke_error("direct full-QUBO metadata must state noisy-simulator feasibility")

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

println("Checking selected VQE optimized-parameter metadata artifact...")
for filename in (
    "vqe_reduced_seed74018_metadata.json",
    "vqe_reduced_top50_sampling_summary.csv",
    "vqe_reduced_top50_sampling_hits.csv",
)
    require_file(joinpath("ds_mfg_vqe_reduced_flow_objective_seed74018_metadata", filename))
end
vqe_metadata_summary = only(parse_csv_rows(
    "ds_mfg_vqe_reduced_flow_objective_seed74018_metadata/vqe_reduced_top50_sampling_summary.csv",
))
require_value(vqe_metadata_summary, "algorithm", "VQE_reduced_surrogate_top50_revisit")
require_int(vqe_metadata_summary, "seed", 74018)
require_int(vqe_metadata_summary, "optimizer_reads", 128)
require_int(vqe_metadata_summary, "final_reads", 32768)
require_int(vqe_metadata_summary, "maximum_iterations", 25)
require_int(vqe_metadata_summary, "total_reads", 32768)
require_int(vqe_metadata_summary, "top50_hits", 5)
require_int(vqe_metadata_summary, "global_hits", 0)
require_value(vqe_metadata_summary, "best_top50_rank", "12")
require_value(vqe_metadata_summary, "best_top50_match", "gurobi_pool")
require_hit_rate_stats(vqe_metadata_summary)

vqe_metadata_text = require_text(
    "ds_mfg_vqe_reduced_flow_objective_seed74018_metadata/vqe_reduced_seed74018_metadata.json",
    (
        "\"artifact_role\":\"ds_mfg_reduced_surrogate_vqe_seed_metadata\"",
        "\"seed\":74018",
        "\"optimized_parameter_artifact_status\":\"retained_from_qiskitopt_sampleset_metadata\"",
        "\"initial_parameters\"",
        "\"optimized_parameters\"",
        "\"parameter_names\"",
        "\"values\"",
        "\"backend_configuration\"",
        "\"matrix_product_state\"",
    ),
)
!occursin("/home/", vqe_metadata_text) ||
    smoke_error("selected VQE metadata artifact must not include absolute home-directory paths")
!occursin("QISKIT_IBM_TOKEN", vqe_metadata_text) ||
    smoke_error("selected VQE metadata artifact must not include token environment names")

println("Checking selected final-budget VQE optimized-parameter metadata artifacts...")
for filename in (
    "vqe_reduced_seed74018_metadata.json",
    "vqe_reduced_seed74007_metadata.json",
    "vqe_reduced_seed74001_metadata.json",
    "vqe_reduced_top50_sampling_summary.csv",
    "vqe_reduced_top50_sampling_hits.csv",
)
    require_file(joinpath("ds_mfg_vqe_reduced_flow_objective_selected_metadata", filename))
end
selected_vqe_rows = parse_csv_rows(
    "ds_mfg_vqe_reduced_flow_objective_selected_metadata/vqe_reduced_top50_sampling_summary.csv",
)
length(selected_vqe_rows) == 3 || smoke_error("expected three selected VQE metadata rows")
all(row["algorithm"] == "VQE_reduced_surrogate_top50_revisit" for row in selected_vqe_rows) ||
    smoke_error("selected VQE metadata rows must use the revisit algorithm label")
all(row["optimizer_reads"] == "128" for row in selected_vqe_rows) ||
    smoke_error("selected VQE metadata rows must record 128 optimizer reads")
all(row["final_reads"] == "524288" for row in selected_vqe_rows) ||
    smoke_error("selected VQE metadata rows must record 524288 final reads")
all(row["maximum_iterations"] == "25" for row in selected_vqe_rows) ||
    smoke_error("selected VQE metadata rows must record 25 maximum iterations")
sum(parse(Int, row["global_hits"]) for row in selected_vqe_rows) == 4 ||
    smoke_error("unexpected selected VQE metadata global-hit total")
sum(parse(Int, row["gurobi_pool_feasible_hits"]) for row in selected_vqe_rows) == 266 ||
    smoke_error("unexpected selected VQE metadata Gurobi-pool feasible-hit total")
selected_by_seed = Dict(row["seed"] => row for row in selected_vqe_rows)
require_value(selected_by_seed["74018"], "best_top50_rank", "3")
require_int(selected_by_seed["74018"], "global_hits", 0)
require_value(selected_by_seed["74007"], "best_top50_match", "global_optimum")
require_int(selected_by_seed["74007"], "global_hits", 1)
require_value(selected_by_seed["74001"], "best_top50_match", "global_optimum")
require_int(selected_by_seed["74001"], "global_hits", 3)
for row in selected_vqe_rows
    require_hit_rate_stats(row)
end
for seed in ("74018", "74007", "74001")
    metadata_text = require_text(
        joinpath(
            "ds_mfg_vqe_reduced_flow_objective_selected_metadata",
            "vqe_reduced_seed$(seed)_metadata.json",
        ),
        (
            "\"artifact_role\":\"ds_mfg_reduced_surrogate_vqe_seed_metadata\"",
            "\"seed\":$(seed)",
            "\"optimized_parameter_artifact_status\":\"retained_from_qiskitopt_sampleset_metadata\"",
            "\"optimized_parameters\"",
            "\"parameter_names\"",
            "\"values\"",
            "\"matrix_product_state\"",
        ),
    )
    !occursin("/home/", metadata_text) ||
        smoke_error("selected VQE metadata artifact for seed $(seed) must not include absolute home-directory paths")
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

println("Checking IBM QAOA hardware pilot artifacts...")
for filename in ("job_manifest.json", "backend_metadata.json", "raw_counts.csv", "scored_counts.csv", "summary.csv")
    require_file(joinpath("ds_mfg_ibm_qaoa_pilot_fez_4096x3", filename))
end

hardware_summary = only(parse_csv_rows("ds_mfg_ibm_qaoa_pilot_fez_4096x3/summary.csv"))
require_value(hardware_summary, "algorithm", "QAOA_reduced_surrogate_JuliQAOA_IBM_pilot")
require_value(hardware_summary, "mode", "hardware")
require_value(hardware_summary, "backend", "ibm_fez")
require_value(hardware_summary, "angle_target", "top10")
require_int(hardware_summary, "p", 5)
require_int(hardware_summary, "final_reads", 4096)
require_int(hardware_summary, "repeats", 3)
require_value(hardware_summary, "transpile_seeds", "92001;92002;92003")
require_int(hardware_summary, "submitted_jobs", 9)
require_int(hardware_summary, "unique_flows", 34573)
require_int(hardware_summary, "total_reads", 36864)
require_int(hardware_summary, "top50_hits", 6)
require_int(hardware_summary, "top10_hits", 1)
require_int(hardware_summary, "global_hits", 0)
require_int(hardware_summary, "gurobi_pool_feasible_hits", 4)
require_value(hardware_summary, "best_top50_rank", "2")
require_float(hardware_summary, "best_top50_exact_repaired_qubo_energy", 11.8105)
require_value(hardware_summary, "best_top50_match", "gurobi_pool")
require_value(hardware_summary, "best_top50_flow_bits", "1001110100111010011")
for (column, expected) in (
    "top50_hit_rate" => "0.000162760416667",
    "top50_hit_rate_wilson95_low" => "7.45966258791e-05",
    "top50_hit_rate_wilson95_high" => "0.000355085671588",
    "top50_tts99_sec" => "56.4001563512",
    "top10_hit_rate" => "2.71267361111e-05",
    "top10_hit_rate_wilson95_low" => "4.78855106823e-06",
    "top10_hit_rate_wilson95_high" => "0.000153654650094",
    "top10_tts99_sec" => "338.422866628",
    "global_hit_rate" => "0",
    "global_hit_rate_wilson95_low" => "0",
    "global_hit_rate_wilson95_high" => "0.000104195381902",
    "global_tts99_sec" => "Inf",
    "gurobi_pool_feasible_hit_rate" => "0.000108506944444",
    "gurobi_pool_feasible_hit_rate_wilson95_low" => "4.21970151315e-05",
    "gurobi_pool_feasible_hit_rate_wilson95_high" => "0.000278989643814",
    "gurobi_pool_feasible_tts99_sec" => "84.6022280288",
)
    require_value(hardware_summary, column, expected)
end

count_csv_data_rows("ds_mfg_ibm_qaoa_pilot_fez_4096x3/raw_counts.csv") == 36539 ||
    smoke_error("unexpected hardware raw-count row count")
count_csv_data_rows("ds_mfg_ibm_qaoa_pilot_fez_4096x3/scored_counts.csv") == 36539 ||
    smoke_error("unexpected hardware scored-count row count")

hardware_manifest = read(require_file("ds_mfg_ibm_qaoa_pilot_fez_4096x3/job_manifest.json"), String)
occursin("\"mode\":\"hardware\"", hardware_manifest) ||
    smoke_error("hardware manifest must record hardware mode")
occursin("\"backend_name\":\"ibm_fez\"", hardware_manifest) ||
    smoke_error("hardware manifest must record ibm_fez backend")
occursin("\"submitted\":true", hardware_manifest) ||
    smoke_error("hardware manifest must record submitted jobs")
occursin("\"status\":\"DONE\"", hardware_manifest) ||
    smoke_error("hardware manifest must record completed jobs")
occursin("\"output_dir\":\"ds_mfg_ibm_qaoa_pilot_fez_4096x3\"", hardware_manifest) ||
    smoke_error("hardware manifest must record a repository-relative output directory")
occursin("\"schema_version\":1", hardware_manifest) ||
    smoke_error("cached hardware manifest must record schema version 1 before rerun")

hardware_backend = read(require_file("ds_mfg_ibm_qaoa_pilot_fez_4096x3/backend_metadata.json"), String)
occursin("\"backend_name_resolved\":\"ibm_fez\"", hardware_backend) ||
    smoke_error("hardware backend metadata must record resolved backend")
occursin("\"num_qubits\":156", hardware_backend) ||
    smoke_error("hardware backend metadata must record ibm_fez qubit count")
occursin("\"simulator\":false", hardware_backend) ||
    smoke_error("hardware backend metadata must record non-simulator backend")

for text in (hardware_manifest, hardware_backend)
    !occursin("QISKIT_IBM_TOKEN", text) ||
        smoke_error("hardware artifacts must not include token environment names")
    !occursin("qiskit-ibm.json", text) ||
        smoke_error("hardware artifacts must not include account file paths")
    !occursin("crn:v1:", lowercase(text)) ||
        smoke_error("hardware artifacts must not include Runtime instance CRNs")
    !occursin("/home/", text) ||
        smoke_error("hardware artifacts must not include absolute home-directory paths")
end

println("Checking simulator-to-hardware comparison artifacts...")
noisy_artifact_dir = "ds_mfg_fake_fez_qaoa_noisy_4096x3"
for filename in ("job_manifest.json", "backend_metadata.json", "raw_counts.csv", "scored_counts.csv", "summary.csv")
    require_file(joinpath(noisy_artifact_dir, filename))
end

noisy_summary = only(parse_csv_rows(joinpath(noisy_artifact_dir, "summary.csv")))
require_value(noisy_summary, "algorithm", "QAOA_reduced_surrogate_JuliQAOA_FakeFez_noisy_simulation")
require_value(noisy_summary, "mode", "model_based_noisy_simulation")
require_value(noisy_summary, "fake_backend", "FakeFez")
require_value(noisy_summary, "angle_target", "top10")
require_int(noisy_summary, "p", 5)
require_int(noisy_summary, "final_reads", 4096)
require_int(noisy_summary, "repeats", 3)
require_value(noisy_summary, "transpile_seeds", "92001;92002;92003")
require_int(noisy_summary, "completed_jobs", 9)
require_int(noisy_summary, "unique_flows", 11172)
require_int(noisy_summary, "total_reads", 36864)
require_int(noisy_summary, "top50_hits", 516)
require_int(noisy_summary, "top10_hits", 108)
require_int(noisy_summary, "global_hits", 12)
require_int(noisy_summary, "gurobi_pool_feasible_hits", 447)
require_value(noisy_summary, "global_hit_rate", "0.000325520833333")
require_value(noisy_summary, "best_top50_match", "global_optimum")
require_value(noisy_summary, "best_top50_flow_bits", GLOBAL_FLOW_BITS)
require_int(noisy_summary, "transpiled_depth_min", 1442)
require_int(noisy_summary, "transpiled_depth_max", 1557)
require_int(noisy_summary, "transpiled_cz_min", 1136)
require_int(noisy_summary, "transpiled_cz_max", 1156)
count_csv_data_rows(joinpath(noisy_artifact_dir, "raw_counts.csv")) == 35636 ||
    smoke_error("unexpected FakeFez noisy raw-count row count")
count_csv_data_rows(joinpath(noisy_artifact_dir, "scored_counts.csv")) == 35636 ||
    smoke_error("unexpected FakeFez noisy scored-count row count")
noisy_manifest = read(require_file(joinpath(noisy_artifact_dir, "job_manifest.json")), String)
occursin("\"mode\":\"model_based_noisy_simulation\"", noisy_manifest) ||
    smoke_error("FakeFez noisy manifest must record model-based simulation mode")
occursin("\"run_simulation\":true", noisy_manifest) ||
    smoke_error("FakeFez noisy manifest must record simulation enabled")
occursin("\"status\":\"DONE\"", noisy_manifest) ||
    smoke_error("FakeFez noisy manifest must record completed jobs")
!occursin("\"status\":\"FAILED\"", noisy_manifest) ||
    smoke_error("FakeFez noisy manifest must not record failed jobs")
occursin("\"output_dir\":\"ds_mfg_fake_fez_qaoa_noisy_4096x3\"", noisy_manifest) ||
    smoke_error("FakeFez noisy manifest must record a repository-relative output directory")
!occursin("/home/", noisy_manifest) ||
    smoke_error("FakeFez noisy manifest must not include absolute home-directory paths")

comparison_rows = parse_csv_rows("ds_mfg_simulator_hardware_comparison/simulator_hardware_comparison_summary.csv")
length(comparison_rows) == 3 || smoke_error("expected three simulator/hardware comparison rows")
ideal_comparison = only(filter(row -> row["evidence_tier"] == "ideal_simulator", comparison_rows))
require_value(ideal_comparison, "row_status", "cached_result")
require_value(ideal_comparison, "source_artifact", "ds_mfg_qaoa_juliqaoa_objective_transfer_highread/qaoa_juliqaoa_transfer_summary.csv")
require_int(ideal_comparison, "total_reads", 262144)
require_int(ideal_comparison, "top50_hits", 62597)
require_int(ideal_comparison, "top10_hits", 14326)
require_int(ideal_comparison, "global_hits", 1007)
require_value(ideal_comparison, "best_top50_flow_bits", GLOBAL_FLOW_BITS)

noisy_comparison = only(filter(row -> row["evidence_tier"] == "model_based_noisy_simulator", comparison_rows))
require_value(noisy_comparison, "row_status", "cached_result")
require_value(noisy_comparison, "source_artifact", "ds_mfg_fake_fez_qaoa_noisy_4096x3/summary.csv")
require_int(noisy_comparison, "total_reads", 36864)
require_int(noisy_comparison, "top50_hits", 516)
require_int(noisy_comparison, "top10_hits", 108)
require_int(noisy_comparison, "global_hits", 12)
require_value(noisy_comparison, "best_top50_flow_bits", GLOBAL_FLOW_BITS)
occursin("Model-based fake-backend simulation", noisy_comparison["notes"]) ||
    smoke_error("noisy comparison row must label cached output as model-based simulation")

hardware_comparison = only(filter(row -> row["evidence_tier"] == "hardware", comparison_rows))
require_value(hardware_comparison, "row_status", "cached_result")
require_value(hardware_comparison, "source_artifact", "ds_mfg_ibm_qaoa_pilot_fez_4096x3/summary.csv")
require_int(hardware_comparison, "total_reads", 36864)
require_int(hardware_comparison, "top50_hits", 6)
require_int(hardware_comparison, "top10_hits", 1)
require_int(hardware_comparison, "global_hits", 0)
require_value(hardware_comparison, "best_top50_flow_bits", "1001110100111010011")

degradation_rows = parse_csv_rows("ds_mfg_simulator_hardware_comparison/simulator_hardware_degradation_summary.csv")
length(degradation_rows) == 4 || smoke_error("expected four simulator/hardware degradation rows")
top50_degradation = only(filter(row -> row["event"] == "top50", degradation_rows))
require_value(top50_degradation, "ideal_aer_hit_rate", "0.238788604736")
require_value(top50_degradation, "noisy_model_hit_rate", "0.0139973958333")
require_value(top50_degradation, "noisy_to_ideal_hit_rate_ratio", "0.0586183576423")
require_value(top50_degradation, "hardware_hit_rate", "0.000162760416667")
require_value(top50_degradation, "hardware_to_ideal_hit_rate_ratio", "0.000681608809794")
global_degradation = only(filter(row -> row["event"] == "global", degradation_rows))
require_value(global_degradation, "noisy_model_hit_rate", "0.000325520833333")
require_value(global_degradation, "ideal_expected_hits_at_hardware_reads", "141.609375")
require_value(global_degradation, "hardware_hits_minus_ideal_expected", "-141.609375")

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
    occursin("\"output_dir\":\"$(basename(pilot_output_dir))\"", manifest_text) ||
        smoke_error("IBM pilot manifest must record only a public output directory")
    !occursin(pilot_output_dir, manifest_text) ||
        smoke_error("IBM pilot manifest must not include absolute output paths")
    occursin("\"schema_version\":2", manifest_text) ||
        smoke_error("IBM pilot dry-run manifest must record schema version 2")
    occursin("\"transpiler_optimization_level\":3", manifest_text) ||
        smoke_error("IBM pilot dry-run manifest must record transpiler optimization level")
    occursin("\"readout_mitigation\":false", manifest_text) ||
        smoke_error("IBM pilot dry-run manifest must record readout mitigation policy")
    occursin("\"error_mitigation\":false", manifest_text) ||
        smoke_error("IBM pilot dry-run manifest must record error mitigation policy")
    occursin("\"transpiled_circuit\":null", manifest_text) ||
        smoke_error("IBM pilot dry-run manifest must include transpiled-circuit schema field")
    occursin("\"queue_timing\"", manifest_text) ||
        smoke_error("IBM pilot dry-run manifest must include queue timing schema field")

    summary_lines = filter(line -> !isempty(strip(line)), readlines(joinpath(pilot_output_dir, "summary.csv")))
    length(summary_lines) == 2 || smoke_error("IBM pilot dry-run summary must contain one data row")
    occursin("QAOA_reduced_surrogate_JuliQAOA_IBM_pilot,dry_run,ibm_brisbane,top10,5,64,1,123,0,", summary_lines[2]) ||
        smoke_error("IBM pilot dry-run summary row has unexpected configuration values")
end

println("Checking FakeFez noisy QAOA dry-run schema...")
mktempdir() do noisy_output_dir
    withenv(
        "DSMFG_NOISY_OUTPUT_DIR" => noisy_output_dir,
        "DSMFG_RUN_NOISY_SIMULATION" => "false",
    ) do
        run(`$(Base.julia_cmd()) --project=$(ROOT) scripts/run_noisy_qaoa_fake_backend.jl`)
    end

    for filename in ("job_manifest.json", "backend_metadata.json", "raw_counts.csv", "scored_counts.csv", "summary.csv")
        path = joinpath(noisy_output_dir, filename)
        isfile(path) || smoke_error("FakeFez noisy dry run did not write $(filename)")
        filesize(path) > 0 || smoke_error("FakeFez noisy dry-run file is empty: $(filename)")
    end

    manifest_text = read(joinpath(noisy_output_dir, "job_manifest.json"), String)
    occursin("\"mode\":\"dry_run\"", manifest_text) ||
        smoke_error("FakeFez noisy manifest must record dry_run mode")
    occursin("\"fake_backend_class\":\"FakeFez\"", manifest_text) ||
        smoke_error("FakeFez noisy manifest must record FakeFez backend class")
    occursin("\"run_simulation\":false", manifest_text) ||
        smoke_error("FakeFez noisy dry-run manifest must not mark simulation enabled")
    occursin("model-based fake-backend simulation", manifest_text) ||
        smoke_error("FakeFez noisy manifest must include model-based interpretation warning")
    occursin("\"output_dir\":\"$(basename(noisy_output_dir))\"", manifest_text) ||
        smoke_error("FakeFez noisy manifest must record only a public output directory")
    !occursin(noisy_output_dir, manifest_text) ||
        smoke_error("FakeFez noisy manifest must not include absolute output paths")

    summary_lines = filter(line -> !isempty(strip(line)), readlines(joinpath(noisy_output_dir, "summary.csv")))
    length(summary_lines) == 2 || smoke_error("FakeFez noisy dry-run summary must contain one data row")
    occursin("QAOA_reduced_surrogate_JuliQAOA_FakeFez_noisy_simulation,dry_run,FakeFez,top10,5,4096,3,92001;92002;92003,93001,3,0,", summary_lines[2]) ||
        smoke_error("FakeFez noisy dry-run summary row has unexpected configuration values")
end

println("Checking direct full-QUBO audit dry-run schema...")
mktempdir() do direct_output_dir
    withenv(
        "DSMFG_DIRECT_FULL_QUBO_OUTPUT_DIR" => direct_output_dir,
        "DSMFG_DIRECT_FULL_QUBO_TRANSPILE" => "false",
    ) do
        run(`$(Base.julia_cmd()) --project=$(ROOT) scripts/run_direct_full_qubo_audit.jl`)
    end

    for filename in (
        "direct_full_qubo_summary.csv",
        "direct_full_qubo_qaoa_hardware_parameters.json",
        "direct_full_qubo_resource_summary.csv",
        "direct_full_qubo_resource_metadata.json",
    )
        path = joinpath(direct_output_dir, filename)
        isfile(path) || smoke_error("direct full-QUBO dry run did not write $(filename)")
        filesize(path) > 0 || smoke_error("direct full-QUBO dry-run file is empty: $(filename)")
    end

    dry_summary_rows = parse_csv_rows(joinpath(direct_output_dir, "direct_full_qubo_summary.csv"))
    length(dry_summary_rows) == 20 || smoke_error("direct full-QUBO dry run must write 20 summary rows")
    dry_resource_row = only(parse_csv_rows(joinpath(direct_output_dir, "direct_full_qubo_resource_summary.csv")))
    require_value(dry_resource_row, "transpile_requested", "false")
    require_value(dry_resource_row, "transpile_status", "not_requested")
    require_int(dry_resource_row, "logical_qubits", 36)
    require_int(dry_resource_row, "logical_rzz_count", 202)
    require_value(dry_resource_row, "density_matrix_elements", "4722366482869645213696")
    require_value(dry_resource_row, "transpile_elapsed_sec", "")
    require_value(dry_resource_row, "simulation_status", "not_run")
end

println("Checking direct full-QUBO QAOA hardware dry-run schema...")
mktempdir() do direct_hardware_output_dir
    withenv(
        "QISKIT_IBM_BACKEND" => "ibm_brisbane",
        "QISKIT_IBM_INSTANCE" => "",
        "DSMFG_DIRECT_HARDWARE_FINAL_READS" => "64",
        "DSMFG_DIRECT_HARDWARE_REPEATS" => "1",
        "DSMFG_DIRECT_HARDWARE_TRANSPILE_SEEDS" => "123",
        "DSMFG_DIRECT_HARDWARE_OUTPUT_DIR" => direct_hardware_output_dir,
        "DSMFG_RUN_DIRECT_FULL_QUBO_HARDWARE" => "false",
    ) do
        run(`$(Base.julia_cmd()) --project=$(ROOT) scripts/run_direct_full_qubo_hardware_pilot.jl`)
    end

    for filename in ("job_manifest.json", "backend_metadata.json", "raw_counts.csv", "scored_counts.csv", "summary.csv")
        path = joinpath(direct_hardware_output_dir, filename)
        isfile(path) || smoke_error("direct full-QUBO hardware dry run did not write $(filename)")
        filesize(path) > 0 || smoke_error("direct full-QUBO hardware dry-run file is empty: $(filename)")
    end

    manifest_text = read(joinpath(direct_hardware_output_dir, "job_manifest.json"), String)
    occursin("\"mode\":\"dry_run\"", manifest_text) ||
        smoke_error("direct full-QUBO hardware manifest must record dry_run mode")
    occursin("\"run_hardware\":false", manifest_text) ||
        smoke_error("direct full-QUBO hardware dry-run manifest must not mark hardware enabled")
    occursin("\"output_dir\":\"$(basename(direct_hardware_output_dir))\"", manifest_text) ||
        smoke_error("direct full-QUBO hardware dry-run manifest must avoid absolute output paths")
    !occursin(direct_hardware_output_dir, manifest_text) ||
        smoke_error("direct full-QUBO hardware dry-run manifest must not include the temp output path")
    occursin("\"n_qubits\":36", manifest_text) ||
        smoke_error("direct full-QUBO hardware manifest must record the 36-qubit problem")
    occursin("\"p\":2", manifest_text) ||
        smoke_error("direct full-QUBO hardware manifest must record p=2")
    occursin("\"source\":\"QiskitOpt.QAOA.fixed_parameter_circuit\"", manifest_text) ||
        smoke_error("direct full-QUBO hardware manifest must record QiskitOpt fixed-parameter circuit source")
    occursin("\"mode\":\"fixed_parameter_circuit\"", manifest_text) ||
        smoke_error("direct full-QUBO hardware manifest must record QiskitOpt fixed-parameter metadata")
    occursin("\"qiskit_minimization_sign\":1", manifest_text) ||
        smoke_error("direct full-QUBO hardware manifest must record QiskitOpt minimization sign")
    occursin("\"values_aligned_to\":\"parameter_names\"", manifest_text) ||
        smoke_error("direct full-QUBO hardware manifest must record QiskitOpt parameter alignment")
    occursin("\"submitted\":false", manifest_text) ||
        smoke_error("direct full-QUBO hardware dry-run manifest must not mark jobs submitted")
    !occursin("QISKIT_IBM_TOKEN", manifest_text) ||
        smoke_error("direct full-QUBO hardware manifest must not include token environment names")
    !occursin("qiskit-ibm.json", manifest_text) ||
        smoke_error("direct full-QUBO hardware manifest must not include account file paths")

    summary_lines = filter(line -> !isempty(strip(line)), readlines(joinpath(direct_hardware_output_dir, "summary.csv")))
    length(summary_lines) == 2 || smoke_error("direct full-QUBO hardware dry-run summary must contain one data row")
    occursin(
        "QAOA_direct_full_qubo_IBM_handoff,dry_run,ibm_brisbane,ds_mfg_direct_full_qubo_qaoa_highread/direct_full_qubo_qaoa_optimized_parameters.json,2,64,1,123,0,",
        summary_lines[2],
    ) || smoke_error("direct full-QUBO hardware dry-run summary row has unexpected configuration values")
end

include(joinpath(@__DIR__, "run_ibm_qaoa_pilot.jl"))

function py_probability_dict_to_julia(py_dict)
    converted = Dict{String,Float64}()
    for item in py_dict.items()
        converted[pyconvert(String, item[0])] = pyconvert(Float64, item[1])
    end
    return converted
end

function validate_ibm_pilot_circuit_probabilities()
    data = load_problem_data()
    circuit = build_qaoa_circuit(data)
    circuit_without_measurements = circuit.remove_final_measurements(; inplace = false)
    statevector = PythonCall.pyimport("qiskit.quantum_info").Statevector.from_instruction(circuit_without_measurements)
    probabilities = py_probability_dict_to_julia(statevector.probabilities_dict())
    observed = Dict(
        "top50_probability" => 0.0,
        "top10_probability" => 0.0,
        "global_probability" => 0.0,
    )

    for (qiskit_bitstring, probability) in probabilities
        flow_bits = qiskit_key_to_flow_bits(qiskit_bitstring, data.scalars.n)
        hit = get(data.top_flows, flow_bits, nothing)
        isnothing(hit) && continue
        observed["top50_probability"] += probability
        hit.rank <= 10 && (observed["top10_probability"] += probability)
        hit.rank == 1 && (observed["global_probability"] += probability)
    end

    for key in keys(observed)
        expected = parse(Float64, data.angle_record[key])
        isapprox(observed[key], expected; atol = 1.0e-3, rtol = 0.0) ||
            smoke_error("IBM pilot circuit $(key) expected $(expected), got $(observed[key])")
    end
end

println("Checking IBM QAOA pilot circuit probabilities...")
validate_ibm_pilot_circuit_probabilities()

println("Checking IBM QAOA pilot manifest durability helpers...")
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
    occursin("\"output_dir\":\"$(basename(pilot_output_dir))\"", manifest_text) ||
        smoke_error("IBM pilot durability helper must use public output directory paths")
    !occursin(pilot_output_dir, manifest_text) ||
        smoke_error("IBM pilot durability helper must not persist absolute output paths")

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
