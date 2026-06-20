#!/usr/bin/env julia

import Pkg

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(joinpath(@__DIR__, ".."))
    Pkg.instantiate()
end

using Printf
using QiskitOpt
using QiskitOpt: QAOA
using QUBODrivers
using QUBOTools
using PythonCall: pyconvert

include(joinpath(@__DIR__, "run_direct_full_qubo_audit.jl"))

const MOI = QUBODrivers.MOI
const SampleReads = QUBOTools.__moi_num_reads()
const DIRECT_HIGHREAD_ALGORITHM = "QAOA_direct_full_qubo_fixed_parameter_highread"
const DEFAULT_HIGHREAD_OUTPUT_DIR = "ds_mfg_direct_full_qubo_qaoa_highread"

struct DirectHighreadConfig
    output_dir::String
    parameter_sources::Vector{String}
    seeds::Vector{Int}
    optimizer_reads::Int
    final_reads::Int
    maximum_iterations::Int
    mps_max_bond_dimension::Int
    mps_truncation_threshold::Float64
    threads::Int
end

struct DirectParameterSource
    name::String
    p::Int
    angles::Vector{Float64}
    description::String
    artifact::String
end

function parse_positive_int_env(name::AbstractString, default::Integer)
    raw = strip(get(ENV, name, string(default)))
    value = parse(Int, raw)
    value > 0 || error("$(name) must be a positive integer, got $(raw)")
    return value
end

function parse_nonnegative_int_env(name::AbstractString, default::Integer)
    raw = strip(get(ENV, name, string(default)))
    value = parse(Int, raw)
    value >= 0 || error("$(name) must be nonnegative, got $(raw)")
    return value
end

function parse_int_list(value::AbstractString)
    items = Int[]
    for part in split(value, ',')
        text = strip(part)
        isempty(text) && continue
        if occursin(':', text)
            bounds = split(text, ':')
            length(bounds) == 2 || error("Invalid integer range: $(text)")
            first = parse(Int, strip(bounds[1]))
            last = parse(Int, strip(bounds[2]))
            first <= last || error("Invalid descending integer range: $(text)")
            append!(items, first:last)
        else
            push!(items, parse(Int, text))
        end
    end
    isempty(items) && error("At least one seed is required")
    return items
end

function parse_string_list(value::AbstractString)
    items = [strip(part) for part in split(value, ',') if !isempty(strip(part))]
    isempty(items) && error("At least one parameter source is required")
    return items
end

function read_highread_config()
    output_name = strip(get(ENV, "DSMFG_DIRECT_HIGHREAD_OUTPUT_DIR", DEFAULT_HIGHREAD_OUTPUT_DIR))
    isempty(output_name) && error("DSMFG_DIRECT_HIGHREAD_OUTPUT_DIR must not be empty")
    output_dir = isabspath(output_name) ? output_name : joinpath(STUDY_ROOT, output_name)
    threads = min(Sys.CPU_THREADS, parse_positive_int_env("DSMFG_DIRECT_HIGHREAD_THREADS", 8))

    return DirectHighreadConfig(
        normpath(abspath(output_dir)),
        parse_string_list(get(ENV, "DSMFG_DIRECT_HIGHREAD_PARAMETER_SOURCES", "direct_p2")),
        parse_int_list(get(ENV, "DSMFG_DIRECT_HIGHREAD_SEEDS", "96001")),
        parse_positive_int_env("DSMFG_DIRECT_HIGHREAD_OPTIMIZER_READS", 32),
        parse_positive_int_env("DSMFG_DIRECT_HIGHREAD_FINAL_READS", 8192),
        parse_nonnegative_int_env("DSMFG_DIRECT_HIGHREAD_MAX_ITER", 0),
        parse_positive_int_env("DSMFG_DIRECT_HIGHREAD_MPS_MAX_BOND_DIMENSION", 64),
        parse(Float64, get(ENV, "DSMFG_DIRECT_HIGHREAD_MPS_TRUNCATION_THRESHOLD", "1.0e-10")),
        max(1, threads),
    )
end

function parse_angles(value::AbstractString)
    angles = [parse(Float64, strip(part)) for part in split(value, ';') if !isempty(strip(part))]
    iseven(length(angles)) || error("Expected beta-then-gamma QAOA angles, got an odd count")
    return angles
end

function angle_row(path::AbstractString, target::AbstractString, p::Integer)
    for row in read_csv_rows(path)
        row["target"] == target && parse(Int, row["p"]) == p && return row
    end
    error("Missing $(target) p=$(p) angle row in $(path)")
end

function source_from_name(name::AbstractString, qubo::QuboData)
    if name == "direct_p2"
        artifact = joinpath(STUDY_ROOT, DEFAULT_OUTPUT_DIR, DIRECT_QAOA_PARAMETER_ARTIFACT)
        params = direct_qaoa_parameter_artifact(qubo)
        return DirectParameterSource(
            name,
            params["p"],
            Float64.(params["qiskit_angles_beta_then_gamma"]),
            params["angle_source"],
            relpath(artifact, STUDY_ROOT),
        )
    end

    if startswith(name, "reduced_")
        parts = split(name, '_')
        length(parts) == 3 || error("Reduced parameter source must look like reduced_top10_p5, got $(name)")
        target = parts[2]
        p_text = parts[3]
        startswith(p_text, "p") || error("Reduced parameter source must include p depth, got $(name)")
        p = parse(Int, p_text[2:end])
        angle_dir = target == "global" ?
            joinpath(STUDY_ROOT, "ds_mfg_qaoa_juliqaoa_global_angle_search_p3") :
            joinpath(STUDY_ROOT, "ds_mfg_qaoa_juliqaoa_objective_angle_search")
        path = joinpath(angle_dir, "juliqaoa_angle_summary.csv")
        row = angle_row(path, target, p)
        angles = parse_angles(row["qiskit_angles_beta_then_gamma"])
        length(angles) == 2p || error("Expected $(2p) angles for $(name), got $(length(angles))")
        return DirectParameterSource(
            name,
            p,
            angles,
            "reduced-surrogate JuliQAOA $(target) p=$(p) transferred onto direct full-QUBO cost Hamiltonian",
            relpath(path, STUDY_ROOT),
        )
    end

    error("Unknown direct high-read parameter source: $(name)")
end

function build_direct_qubo_model(qubo::QuboData)
    model = MOI.instantiate(QAOA.Optimizer; with_cache_type = Float64)
    x = MOI.add_variables(model, qubo.n)
    for variable in x
        MOI.add_constraint(model, variable, MOI.ZeroOne())
    end

    quadratic_terms = MOI.ScalarQuadraticTerm{Float64}[]
    for i in 1:qubo.n, j in 1:qubo.n
        coefficient = qubo.scale * qubo.quadratic[i, j]
        iszero(coefficient) && continue
        push!(quadratic_terms, MOI.ScalarQuadraticTerm(i == j ? 2coefficient : coefficient, x[i], x[j]))
    end

    linear_terms = MOI.ScalarAffineTerm{Float64}[]
    for i in 1:qubo.n
        coefficient = qubo.scale * qubo.linear[i]
        iszero(coefficient) && continue
        push!(linear_terms, MOI.ScalarAffineTerm(coefficient, x[i]))
    end

    objective = MOI.ScalarQuadraticFunction(quadratic_terms, linear_terms, qubo.scale * qubo.offset)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{typeof(objective)}(), objective)
    return model, x
end

function configure_direct_qaoa!(
    model,
    source::DirectParameterSource,
    seed::Integer,
    config::DirectHighreadConfig,
)
    MOI.set(model, QAOA.NumberOfLayers(), source.p)
    MOI.set(model, QAOA.NumberOfReads(), config.optimizer_reads)
    MOI.set(model, QUBODrivers.RandomSeed(), seed)
    MOI.set(model, QUBODrivers.FinalNumberOfReads(), config.final_reads)
    MOI.set(model, QAOA.MaximumIterations(), config.maximum_iterations)
    MOI.set(model, QAOA.InitialParameters(), source.angles)
    MOI.set(model, QAOA.InitialParameterSource(), source.description)
    MOI.set(model, QAOA.AerBackendMethod(), "matrix_product_state")
    MOI.set(model, QAOA.AerPrecision(), "single")
    MOI.set(model, QAOA.AerMaxParallelThreads(), config.threads)
    MOI.set(model, QAOA.AerMPSOmpThreads(), config.threads)
    MOI.set(model, QAOA.AerMPSTruncationThreshold(), config.mps_truncation_threshold)
    MOI.set(model, QAOA.AerMPSMaxBondDimension(), config.mps_max_bond_dimension)
    MOI.set(model, QAOA.AerMPSSampleMeasureAlgorithm(), "mps_heuristic")
    return model
end

function summary_header()
    return [
        "algorithm", "parameter_source", "parameter_artifact", "optimized_parameter_artifact",
        "p", "seed", "optimizer_reads", "final_reads", "maximum_iterations", "solve_time_sec",
        "unique_states", "total_reads", "top50_hits", "top10_hits", "global_hits",
        "gurobi_pool_feasible_hits", hit_rate_stat_headers()..., "best_top50_rank",
        "best_top50_exact_repaired_qubo_energy", "best_top50_match", "best_top50_flow_bits",
        "best_repaired_qubo_energy", "best_repaired_match", "best_repaired_flow_bits",
        "distribution_csv", "parameter_order", "optimized_qiskit_angles_beta_then_gamma",
    ]
end

function distribution_header()
    return [
        "algorithm", "parameter_source", "p", "seed", "result", "reads", "probability",
        "qubo_energy", "top50_rank", "exact_repaired_qubo_energy", "match",
        "matched_ip_obj", "flow_bits", "full_bits", "repaired_aux_bits",
        "repaired_full_bits", "encoded_repaired_match",
    ]
end

function full_qubo_energy(bits::Vector{Int}, qubo::QuboData)
    length(bits) == qubo.n || error("Expected $(qubo.n) bits, got $(length(bits))")
    total = qubo.offset
    for i in eachindex(bits)
        bits[i] == 1 && (total += qubo.linear[i])
    end
    for i in eachindex(bits), j in eachindex(bits)
        bits[i] == 1 && bits[j] == 1 && (total += qubo.quadratic[i, j])
    end
    return qubo.scale * total
end

function sampler_count_rows(py_counts, n::Integer)
    rows = Vector{Tuple{String,Vector{Int},Int}}()
    for key in py_counts.keys()
        key_text = pyconvert(String, key)
        bits = Int.(QAOA.count_key_bits(key_text))
        length(bits) == n || error("Expected $(n) bits from Qiskit count key $(key_text)")
        push!(rows, (key_text, bits, pyconvert(Int, py_counts[key])))
    end
    sort!(rows; by = row -> bits_signature(row[2]))
    return rows
end

function write_optimized_parameter_artifact(
    path::AbstractString,
    source::DirectParameterSource,
    seed::Integer,
    config::DirectHighreadConfig,
    qubo::QuboData,
    optimized_angles::Vector{Float64},
    summary,
    distribution_path::AbstractString,
)
    write_json_file(
        path,
        Dict{String,Any}(
            "schema_version" => 1,
            "artifact_role" => "direct_full_qubo_qaoa_hardware_parameter_handoff",
            "algorithm" => "QAOA_direct_full_qubo_optimized_parameter_handoff",
            "problem_name" => "DS-MFG direct original 36-variable QUBO",
            "data_source" => ZIP_NAME,
            "n_qubits" => qubo.n,
            "n_flow_variables" => 19,
            "n_auxiliary_variables" => qubo.n - 19,
            "scale" => qubo.scale,
            "offset" => qubo.offset,
            "angle_source" => "optimized direct full-QUBO p=$(source.p) QAOA from $(source.name)",
            "initial_angle_source" => source.description,
            "initial_parameter_artifact" => source.artifact,
            "source_summary_artifact" => relpath(joinpath(config.output_dir, "direct_full_qubo_qaoa_highread_summary.csv"), STUDY_ROOT),
            "source_distribution_artifact" => relpath(distribution_path, STUDY_ROOT),
            "p" => source.p,
            "parameter_order" => "beta_then_gamma",
            "qiskit_angles_beta_then_gamma" => optimized_angles,
            "initial_qiskit_angles_beta_then_gamma" => source.angles,
            "optimizer" => Dict{String,Any}(
                "maximum_iterations" => config.maximum_iterations,
                "optimizer_reads" => config.optimizer_reads,
                "seed" => seed,
                "backend_method" => "matrix_product_state",
                "mps_truncation_threshold" => config.mps_truncation_threshold,
                "mps_max_bond_dimension" => config.mps_max_bond_dimension,
            ),
            "validation_sample" => Dict{String,Any}(
                "final_reads" => config.final_reads,
                "top50_hits" => summary.top50_hits,
                "top10_hits" => summary.top10_hits,
                "global_hits" => summary.global_hits,
                "gurobi_pool_feasible_hits" => summary.gurobi_pool_feasible_hits,
                "best_repaired_qubo_energy" => summary.best_repaired.exact_repaired_qubo_energy,
                "best_repaired_flow_bits" => summary.best_repaired_flow_bits,
            ),
            "hardware_handoff_policy" => "fixed-parameter sampling only; use the optimized QAOA vector without running an optimizer loop on hardware",
            "measurement_bit_order" => "Qiskit count keys must be reversed before scoring as x1..x36.",
            "producer_script" => relpath(@__FILE__, STUDY_ROOT),
        ),
    )
end

function score_run!(
    summary_io,
    source::DirectParameterSource,
    seed::Integer,
    config::DirectHighreadConfig,
    qubo::QuboData,
    components::Vector{AuxComponent},
    top_flows::Dict{String,NamedTuple},
)
    model, x = build_direct_qubo_model(qubo)
    configure_direct_qaoa!(model, source, seed, config)

    @info "Running direct full-QUBO QAOA high-read sample" source=source.name p=source.p seed final_reads=config.final_reads max_iter=config.maximum_iterations
    MOI.Utilities.attach_optimizer(model)
    sampler = model.optimizer
    optimized_angles = Float64[]
    count_rows = Vector{Tuple{String,Vector{Int},Int}}()
    solve_time = @elapsed QAOA.retrieve(sampler) do result, py_counts
        optimized_angles = pyconvert(Vector{Float64}, result.x)
        count_rows = sampler_count_rows(py_counts, qubo.n)
    end

    result_count = length(count_rows)
    total_reads = sum(row[3] for row in count_rows; init = 0)
    distribution_name = "qaoa_$(source.name)_seed$(seed)_reads$(config.final_reads)_distribution.csv"
    distribution_path = joinpath(config.output_dir, distribution_name)

    top50_hits = 0
    top10_hits = 0
    global_hits = 0
    gurobi_pool_feasible_hits = 0
    best_top50 = nothing
    best_repaired = nothing
    best_repaired_flow_bits = ""

    open(distribution_path, "w") do distribution_io
        println(distribution_io, csv_row(distribution_header()))
        for (result, (_, bits, reads)) in enumerate(count_rows)
            full_bits = bits_signature(bits)
            flow_bits = bits_signature(bits[1:19])
            probability = total_reads == 0 ? 0.0 : reads / total_reads
            qubo_energy = full_qubo_energy(bits, qubo)
            repaired = repair_flow(bits[1:19], qubo, components)
            info = rank_info(flow_bits, top_flows)
            encoded = !ismissing(info.rank) && full_bits == repaired.repaired_full_bits

            if isnothing(best_repaired) ||
                repaired.exact_repaired_qubo_energy < best_repaired.exact_repaired_qubo_energy
                best_repaired = repaired
                best_repaired_flow_bits = flow_bits
            end

            if !ismissing(info.rank)
                top50_hits += reads
                info.rank <= 10 && (top10_hits += reads)
                info.rank == 1 && (global_hits += reads)
                is_gurobi_pool_feasible_match(info.match) && (gurobi_pool_feasible_hits += reads)
                if isnothing(best_top50) || info.rank < best_top50.rank
                    best_top50 = merge(
                        info,
                        (
                            flow_bits = flow_bits,
                            exact_repaired_qubo_energy = repaired.exact_repaired_qubo_energy,
                        ),
                    )
                end
            end

            println(
                distribution_io,
                csv_row([
                    DIRECT_HIGHREAD_ALGORITHM,
                    source.name,
                    source.p,
                    seed,
                    result,
                    reads,
                    probability,
                    qubo_energy,
                    info.rank,
                    repaired.exact_repaired_qubo_energy,
                    info.match,
                    info.matched_ip_obj,
                    flow_bits,
                    full_bits,
                    repaired.repaired_aux_bits,
                    repaired.repaired_full_bits,
                    encoded,
                ]),
            )
        end
    end

    elapsed_text = @sprintf("%.6f", solve_time)
    elapsed_sec = parse(Float64, elapsed_text)
    hit_stats = hit_rate_stat_values(
        total_reads,
        top50_hits,
        top10_hits,
        global_hits,
        gurobi_pool_feasible_hits,
        elapsed_sec,
    )
    best_top50_rank = isnothing(best_top50) ? "" : best_top50.rank
    best_top50_energy = isnothing(best_top50) ? "" : best_top50.exact_repaired_qubo_energy
    best_top50_match = isnothing(best_top50) ? "" : best_top50.match
    best_top50_flow = isnothing(best_top50) ? "" : best_top50.flow_bits
    best_repaired_info = rank_info(best_repaired_flow_bits, top_flows)
    optimized_artifact_name = "qaoa_$(source.name)_seed$(seed)_reads$(config.final_reads)_optimized_parameters.json"
    optimized_artifact_path = joinpath(config.output_dir, optimized_artifact_name)
    run_summary = (
        top50_hits = top50_hits,
        top10_hits = top10_hits,
        global_hits = global_hits,
        gurobi_pool_feasible_hits = gurobi_pool_feasible_hits,
        best_repaired = best_repaired,
        best_repaired_flow_bits = best_repaired_flow_bits,
    )
    write_optimized_parameter_artifact(
        optimized_artifact_path,
        source,
        seed,
        config,
        qubo,
        optimized_angles,
        run_summary,
        distribution_path,
    )
    if best_repaired_info.match == "global_optimum"
        write_optimized_parameter_artifact(
            joinpath(config.output_dir, "direct_full_qubo_qaoa_optimized_parameters.json"),
            source,
            seed,
            config,
            qubo,
            optimized_angles,
            run_summary,
            distribution_path,
        )
    end

    println(
        summary_io,
        csv_row([
            DIRECT_HIGHREAD_ALGORITHM,
            source.name,
            source.artifact,
            relpath(optimized_artifact_path, STUDY_ROOT),
            source.p,
            seed,
            config.optimizer_reads,
            config.final_reads,
            config.maximum_iterations,
            elapsed_text,
            result_count,
            total_reads,
            top50_hits,
            top10_hits,
            global_hits,
            gurobi_pool_feasible_hits,
            hit_stats...,
            best_top50_rank,
            best_top50_energy,
            best_top50_match,
            best_top50_flow,
            isnothing(best_repaired) ? "" : best_repaired.exact_repaired_qubo_energy,
            best_repaired_info.match,
            best_repaired_flow_bits,
            relpath(distribution_path, STUDY_ROOT),
            "beta_then_gamma",
            join(optimized_angles, ';'),
        ]),
    )
    flush(summary_io)

    return (
        top50_hits = top50_hits,
        top10_hits = top10_hits,
        global_hits = global_hits,
        best_repaired = best_repaired,
        best_repaired_flow_bits = best_repaired_flow_bits,
        distribution_path = distribution_path,
    )
end

function main()
    config = read_highread_config()
    mkpath(config.output_dir)
    qubo = read_qubo_data(STUDY_ROOT)
    components = read_aux_components(joinpath(REDUCED_DIR, "auxiliary_components.csv"), qubo)
    top_flows = read_exact_top_flows(joinpath(REDUCED_DIR, "reduced_exact_top_flows.csv"))
    sources = [source_from_name(name, qubo) for name in config.parameter_sources]
    summary_path = joinpath(config.output_dir, "direct_full_qubo_qaoa_highread_summary.csv")

    open(summary_path, "w") do summary_io
        println(summary_io, csv_row(summary_header()))
        for source in sources, seed in config.seeds
            result = score_run!(summary_io, source, seed, config, qubo, components, top_flows)
            @info "Finished direct full-QUBO high-read sample" source=source.name seed global_hits=result.global_hits top10_hits=result.top10_hits top50_hits=result.top50_hits best_repaired=result.best_repaired.exact_repaired_qubo_energy best_flow=result.best_repaired_flow_bits
        end
    end

    println("Direct full-QUBO QAOA high-read summary written to: ", summary_path)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
