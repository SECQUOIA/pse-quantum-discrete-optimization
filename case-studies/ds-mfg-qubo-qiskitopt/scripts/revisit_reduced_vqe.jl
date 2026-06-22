#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using DelimitedFiles
using Printf
using QiskitOpt
using QiskitOpt: VQE
using QUBODrivers
using QUBOTools

const MOI = QUBODrivers.MOI
const SampleReads = QUBOTools.__moi_num_reads()

include(joinpath(@__DIR__, "hit_rate_stats.jl"))

function parse_int_list(value::AbstractString)
    parts = split(value, ',')
    seeds = Int[]
    for part in parts
        text = strip(part)
        isempty(text) && continue
        if occursin(':', text)
            bounds = split(text, ':')
            length(bounds) == 2 || error("Invalid range in seed list: $(text)")
            append!(seeds, parse(Int, bounds[1]):parse(Int, bounds[2]))
        else
            push!(seeds, parse(Int, text))
        end
    end
    return seeds
end

function read_reduced_scalars(path)
    rows = readlines(path)
    values = split(rows[2], ',')
    return (
        n = round(Int, parse(Float64, values[1])),
        scale = parse(Float64, values[2]),
        offset = parse(Float64, values[3]),
    )
end

function read_reduced_vector(path)
    rows = readlines(path)[2:end]
    values = Float64[]
    for row in rows
        parts = split(row, ',')
        push!(values, parse(Float64, parts[2]))
    end
    return values
end

function read_reduced_matrix(path)
    rows = readlines(path)[2:end]
    matrix = Matrix{Float64}(undef, length(rows), length(rows))
    for (i, row) in enumerate(rows)
        parts = split(row, ',')
        for j in 1:length(rows)
            matrix[i, j] = parse(Float64, parts[j + 1])
        end
    end
    return matrix
end

function read_exact_top_flows(path)
    rows = readlines(path)[2:end]
    top = Dict{String,NamedTuple}()
    for row in rows
        parts = split(row, ',')
        top[parts[3]] = (
            rank = parse(Int, parts[1]),
            exact_repaired_qubo_energy = parse(Float64, parts[2]),
            repaired_aux_bits = parts[4],
            repaired_full_bits = parts[5],
            match = parts[6],
            matched_ip_obj = length(parts) >= 7 && !isempty(parts[7]) ? parse(Float64, parts[7]) : missing,
        )
    end
    return top
end

bits_signature(bits) = join(string.(bits))

function build_qubo_model(L, Q; scale = 1.0, offset = 0.0)
    model = MOI.instantiate(VQE.Optimizer; with_cache_type = Float64)
    n = length(L)
    x = MOI.add_variables(model, n)
    for variable in x
        MOI.add_constraint(model, variable, MOI.ZeroOne())
    end

    quadratic_terms = MOI.ScalarQuadraticTerm{Float64}[]
    for i in 1:n, j in 1:n
        coefficient = scale * Q[i, j]
        iszero(coefficient) && continue
        push!(quadratic_terms, MOI.ScalarQuadraticTerm(i == j ? 2coefficient : coefficient, x[i], x[j]))
    end

    linear_terms = MOI.ScalarAffineTerm{Float64}[]
    for i in 1:n
        coefficient = scale * L[i]
        iszero(coefficient) && continue
        push!(linear_terms, MOI.ScalarAffineTerm(coefficient, x[i]))
    end

    objective = MOI.ScalarQuadraticFunction(quadratic_terms, linear_terms, scale * offset)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(model, MOI.ObjectiveFunction{typeof(objective)}(), objective)
    return model, x
end

function configure_vqe!(model; seed, optimizer_reads, final_reads, maximum_iterations)
    threads = max(1, min(Sys.CPU_THREADS, 8))
    MOI.set(model, VQE.NumberOfReads(), optimizer_reads)
    MOI.set(model, QUBODrivers.RandomSeed(), seed)
    MOI.set(model, QUBODrivers.FinalNumberOfReads(), final_reads)
    MOI.set(model, VQE.MaximumIterations(), maximum_iterations)
    MOI.set(model, VQE.InitialParameters(), VQE.random_initial_parameters(n_variables = 19, seed = seed))
    MOI.set(model, VQE.InitialParameterSource(), "random_seed_$(seed)")
    MOI.set(model, VQE.AerBackendMethod(), "matrix_product_state")
    MOI.set(model, VQE.AerPrecision(), "single")
    MOI.set(model, VQE.AerMaxParallelThreads(), threads)
    MOI.set(model, VQE.AerMPSOmpThreads(), threads)
    MOI.set(model, VQE.AerMPSTruncationThreshold(), 1.0e-10)
    MOI.set(model, VQE.AerMPSMaxBondDimension(), 64)
    MOI.set(model, VQE.AerMPSSampleMeasureAlgorithm(), "mps_heuristic")
    return model
end

function csv_value(value)
    ismissing(value) && return ""
    return string(value)
end

function json_escape(value::AbstractString)
    escaped = replace(value, "\\" => "\\\\", "\"" => "\\\"")
    escaped = replace(escaped, "\n" => "\\n", "\r" => "\\r", "\t" => "\\t")
    return escaped
end

function write_json_value(io::IO, value)
    if value === nothing || ismissing(value)
        print(io, "null")
    elseif value isa Bool
        print(io, value ? "true" : "false")
    elseif value isa AbstractString
        print(io, "\"", json_escape(value), "\"")
    elseif value isa Integer
        print(io, value)
    elseif value isa AbstractFloat
        isfinite(value) ? print(io, @sprintf("%.17g", value)) : print(io, "null")
    elseif value isa AbstractDict
        print(io, "{")
        first_item = true
        for key in sort(collect(keys(value)); by = string)
            first_item || print(io, ",")
            first_item = false
            write_json_value(io, string(key))
            print(io, ":")
            write_json_value(io, value[key])
        end
        print(io, "}")
    elseif value isa Tuple || value isa AbstractVector
        print(io, "[")
        for (index, item) in enumerate(value)
            index == 1 || print(io, ",")
            write_json_value(io, item)
        end
        print(io, "]")
    else
        write_json_value(io, string(value))
    end
end

function write_json_file(path::AbstractString, value)
    mkpath(dirname(path))
    open(path, "w") do io
        write_json_value(io, value)
        println(io)
    end
end

function sampleset_metadata(model)
    raw = MOI.get(model, MOI.RawSolver())
    return deepcopy(QUBOTools.metadata(QUBOTools.solution(raw)))
end

function vqe_metadata_artifact(seed, optimizer_reads, final_reads, maximum_iterations, solve_time, result_count, total_reads, metadata)
    return Dict{String,Any}(
        "schema_version" => 1,
        "artifact_role" => "ds_mfg_reduced_surrogate_vqe_seed_metadata",
        "algorithm" => "VQE_reduced_surrogate_top50_revisit",
        "seed" => seed,
        "configuration" => Dict{String,Any}(
            "optimizer_reads" => optimizer_reads,
            "final_reads" => final_reads,
            "maximum_iterations" => maximum_iterations,
            "ansatz" => "Qiskit EfficientSU2",
            "problem" => "DS-MFG reduced 19-flow quadratic surrogate",
        ),
        "result_summary" => Dict{String,Any}(
            "solve_time_sec" => solve_time,
            "unique_states" => result_count,
            "total_reads" => total_reads,
        ),
        "qiskitopt_sampleset_metadata" => metadata,
        "optimized_parameter_artifact_status" => if haskey(metadata, "optimized_parameters")
            "retained_from_qiskitopt_sampleset_metadata"
        else
            "not_exposed_by_qiskitopt_sampleset_metadata"
        end,
        "credential_policy" => "No IBM tokens, account files, credential paths, or machine-local private output paths are written.",
    )
end

function main()
    base = abspath(joinpath(@__DIR__, ".."))
    reduced_dir = joinpath(base, "ds_mfg_reduced_flow_objective")
    output_dir_name = get(ENV, "DSMFG_REVISIT_OUTPUT_DIR", "ds_mfg_vqe_reduced_flow_objective_v3")
    output_dir = isabspath(output_dir_name) ? output_dir_name : joinpath(base, output_dir_name)
    mkpath(output_dir)

    optimizer_reads = parse(Int, get(ENV, "DSMFG_REVISIT_OPTIMIZER_READS", "128"))
    final_reads = parse(Int, get(ENV, "DSMFG_REVISIT_FINAL_READS", "32768"))
    maximum_iterations = parse(Int, get(ENV, "DSMFG_REVISIT_MAX_ITER", "25"))
    seeds = parse_int_list(get(ENV, "DSMFG_REVISIT_SEEDS", "74001:74020"))

    scalars = read_reduced_scalars(joinpath(reduced_dir, "reduced_scalars.csv"))
    L = read_reduced_vector(joinpath(reduced_dir, "reduced_L_vector.csv"))
    Q = read_reduced_matrix(joinpath(reduced_dir, "reduced_Q_matrix.csv"))
    top_flows = read_exact_top_flows(joinpath(reduced_dir, "reduced_exact_top_flows.csv"))

    summary_path = joinpath(output_dir, "vqe_reduced_top50_sampling_summary.csv")
    hits_path = joinpath(output_dir, "vqe_reduced_top50_sampling_hits.csv")

    summary_header = [
        "algorithm", "seed", "optimizer_reads", "final_reads", "maximum_iterations",
        "solve_time_sec", "unique_states", "total_reads", "top50_hits", "top10_hits",
        "global_hits", "gurobi_pool_feasible_hits", hit_rate_stat_headers()..., "best_top50_rank", "best_top50_exact_repaired_qubo_energy",
        "best_top50_match", "best_top50_flow_bits",
    ]
    hits_header = [
        "algorithm", "seed", "result", "reads", "probability", "surrogate_qubo_energy",
        "top50_rank", "exact_repaired_qubo_energy", "match", "matched_ip_obj",
        "flow_bits", "repaired_aux_bits", "repaired_full_bits",
    ]

    open(summary_path, "w") do summary_io
        println(summary_io, join(summary_header, ','))
        open(hits_path, "w") do hits_io
            println(hits_io, join(hits_header, ','))

            for seed in seeds
                @info "Running reduced-surrogate VQE" seed optimizer_reads final_reads maximum_iterations
                model, x = build_qubo_model(L, Q; scale = scalars.scale, offset = scalars.offset)
                configure_vqe!(model; seed, optimizer_reads, final_reads, maximum_iterations)
                MOI.optimize!(model)

                result_count = MOI.get(model, MOI.ResultCount())
                total_reads = sum(MOI.get(model, SampleReads(result)) for result in 1:result_count; init = 0)
                solve_time = MOI.get(model, MOI.SolveTimeSec())
                metadata = sampleset_metadata(model)
                metadata_path = joinpath(output_dir, "vqe_reduced_seed$(seed)_metadata.json")
                write_json_file(
                    metadata_path,
                    vqe_metadata_artifact(
                        seed,
                        optimizer_reads,
                        final_reads,
                        maximum_iterations,
                        solve_time,
                        result_count,
                        total_reads,
                        metadata,
                    ),
                )
                top50_hits = 0
                top10_hits = 0
                global_hits = 0
                gurobi_pool_feasible_hits = 0
                best_hit = nothing

                for result in 1:result_count
                    bits = round.(Int, MOI.get.(model, MOI.VariablePrimal(result), x))
                    flow_bits = bits_signature(bits)
                    hit = get(top_flows, flow_bits, nothing)
                    isnothing(hit) && continue

                    reads = MOI.get(model, SampleReads(result))
                    surrogate_energy = MOI.get(model, MOI.ObjectiveValue(result))
                    probability = total_reads == 0 ? 0.0 : reads / total_reads
                    top50_hits += reads
                    hit.rank <= 10 && (top10_hits += reads)
                    hit.rank == 1 && (global_hits += reads)
                    is_gurobi_pool_feasible_match(hit.match) && (gurobi_pool_feasible_hits += reads)
                    if isnothing(best_hit) || hit.rank < best_hit.rank
                        best_hit = merge(hit, (flow_bits = flow_bits,))
                    end

                    println(
                        hits_io,
                        join(
                            csv_value.([
                                "VQE_reduced_surrogate_top50_revisit", seed, result, reads, probability,
                                surrogate_energy, hit.rank, hit.exact_repaired_qubo_energy, hit.match,
                                hit.matched_ip_obj, flow_bits, hit.repaired_aux_bits, hit.repaired_full_bits,
                            ]),
                            ',',
                        ),
                    )
                end

                hit_stats = hit_rate_stat_values(
                    total_reads,
                    top50_hits,
                    top10_hits,
                    global_hits,
                    gurobi_pool_feasible_hits,
                    solve_time,
                )
                summary_row = if isnothing(best_hit)
                    Any[
                        "VQE_reduced_surrogate_top50_revisit", seed, optimizer_reads, final_reads,
                        maximum_iterations, solve_time, result_count, total_reads, top50_hits, top10_hits,
                        global_hits, gurobi_pool_feasible_hits, hit_stats..., "", "", "", "",
                    ]
                else
                    Any[
                        "VQE_reduced_surrogate_top50_revisit", seed, optimizer_reads, final_reads,
                        maximum_iterations, solve_time, result_count, total_reads, top50_hits, top10_hits,
                        global_hits, gurobi_pool_feasible_hits, hit_stats..., best_hit.rank, best_hit.exact_repaired_qubo_energy,
                        best_hit.match, best_hit.flow_bits,
                    ]
                end
                println(summary_io, join(csv_value.(summary_row), ','))
                flush(summary_io)
                flush(hits_io)
                if !isnothing(best_hit)
                    @info "Top-50 hit" seed rank=best_hit.rank energy=best_hit.exact_repaired_qubo_energy flow=best_hit.flow_bits global_hits
                else
                    @info "No top-50 hit" seed
                end
            end
        end
    end

    println("Summary written to: ", summary_path)
    println("Hits written to: ", hits_path)
end

main()
