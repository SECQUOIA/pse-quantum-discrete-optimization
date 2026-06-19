#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Printf
using QiskitOpt
using QiskitOpt: QAOA
using QUBODrivers
using QUBOTools

const MOI = QUBODrivers.MOI
const SampleReads = QUBOTools.__moi_num_reads()

include(joinpath(@__DIR__, "hit_rate_stats.jl"))

function parse_int_list(value::AbstractString)
    items = Int[]
    for part in split(value, ',')
        text = strip(part)
        isempty(text) && continue
        if occursin(':', text)
            bounds = split(text, ':')
            length(bounds) == 2 || error("Invalid integer range: $(text)")
            append!(items, parse(Int, bounds[1]):parse(Int, bounds[2]))
        else
            push!(items, parse(Int, text))
        end
    end
    return items
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
    values = Float64[]
    for row in readlines(path)[2:end]
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

function csv_value(value)
    ismissing(value) && return ""
    return string(value)
end

bits_signature(bits) = join(string.(bits))

function parse_angle_summary(path)
    rows = readlines(path)
    header = split(rows[1], ',')
    column = Dict(name => index for (index, name) in enumerate(header))
    parsed = Dict{Tuple{String,Int},Dict{String,String}}()
    for row in rows[2:end]
        isempty(strip(row)) && continue
        values = split(row, ',')
        record = Dict(name => get(values, index, "") for (name, index) in column)
        parsed[(record["target"], parse(Int, record["p"]))] = record
    end
    return parsed
end

function parse_angles(value)
    return [parse(Float64, part) for part in split(value, ';') if !isempty(strip(part))]
end

function parse_selection(value)
    selections = Tuple{String,Int}[]
    for item in split(value, ',')
        text = strip(item)
        isempty(text) && continue
        parts = split(text, ':')
        length(parts) == 2 || error("Invalid selection $(text). Use target:p, e.g. energy:5.")
        push!(selections, (strip(parts[1]), parse(Int, parts[2])))
    end
    return selections
end

function build_qubo_model(L, Q; scale = 1.0, offset = 0.0)
    model = MOI.instantiate(QAOA.Optimizer; with_cache_type = Float64)
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

function configure_qaoa!(model; p, seed, optimizer_reads, final_reads, maximum_iterations, initial_parameters, source)
    threads = max(1, min(Sys.CPU_THREADS, 8))
    MOI.set(model, QAOA.NumberOfLayers(), p)
    MOI.set(model, QAOA.NumberOfReads(), optimizer_reads)
    MOI.set(model, QUBODrivers.RandomSeed(), seed)
    MOI.set(model, QUBODrivers.FinalNumberOfReads(), final_reads)
    MOI.set(model, QAOA.MaximumIterations(), maximum_iterations)
    MOI.set(model, QAOA.InitialParameters(), initial_parameters)
    MOI.set(model, QAOA.InitialParameterSource(), source)
    MOI.set(model, QAOA.AerBackendMethod(), "matrix_product_state")
    MOI.set(model, QAOA.AerPrecision(), "single")
    MOI.set(model, QAOA.AerMaxParallelThreads(), threads)
    MOI.set(model, QAOA.AerMPSOmpThreads(), threads)
    MOI.set(model, QAOA.AerMPSTruncationThreshold(), 1.0e-10)
    MOI.set(model, QAOA.AerMPSMaxBondDimension(), 64)
    MOI.set(model, QAOA.AerMPSSampleMeasureAlgorithm(), "mps_heuristic")
    return model
end

function main()
    base = abspath(joinpath(@__DIR__, ".."))
    reduced_dir = joinpath(base, "ds_mfg_reduced_flow_objective")
    angle_dir = joinpath(base, get(ENV, "DSMFG_QAOA_JULIQAOA_ANGLE_DIR", "ds_mfg_qaoa_juliqaoa_angle_search"))
    angle_path = joinpath(angle_dir, "juliqaoa_angle_summary.csv")
    output_dir = joinpath(base, get(ENV, "DSMFG_QAOA_JULIQAOA_OUTPUT_DIR", "ds_mfg_qaoa_juliqaoa_transfer"))
    mkpath(output_dir)

    optimizer_reads = parse(Int, get(ENV, "DSMFG_QAOA_JULIQAOA_OPTIMIZER_READS", "32"))
    final_reads = parse(Int, get(ENV, "DSMFG_QAOA_JULIQAOA_FINAL_READS", "32768"))
    maximum_iterations = parse(Int, get(ENV, "DSMFG_QAOA_JULIQAOA_MAX_ITER", "0"))
    seeds = parse_int_list(get(ENV, "DSMFG_QAOA_JULIQAOA_SEEDS", "92001"))
    selections = parse_selection(get(ENV, "DSMFG_QAOA_JULIQAOA_SELECTIONS", "energy:5"))

    scalars = read_reduced_scalars(joinpath(reduced_dir, "reduced_scalars.csv"))
    L = read_reduced_vector(joinpath(reduced_dir, "reduced_L_vector.csv"))
    Q = read_reduced_matrix(joinpath(reduced_dir, "reduced_Q_matrix.csv"))
    top_flows = read_exact_top_flows(joinpath(reduced_dir, "reduced_exact_top_flows.csv"))
    angle_rows = parse_angle_summary(angle_path)

    summary_path = joinpath(output_dir, "qaoa_juliqaoa_transfer_summary.csv")
    hits_path = joinpath(output_dir, "qaoa_juliqaoa_transfer_hits.csv")

    summary_header = [
        "algorithm", "angle_target", "p", "seed", "optimizer_reads", "final_reads",
        "maximum_iterations", "solve_time_sec", "unique_states", "total_reads",
        "top50_hits", "top10_hits", "global_hits", "gurobi_pool_feasible_hits", hit_rate_stat_headers()..., "best_top50_rank",
        "best_top50_exact_repaired_qubo_energy", "best_top50_match", "best_top50_flow_bits",
        "juliqaoa_top50_probability", "juliqaoa_top10_probability", "juliqaoa_global_probability",
    ]
    hits_header = [
        "algorithm", "angle_target", "p", "seed", "result", "reads", "probability",
        "surrogate_qubo_energy", "top50_rank", "exact_repaired_qubo_energy", "match",
        "matched_ip_obj", "flow_bits", "repaired_aux_bits", "repaired_full_bits",
    ]

    open(summary_path, "w") do summary_io
        println(summary_io, join(summary_header, ','))
        open(hits_path, "w") do hits_io
            println(hits_io, join(hits_header, ','))

            for (target, p) in selections, seed in seeds
                angle_record = angle_rows[(target, p)]
                initial_parameters = parse_angles(angle_record["qiskit_angles_beta_then_gamma"])
                source = "JuliQAOA $(target) p=$(p) seed=$(angle_record["seed"]) niter=$(angle_record["basinhopping_niter"])"

                @info "Running QiskitOpt QAOA with JuliQAOA angles" target p seed optimizer_reads final_reads maximum_iterations
                model, x = build_qubo_model(L, Q; scale = scalars.scale, offset = scalars.offset)
                configure_qaoa!(
                    model;
                    p = p,
                    seed = seed,
                    optimizer_reads = optimizer_reads,
                    final_reads = final_reads,
                    maximum_iterations = maximum_iterations,
                    initial_parameters = initial_parameters,
                    source = source,
                )
                MOI.optimize!(model)

                result_count = MOI.get(model, MOI.ResultCount())
                total_reads = sum(MOI.get(model, SampleReads(result)) for result in 1:result_count; init = 0)
                solve_time = MOI.get(model, MOI.SolveTimeSec())
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
                                "QAOA_reduced_surrogate_JuliQAOA_transfer", target, p, seed, result,
                                reads, probability, surrogate_energy, hit.rank,
                                hit.exact_repaired_qubo_energy, hit.match, hit.matched_ip_obj,
                                flow_bits, hit.repaired_aux_bits, hit.repaired_full_bits,
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
                        "QAOA_reduced_surrogate_JuliQAOA_transfer", target, p, seed,
                        optimizer_reads, final_reads, maximum_iterations, solve_time,
                        result_count, total_reads, top50_hits, top10_hits, global_hits, gurobi_pool_feasible_hits,
                        hit_stats..., "", "", "", "", angle_record["top50_probability"],
                        angle_record["top10_probability"], angle_record["global_probability"],
                    ]
                else
                    Any[
                        "QAOA_reduced_surrogate_JuliQAOA_transfer", target, p, seed,
                        optimizer_reads, final_reads, maximum_iterations, solve_time,
                        result_count, total_reads, top50_hits, top10_hits, global_hits, gurobi_pool_feasible_hits,
                        hit_stats..., best_hit.rank, best_hit.exact_repaired_qubo_energy, best_hit.match,
                        best_hit.flow_bits, angle_record["top50_probability"],
                        angle_record["top10_probability"], angle_record["global_probability"],
                    ]
                end
                println(summary_io, join(csv_value.(summary_row), ','))
                flush(summary_io)
                flush(hits_io)
                if isnothing(best_hit)
                    @info "No top-50 hit" target p seed
                else
                    @info "Top-50 hit" target p seed rank=best_hit.rank energy=best_hit.exact_repaired_qubo_energy global_hits
                end
            end
        end
    end

    println("Summary written to: ", summary_path)
    println("Hits written to: ", hits_path)
end

main()
