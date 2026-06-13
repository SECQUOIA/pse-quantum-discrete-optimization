#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Printf
using Random

const STUDY_ROOT = abspath(joinpath(@__DIR__, ".."))
const ZIP_NAME = "Fw_ DS mfg case qubo information.zip"

include(joinpath(@__DIR__, "hit_rate_stats.jl"))

struct QuboData
    scale::Float64
    offset::Float64
    linear::Vector{Float64}
    quadratic::Matrix{Float64}
end

struct AuxComponent
    aux_indices::Vector{Int}
    quadratic_terms::Vector{Tuple{Int,Int,Float64}}
end

bits_signature(bits) = join(string.(bits))

function csv_value(value)
    ismissing(value) && return ""
    return string(value)
end

function parse_int_list(value::AbstractString)
    items = Int[]
    for part in split(value, ',')
        text = strip(part)
        isempty(text) && continue
        if occursin(':', text)
            bounds = split(text, ':')
            length(bounds) == 2 || error("Invalid integer range: $(text)")
            append!(items, parse(Int, strip(bounds[1])):parse(Int, strip(bounds[2])))
        else
            push!(items, parse(Int, text))
        end
    end
    return items
end

function read_member_text(base::AbstractString, member::AbstractString)
    candidates = [
        joinpath(base, member),
        joinpath(base, "julia_exports", member),
        joinpath(base, "result_gurobi", member),
    ]
    for path in candidates
        isfile(path) && return read(path, String)
    end

    archive = joinpath(base, ZIP_NAME)
    isfile(archive) || error("Could not find $(member) or $(ZIP_NAME)")
    return read(`unzip -p $archive $member`, String)
end

function gzip_file!(path::AbstractString)
    gzip = Sys.which("gzip")
    isnothing(gzip) && error("gzip is required to compress $(path)")
    run(`$gzip -n -f $path`)
    return path * ".gz"
end

function nonempty_lines(text::AbstractString)
    return [strip(line) for line in split(text, '\n') if !isempty(strip(line))]
end

function parse_scalar_csv(text::AbstractString)
    rows = nonempty_lines(text)
    length(rows) >= 2 || error("Expected header and value row in scalars.csv")
    values = split(rows[2], ',')
    length(values) >= 3 || error("Expected n, scale, offset in scalars.csv")
    return (
        n = round(Int, parse(Float64, strip(values[1]))),
        scale = parse(Float64, strip(values[2])),
        offset = parse(Float64, strip(values[3])),
    )
end

function parse_vector_csv(text::AbstractString)
    return [parse(Float64, strip(row)) for row in nonempty_lines(text)]
end

function parse_matrix_csv(text::AbstractString)
    rows = nonempty_lines(text)
    matrix = Matrix{Float64}(undef, length(rows), length(split(rows[1], ',')))
    for (i, row) in enumerate(rows)
        values = split(row, ',')
        for (j, value) in enumerate(values)
            matrix[i, j] = parse(Float64, strip(value))
        end
    end
    return matrix
end

function read_qubo_data(base::AbstractString)
    scalars = parse_scalar_csv(read_member_text(base, "scalars.csv"))
    linear = parse_vector_csv(read_member_text(base, "L_vector.csv"))
    quadratic = parse_matrix_csv(read_member_text(base, "Q_matrix.csv"))
    length(linear) == scalars.n || error("Expected $(scalars.n) linear coefficients")
    size(quadratic) == (scalars.n, scalars.n) || error("Expected $(scalars.n)x$(scalars.n) Q matrix")
    return QuboData(scalars.scale, scalars.offset, linear, quadratic)
end

function read_aux_components(path::AbstractString, qubo::QuboData; n_flow::Int = 19)
    rows = readlines(path)[2:end]
    components = AuxComponent[]
    for row in rows
        isempty(strip(row)) && continue
        parts = split(row, ',')
        aux_indices = [parse(Int, text) for text in split(parts[2], ';')]
        aux_set = Set(aux_indices)
        terms = Tuple{Int,Int,Float64}[]
        for i in eachindex(qubo.linear), j in eachindex(qubo.linear)
            iszero(qubo.quadratic[i, j]) && continue
            (i in aux_set || j in aux_set) || continue
            i > n_flow && !(i in aux_set) && continue
            j > n_flow && !(j in aux_set) && continue
            push!(terms, (i, j, qubo.quadratic[i, j]))
        end
        push!(components, AuxComponent(aux_indices, terms))
    end
    return components
end

function read_exact_top_flows(path::AbstractString)
    rows = readlines(path)[2:end]
    top = Dict{String,NamedTuple}()
    for row in rows
        isempty(strip(row)) && continue
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

function flow_constant(flow::Vector{Int}, qubo::QuboData)
    total = qubo.offset
    for i in eachindex(flow)
        flow[i] == 1 && (total += qubo.linear[i])
    end
    n_flow = length(flow)
    for i in 1:n_flow, j in 1:n_flow
        flow[i] == 1 && flow[j] == 1 && (total += qubo.quadratic[i, j])
    end
    return total
end

function component_contribution(
    flow::Vector{Int},
    qubo::QuboData,
    component::AuxComponent,
    assignment::Int,
)
    n_flow = length(flow)
    function value_at(index::Int)
        index <= n_flow && return flow[index]
        for (local_index, aux_index) in enumerate(component.aux_indices)
            aux_index == index && return (assignment >> (local_index - 1)) & 1
        end
        return 0
    end

    total = 0.0
    for (local_index, aux_index) in enumerate(component.aux_indices)
        ((assignment >> (local_index - 1)) & 1) == 1 && (total += qubo.linear[aux_index])
    end

    for (i, j, coefficient) in component.quadratic_terms
        value_at(i) == 1 && value_at(j) == 1 && (total += coefficient)
    end
    return total
end

function repair_flow(flow::Vector{Int}, qubo::QuboData, components::Vector{AuxComponent})
    n_flow = length(flow)
    aux_bits = zeros(Int, length(qubo.linear) - n_flow)
    total = flow_constant(flow, qubo)

    for component in components
        best_value = Inf
        best_assignment = 0
        for assignment in 0:(2^length(component.aux_indices)-1)
            contribution = component_contribution(flow, qubo, component, assignment)
            if contribution < best_value
                best_value = contribution
                best_assignment = assignment
            end
        end
        total += best_value
        for (local_index, aux_index) in enumerate(component.aux_indices)
            aux_bits[aux_index - n_flow] = (best_assignment >> (local_index - 1)) & 1
        end
    end

    repaired_energy = qubo.scale * total
    full_bits = vcat(flow, aux_bits)
    return (
        exact_repaired_qubo_energy = repaired_energy,
        repaired_aux_bits = bits_signature(aux_bits),
        repaired_full_bits = bits_signature(full_bits),
    )
end

function flow_from_index(index::Integer, n_flow::Integer)
    return [((index >> (n_flow - i)) & 1) for i in 1:n_flow]
end

function random_flow(rng::AbstractRNG, n_flow::Integer)
    return [rand(rng, Bool) ? 1 : 0 for _ in 1:n_flow]
end

function rank_info(flow_bits::String, top_flows::Dict{String,NamedTuple})
    hit = get(top_flows, flow_bits, nothing)
    isnothing(hit) && return (
        rank = missing,
        match = "not_top50",
        matched_ip_obj = missing,
    )
    return (
        rank = hit.rank,
        match = hit.match,
        matched_ip_obj = hit.matched_ip_obj,
    )
end

function score_flow!(
    distribution::Dict{String,Int},
    repair_cache::Dict{String,NamedTuple},
    flow::Vector{Int},
    qubo::QuboData,
    components::Vector{AuxComponent},
)
    flow_bits = bits_signature(flow)
    distribution[flow_bits] = get(distribution, flow_bits, 0) + 1
    return get!(repair_cache, flow_bits) do
        repair_flow(flow, qubo, components)
    end
end

function uniform_random_distribution(
    rng::AbstractRNG,
    sample_budget::Int,
    n_flow::Int,
    qubo::QuboData,
    components::Vector{AuxComponent},
    repair_cache::Dict{String,NamedTuple},
)
    distribution = Dict{String,Int}()
    for _ in 1:sample_budget
        flow = random_flow(rng, n_flow)
        score_flow!(distribution, repair_cache, flow, qubo, components)
    end
    return distribution, sample_budget, 0
end

function hill_climb_distribution(
    rng::AbstractRNG,
    evaluation_budget::Int,
    n_flow::Int,
    qubo::QuboData,
    components::Vector{AuxComponent},
    repair_cache::Dict{String,NamedTuple};
    max_steps_per_restart::Int,
)
    distribution = Dict{String,Int}()
    evaluations = 0
    restarts = 0

    while evaluations < evaluation_budget
        restarts += 1
        flow = random_flow(rng, n_flow)
        current = score_flow!(distribution, repair_cache, flow, qubo, components)
        evaluations += 1

        for _ in 1:max_steps_per_restart
            evaluations >= evaluation_budget && break
            best_energy = current.exact_repaired_qubo_energy
            best_flow = copy(flow)
            best_record = current

            for bit_index in 1:n_flow
                evaluations >= evaluation_budget && break
                candidate = copy(flow)
                candidate[bit_index] = 1 - candidate[bit_index]
                record = score_flow!(distribution, repair_cache, candidate, qubo, components)
                evaluations += 1
                if record.exact_repaired_qubo_energy < best_energy
                    best_energy = record.exact_repaired_qubo_energy
                    best_flow = candidate
                    best_record = record
                end
            end

            best_flow == flow && break
            flow = best_flow
            current = best_record
        end
    end

    return distribution, evaluations, restarts
end

function summarize_distribution(
    algorithm::String,
    seed::Int,
    sample_budget::Int,
    evaluations::Int,
    restarts::Int,
    elapsed_sec::Float64,
    distribution::Dict{String,Int},
    repair_cache::Dict{String,NamedTuple},
    top_flows::Dict{String,NamedTuple};
    max_steps_per_restart::Union{Missing,Int} = missing,
)
    total_samples = sum(values(distribution); init = 0)
    top50_hits = 0
    top10_hits = 0
    global_hits = 0
    feasible_hits = 0
    best_flow = ""
    best_repair = nothing
    best_rank = missing
    best_match = ""
    best_matched_ip_obj = missing

    for (flow_bits, visits) in distribution
        repair = repair_cache[flow_bits]
        info = rank_info(flow_bits, top_flows)
        !ismissing(info.rank) && (top50_hits += visits)
        !ismissing(info.rank) && info.rank <= 10 && (top10_hits += visits)
        !ismissing(info.rank) && info.rank == 1 && (global_hits += visits)
        is_feasible_match(info.match) && (feasible_hits += visits)

        if isnothing(best_repair) || repair.exact_repaired_qubo_energy < best_repair.exact_repaired_qubo_energy
            best_flow = flow_bits
            best_repair = repair
            best_rank = info.rank
            best_match = info.match
            best_matched_ip_obj = info.matched_ip_obj
        end
    end

    settings = ismissing(max_steps_per_restart) ?
        "independent_uniform_samples" :
        "steepest_descent_restarts;max_steps_per_restart=$(max_steps_per_restart)"
    hit_stats = hit_rate_stat_values(total_samples, top50_hits, top10_hits, global_hits, feasible_hits, elapsed_sec)

    return Any[
        algorithm,
        seed,
        sample_budget,
        evaluations,
        restarts,
        @sprintf("%.6f", elapsed_sec),
        length(distribution),
        total_samples,
        top50_hits,
        top10_hits,
        global_hits,
        feasible_hits,
        hit_stats...,
        best_rank,
        isnothing(best_repair) ? "" : best_repair.exact_repaired_qubo_energy,
        best_match,
        best_flow,
        isnothing(best_repair) ? "" : best_repair.repaired_aux_bits,
        isnothing(best_repair) ? "" : best_repair.repaired_full_bits,
        best_matched_ip_obj,
        settings,
    ]
end

function write_full_distribution!(
    io,
    algorithm::String,
    seed::Int,
    sample_budget::Int,
    distribution::Dict{String,Int},
    repair_cache::Dict{String,NamedTuple},
    top_flows::Dict{String,NamedTuple},
)
    total_samples = sum(values(distribution); init = 0)

    for flow_bits in sort(collect(keys(distribution)))
        info = rank_info(flow_bits, top_flows)
        visits = distribution[flow_bits]
        repair = repair_cache[flow_bits]
        println(
            io,
            join(
                csv_value.([
                    algorithm,
                    seed,
                    sample_budget,
                    flow_bits,
                    visits,
                    total_samples == 0 ? 0.0 : visits / total_samples,
                    info.rank,
                    repair.exact_repaired_qubo_energy,
                    info.match,
                    info.matched_ip_obj,
                    repair.repaired_aux_bits,
                    repair.repaired_full_bits,
                ]),
                ',',
            ),
        )
    end
end

function write_retained_flows!(
    io,
    algorithm::String,
    seed::Int,
    sample_budget::Int,
    distribution::Dict{String,Int},
    repair_cache::Dict{String,NamedTuple},
    top_flows::Dict{String,NamedTuple},
)
    total_samples = sum(values(distribution); init = 0)
    best_energy = Inf
    best_flow_bits = ""
    for (flow_bits, _) in distribution
        repair = repair_cache[flow_bits]
        if repair.exact_repaired_qubo_energy < best_energy
            best_energy = repair.exact_repaired_qubo_energy
            best_flow_bits = flow_bits
        end
    end

    for flow_bits in sort(collect(keys(distribution)))
        info = rank_info(flow_bits, top_flows)
        retain_reasons = String[]
        !ismissing(info.rank) && push!(retain_reasons, "top50")
        flow_bits == best_flow_bits && push!(retain_reasons, "best_sampled_flow")
        isempty(retain_reasons) && continue

        visits = distribution[flow_bits]
        repair = repair_cache[flow_bits]
        println(
            io,
            join(
                csv_value.([
                    algorithm,
                    seed,
                    sample_budget,
                    flow_bits,
                    visits,
                    total_samples == 0 ? 0.0 : visits / total_samples,
                    info.rank,
                    repair.exact_repaired_qubo_energy,
                    info.match,
                    info.matched_ip_obj,
                    repair.repaired_aux_bits,
                    repair.repaired_full_bits,
                    join(retain_reasons, ';'),
                ]),
                ',',
            ),
        )
    end
end

function validate_repair(
    qubo::QuboData,
    components::Vector{AuxComponent},
    top_flows::Dict{String,NamedTuple},
)
    for (flow_bits, hit) in top_flows
        flow = [parse(Int, char) for char in collect(flow_bits)]
        repair = repair_flow(flow, qubo, components)
        isapprox(
            repair.exact_repaired_qubo_energy,
            hit.exact_repaired_qubo_energy;
            atol = 1.0e-8,
            rtol = 1.0e-10,
        ) || error(
            "Repair validation failed for $(flow_bits): got $(repair.exact_repaired_qubo_energy), expected $(hit.exact_repaired_qubo_energy)",
        )
    end
end

function run_baseline!(
    summary_io,
    distribution_io,
    retained_io,
    algorithm::String,
    seed::Int,
    sample_budget::Int,
    qubo::QuboData,
    components::Vector{AuxComponent},
    top_flows::Dict{String,NamedTuple},
    repair_cache::Dict{String,NamedTuple};
    max_steps_per_restart::Union{Missing,Int} = missing,
)
    rng = MersenneTwister(seed)
    start_ns = time_ns()

    distribution, evaluations, restarts = if algorithm == "uniform_random_repaired_flow"
        uniform_random_distribution(rng, sample_budget, 19, qubo, components, repair_cache)
    elseif algorithm == "hill_climb_restarts_repaired_flow"
        hill_climb_distribution(
            rng,
            sample_budget,
            19,
            qubo,
            components,
            repair_cache;
            max_steps_per_restart = max_steps_per_restart,
        )
    else
        error("Unknown algorithm $(algorithm)")
    end

    elapsed_sec = (time_ns() - start_ns) / 1.0e9
    summary_row = summarize_distribution(
        algorithm,
        seed,
        sample_budget,
        evaluations,
        restarts,
        elapsed_sec,
        distribution,
        repair_cache,
        top_flows;
        max_steps_per_restart = max_steps_per_restart,
    )
    println(summary_io, join(csv_value.(summary_row), ','))
    write_full_distribution!(distribution_io, algorithm, seed, sample_budget, distribution, repair_cache, top_flows)
    write_retained_flows!(retained_io, algorithm, seed, sample_budget, distribution, repair_cache, top_flows)
    flush(summary_io)
    flush(distribution_io)
    flush(retained_io)
    @info "Completed classical baseline" algorithm seed sample_budget elapsed_sec
end

function main()
    reduced_dir = joinpath(STUDY_ROOT, "ds_mfg_reduced_flow_objective")
    output_dir_name = get(ENV, "DSMFG_CLASSICAL_OUTPUT_DIR", "ds_mfg_classical_baselines")
    output_dir = isabspath(output_dir_name) ? output_dir_name : joinpath(STUDY_ROOT, output_dir_name)
    mkpath(output_dir)

    qubo = read_qubo_data(STUDY_ROOT)
    components = read_aux_components(joinpath(reduced_dir, "auxiliary_components.csv"), qubo)
    top_flows = read_exact_top_flows(joinpath(reduced_dir, "reduced_exact_top_flows.csv"))
    validate_repair(qubo, components, top_flows)

    uniform_budgets = parse_int_list(get(ENV, "DSMFG_CLASSICAL_UNIFORM_BUDGETS", "262144,524288"))
    uniform_seeds = parse_int_list(get(ENV, "DSMFG_CLASSICAL_UNIFORM_SEEDS", "81001,81002"))
    length(uniform_seeds) == length(uniform_budgets) ||
        error("DSMFG_CLASSICAL_UNIFORM_SEEDS must match DSMFG_CLASSICAL_UNIFORM_BUDGETS")

    hill_budgets = parse_int_list(get(ENV, "DSMFG_CLASSICAL_HILL_BUDGETS", "262144"))
    hill_seeds = parse_int_list(get(ENV, "DSMFG_CLASSICAL_HILL_SEEDS", "82001"))
    length(hill_seeds) == length(hill_budgets) ||
        error("DSMFG_CLASSICAL_HILL_SEEDS must match DSMFG_CLASSICAL_HILL_BUDGETS")
    hill_max_steps = parse(Int, get(ENV, "DSMFG_CLASSICAL_HILL_MAX_STEPS_PER_RESTART", "128"))

    summary_path = joinpath(output_dir, "classical_baseline_summary.csv")
    distribution_path = joinpath(output_dir, "classical_baseline_distribution.csv")
    retained_path = joinpath(output_dir, "classical_baseline_retained_flows.csv")

    summary_header = [
        "algorithm",
        "seed",
        "sample_budget",
        "objective_evaluations",
        "restart_count",
        "wall_time_sec",
        "unique_flows",
        "total_samples",
        "top50_hits",
        "top10_hits",
        "global_hits",
        "feasible_hits",
        hit_rate_stat_headers()...,
        "best_top50_rank",
        "best_repaired_qubo_energy",
        "best_match",
        "best_flow_bits",
        "best_repaired_aux_bits",
        "best_repaired_full_bits",
        "best_matched_ip_obj",
        "settings",
    ]
    distribution_header = [
        "algorithm",
        "seed",
        "sample_budget",
        "flow_bits",
        "visits",
        "probability",
        "top50_rank",
        "exact_repaired_qubo_energy",
        "match",
        "matched_ip_obj",
        "repaired_aux_bits",
        "repaired_full_bits",
    ]
    retained_header = [
        "algorithm",
        "seed",
        "sample_budget",
        "flow_bits",
        "visits",
        "probability",
        "top50_rank",
        "exact_repaired_qubo_energy",
        "match",
        "matched_ip_obj",
        "repaired_aux_bits",
        "repaired_full_bits",
        "retained_reason",
    ]

    repair_cache = Dict{String,NamedTuple}()

    open(summary_path, "w") do summary_io
        println(summary_io, join(summary_header, ','))
        open(distribution_path, "w") do distribution_io
            println(distribution_io, join(distribution_header, ','))
            open(retained_path, "w") do retained_io
                println(retained_io, join(retained_header, ','))

                for (seed, budget) in zip(uniform_seeds, uniform_budgets)
                    run_baseline!(
                        summary_io,
                        distribution_io,
                        retained_io,
                        "uniform_random_repaired_flow",
                        seed,
                        budget,
                        qubo,
                        components,
                        top_flows,
                        repair_cache,
                    )
                end

                for (seed, budget) in zip(hill_seeds, hill_budgets)
                    run_baseline!(
                        summary_io,
                        distribution_io,
                        retained_io,
                        "hill_climb_restarts_repaired_flow",
                        seed,
                        budget,
                        qubo,
                        components,
                        top_flows,
                        repair_cache;
                        max_steps_per_restart = hill_max_steps,
                    )
                end
            end
        end
    end
    compressed_distribution_path = gzip_file!(distribution_path)

    println("Summary written to: ", summary_path)
    println("Distribution written to: ", compressed_distribution_path)
    println("Retained flows written to: ", retained_path)
end

main()
