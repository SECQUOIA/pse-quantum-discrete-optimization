#!/usr/bin/env julia
"""Reproducible, CPU-only robustness audit; never submits hardware jobs.

Run with Julia 1.10, one Julia/BLAS thread. Each classical run starts with an
empty repair cache; enumeration and surrogate information are not search inputs.
"""
module RevisionSensitivity

using LinearAlgebra, Random, Statistics, Printf, Dates, TOML, SHA
module Baselines
include(joinpath(@__DIR__, "run_classical_baselines.jl"))
end
const B = Baselines
const ROOT = B.STUDY_ROOT
const OUTPUT = joinpath(ROOT, "ds_mfg_revision_audit")
# Original IP capital costs, in the f00..f18 order used by
# scripts/enumerate_ip_provenance.py. Tests verify repair and scaling against
# brute-force enumeration of the original QUBO with these fixed capital costs.
const COST = [0., .0289, .0399, .1225, 0., 0., 4.959, 0., 6.193,
              9.717, 3.352, 6.163, 0., .739, .84, 1.271, 6.99, 1.333, 0.]
const SEEDS = 97001:97005
const FIT_SEEDS = 97101:97103
const BUDGET = 262144
const N = 19
const REFERENCE = "1001110100111100011"

function write_rows(path, rows)
    isempty(rows) && error("No rows for $path")
    names = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(names), ','))
        for row in rows
            # Values in these numeric/provenance tables contain no CSV delimiters.
            values = [string(getproperty(row, name)) for name in names]
            any(v -> occursin(r"[,\n\r]", v), values) && error("CSV delimiter in value")
            println(io, join(values, ','))
        end
    end
end

function load_problem()
    qubo = B.read_qubo_data(ROOT)
    components = B.read_aux_components(joinpath(ROOT,
        "ds_mfg_reduced_flow_objective", "auxiliary_components.csv"), qubo)
    top = B.read_exact_top_flows(joinpath(ROOT,
        "ds_mfg_reduced_flow_objective", "reduced_exact_top_flows.csv"))
    B.validate_repair(qubo, components, top)
    return qubo, components, top
end

function enumerate_objective(qubo, components)
    energy = Vector{Float64}(undef, 1 << N)
    cost = similar(energy)
    for mask in 0:((1 << N)-1)
        y = B.flow_from_index(mask, N) # f00 is the most significant bit.
        energy[mask+1] = B.repair_flow(y, qubo, components).exact_repaired_qubo_energy
        cost[mask+1] = dot(COST, y)
    end
    return energy, cost
end

function penalty_audit(energy, cost, top)
    feasible = falses(length(energy))
    for (bits, row) in top
        row.match in ("gurobi_pool", "global_optimum") || continue
        feasible[parse(Int, bits; base=2)+1] = true
    end
    count(feasible) == 36 || error("Expected the certified feasible set")
    @assert maximum(abs.(energy[feasible] .- cost[feasible])) < 1e-8
    penalty = energy .- cost
    minimum(penalty) >= -1e-8 || error("Negative penalty residual")
    rows = NamedTuple[]
    for factor in (.01, .05, .1, .25, .5, 1., 2., 4.)
        # For factor > 0, min_z [c(y) + factor P(y,z)] = c(y) + factor min_z P.
        score = cost .+ factor .* penalty
        order = sortperm(score)
        best = first(order)
        push!(rows, (penalty_multiplier=factor, optimum=score[best],
            optimum_flow=bitstring(UInt32(best-1))[end-N+1:end],
            optimum_feasible=feasible[best],
            feasible_in_top50=count(feasible[order[1:50]]),
            min_infeasible_minus_max_feasible=minimum(score[.!feasible])-maximum(score[feasible])))
    end
    return rows
end

function walsh!(a)
    h = 1
    while h < length(a)
        for start in 1:2h:length(a), j in start:(start+h-1)
            x, y = a[j], a[j+h]
            a[j], a[j+h] = x+y, x-y
        end
        h *= 2
    end
    return a
end

feature_masks(n) = vcat(0, [1 << i for i in 0:n-1],
    [(1 << i) | (1 << j) for i in 0:n-1 for j in i+1:n-1])

function design_matrix(indices, masks)
    X = Matrix{Float64}(undef, length(indices), length(masks))
    for (j, mask) in enumerate(masks), (i, index) in enumerate(indices)
        X[i,j] = isodd(count_ones((index-1) & mask)) ? -1. : 1.
    end
    return X
end

function predict(coefficients, masks, n)
    spectrum = zeros(1 << n)
    spectrum[masks .+ 1] .= coefficients
    return walsh!(spectrum)
end

function fit_full(energy, n)
    masks = feature_masks(n)
    spectrum = walsh!(copy(energy)) ./ length(energy)
    # Walsh features are orthonormal under the uniform Boolean-cube measure.
    # Keeping degree <= 2 is exactly the unweighted quadratic least-squares fit.
    return spectrum[masks .+ 1], masks
end

function fit_metrics(prediction, energy, order)
    best = argmin(prediction)
    reference = parse(Int, REFERENCE; base=2)+1
    return (global_r2=1-sum(abs2, prediction-energy)/sum(abs2, energy.-mean(energy)),
        optimum_repaired_energy=energy[best],
        optimum_flow=bitstring(UInt32(best-1))[end-N+1:end],
        reference_surrogate_rank=1+count(<(prediction[reference]-1e-8), prediction),
        top50_max_absolute_error=maximum(abs.(prediction[order[1:50]]-energy[order[1:50]])))
end

function surrogate_audit(energy)
    order = sortperm(energy)
    rows = NamedTuple[]
    start = time_ns()
    coefficients, masks = fit_full(energy, N)
    fit_seconds = (time_ns()-start)/1e9
    prediction = predict(coefficients, masks, N)
    push!(rows, merge((construction="full_uniform_walsh", training_flows=length(energy),
        seed=0, fit_seconds=fit_seconds), fit_metrics(prediction, energy, order)))
    for size in (1024, 8192, 65536), seed in FIT_SEEDS
        rng = MersenneTwister(seed)
        start = time_ns()
        indices = randperm(rng, length(energy))[1:size]
        X = design_matrix(indices, masks)
        coefficients = X \ energy[indices]
        elapsed = (time_ns()-start)/1e9
        prediction = predict(coefficients, masks, N)
        push!(rows, merge((construction="uniform_subset_least_squares", training_flows=size,
            seed=seed, fit_seconds=elapsed), fit_metrics(prediction, energy, order)))
    end
    return rows
end

function annealing_distribution(rng, budget, qubo, components, cache;
        restart_evaluations=4096, temperature_start=100., temperature_end=.01)
    restart_evaluations >= 3 || error("Restart budget must be >= 3")
    temperature_start >= temperature_end > 0 || error("Invalid temperature schedule")
    distribution = Dict{String,Int}()
    evaluations, restarts = 0, 0
    while evaluations < budget
        restarts += 1
        flow = B.random_flow(rng, N)
        current = B.score_flow!(distribution, cache, flow, qubo, components)
        evaluations += 1
        for step in 1:restart_evaluations-1
            evaluations >= budget && break
            candidate = copy(flow)
            bit = rand(rng, 1:N)
            candidate[bit] = 1-candidate[bit]
            record = B.score_flow!(distribution, cache, candidate, qubo, components)
            evaluations += 1
            temperature = temperature_start * (temperature_end/temperature_start)^(
                (step-1)/(restart_evaluations-2))
            delta = record.exact_repaired_qubo_energy-current.exact_repaired_qubo_energy
            if delta <= 0 || rand(rng) < exp(-delta/temperature)
                flow, current = candidate, record
            end
        end
    end
    return distribution, evaluations, restarts
end

function baseline_run(method, seed, budget, qubo, components)
    cache = Dict{String,NamedTuple}()
    rng = MersenneTwister(seed)
    start = time_ns()
    distribution, evaluations, restarts = if method == "uniform"
        B.uniform_random_distribution(rng, budget, N, qubo, components, cache)
    elseif method == "hill_climb"
        B.hill_climb_distribution(rng, budget, N, qubo, components, cache; max_steps_per_restart=100)
    elseif method == "simulated_annealing"
        annealing_distribution(rng, budget, qubo, components, cache)
    else
        error("Unknown method: $method")
    end
    elapsed = (time_ns()-start)/1e9
    @assert evaluations == sum(values(distribution)) == budget
    best = minimum(r.exact_repaired_qubo_energy for r in values(cache))
    return (algorithm=method, seed=seed, evaluations=evaluations, restarts=restarts,
        search_seconds=elapsed, unique_evaluations=length(cache),
        optimum_evaluations=get(distribution, REFERENCE, 0),
        found_optimum=haskey(distribution, REFERENCE), best_repaired_energy=best)
end

function main()
    BLAS.set_num_threads(1)
    Threads.nthreads() == 1 || error("Run with one Julia thread")
    output = isempty(ARGS) ? OUTPUT : abspath(ARGS[1])
    mkpath(output)
    start = time_ns()
    qubo, components, top = load_problem()
    setup_seconds = (time_ns()-start)/1e9
    # Compile every search path with a throwaway cache before timing any run.
    methods = ("uniform", "hill_climb", "simulated_annealing")
    for method in methods
        baseline_run(method, 1, 256, qubo, components)
    end
    # Search before enumeration, with a new empty cache for every method/seed.
    baselines = NamedTuple[]
    for seed in SEEDS, method in methods
        row = baseline_run(method, seed, BUDGET, qubo, components)
        push!(baselines, row)
        println(row)
        write_rows(joinpath(output, "classical_runs.csv"), baselines)
    end
    start = time_ns()
    energy, cost = enumerate_objective(qubo, components)
    enumeration_seconds = (time_ns()-start)/1e9
    @assert isapprox(minimum(energy), 11.7095; atol=1e-8)
    write_rows(joinpath(output, "penalty_sensitivity.csv"), penalty_audit(energy, cost, top))
    # Warm the transform/QR code on a small, independent toy problem.
    fit_full(collect(1.:8.), 3)
    design_matrix(1:8, feature_masks(3)) \ collect(1.:8.)
    write_rows(joinpath(output, "surrogate_sensitivity.csv"), surrogate_audit(energy))
    write_rows(joinpath(output, "preprocessing_times.csv"), [
        (stage="load_and_validate_original_qubo", seconds=setup_seconds, scope="includes_first_call_compilation"),
        (stage="enumerate_and_repair_all_524288_flows", seconds=enumeration_seconds, scope="warm_repair_kernel_includes_flow_and_cost_construction")])
    metadata = Dict("julia_version"=>string(VERSION), "julia_threads"=>Threads.nthreads(),
        "blas_threads"=>BLAS.get_num_threads(), "cpu"=>Sys.cpu_info()[1].model,
        "created_utc"=>string(now(UTC)), "baseline_seeds"=>collect(SEEDS),
        "fit_seeds"=>collect(FIT_SEEDS), "evaluation_budget"=>BUDGET,
        "annealing_restart_evaluations"=>4096, "annealing_temperature_start"=>100.,
        "annealing_temperature_end"=>.01, "hill_climb_max_steps"=>100,
        "timing_policy"=>"warm search with per-run empty cache; no enumerated lookup; timers exclude process startup and artifact writing",
        "script_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "original_qubo_sha256"=>bytes2hex(sha256(read(joinpath(ROOT,B.ZIP_NAME)))))
    open(joinpath(output, "run_metadata.toml"), "w") do io
        TOML.print(io, metadata)
    end
    println("Revision audit written to $output")
end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    RevisionSensitivity.main()
end
