#!/usr/bin/env julia

import Pkg

const STUDENT_ROOT = normpath(abspath(joinpath(@__DIR__, "..")))
const DEFAULT_JULIQAOA_PROJECT = abspath(joinpath(@__DIR__, "..", "..", "..", "..", "JuliQAOA.jl"))
const JULIQAOA_PROJECT = get(ENV, "JULIQAOA_PROJECT", DEFAULT_JULIQAOA_PROJECT)

Pkg.activate(JULIQAOA_PROJECT)
Pkg.instantiate()

using JuliQAOA
using Printf
using Random
using Statistics

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

bits_signature(bits) = join(string.(bits))

function mask_from_bits(bits::AbstractString)
    mask = 0
    for (index, char) in enumerate(bits)
        char == '1' && (mask |= 1 << (index - 1))
    end
    return mask
end

function enumerate_surrogate_energies(L, Q; scale, offset)
    n = length(L)
    energies = Vector{Float64}(undef, 1 << n)
    bitstrings = Vector{String}(undef, 1 << n)
    for mask in 0:((1 << n) - 1)
        bits = [(mask >> (i - 1)) & 1 for i in 1:n]
        energy = offset
        for i in 1:n
            energy += L[i] * bits[i]
        end
        for i in 1:n, j in 1:n
            energy += Q[i, j] * bits[i] * bits[j]
        end
        energies[mask + 1] = scale * energy
        bitstrings[mask + 1] = bits_signature(bits)
    end
    return energies, bitstrings
end

function csv_join(values)
    return join(map(string, values), ',')
end

function angles_string(values)
    return join((@sprintf("%.17g", value) for value in values), ';')
end

function qiskit_angle_string(angles, energy_scale)
    p = length(angles) ÷ 2
    values = vcat(angles[1:p], angles[(p + 1):end] ./ energy_scale)
    return angles_string(values)
end

function target_measure(name, normalized_energy, masks_by_rank)
    if name == "energy"
        return (measure = normalized_energy, maximize = false)
    elseif name == "top50"
        measure = zeros(Float64, length(normalized_energy))
        for mask in values(masks_by_rank)
            measure[mask + 1] = 1.0
        end
        return (measure = measure, maximize = true)
    elseif name == "top10"
        measure = zeros(Float64, length(normalized_energy))
        for rank in 1:10
            mask = get(masks_by_rank, rank, nothing)
            isnothing(mask) || (measure[mask + 1] = 1.0)
        end
        return (measure = measure, maximize = true)
    elseif name == "global"
        measure = zeros(Float64, length(normalized_energy))
        measure[masks_by_rank[1] + 1] = 1.0
        return (measure = measure, maximize = true)
    else
        error("Unknown target $(name). Use one of: energy, top50, top10, global.")
    end
end

function probability_summary(probabilities, surrogate_energy, top_flows, bitstrings)
    total_top50 = 0.0
    total_top10 = 0.0
    total_global = 0.0
    best_rank = 0
    best_rank_probability = 0.0
    best_rank_energy = NaN
    best_rank_flow = ""

    for (flow, hit) in top_flows
        mask = mask_from_bits(flow)
        probability = probabilities[mask + 1]
        total_top50 += probability
        hit.rank <= 10 && (total_top10 += probability)
        hit.rank == 1 && (total_global += probability)
        if hit.rank == 1 || probability > best_rank_probability
            best_rank = hit.rank
            best_rank_probability = probability
            best_rank_energy = hit.exact_repaired_qubo_energy
            best_rank_flow = flow
        end
    end

    max_probability_index = argmax(probabilities)
    min_surrogate_index = argmin(surrogate_energy)
    return (
        top50_probability = total_top50,
        top10_probability = total_top10,
        global_probability = total_global,
        best_top50_rank_by_probability = best_rank,
        best_top50_probability = best_rank_probability,
        best_top50_exact_repaired_qubo_energy = best_rank_energy,
        best_top50_flow_bits = best_rank_flow,
        most_probable_flow_bits = bitstrings[max_probability_index],
        most_probable_probability = probabilities[max_probability_index],
        most_probable_surrogate_energy = surrogate_energy[max_probability_index],
        surrogate_min_flow_bits = bitstrings[min_surrogate_index],
        surrogate_min_probability = probabilities[min_surrogate_index],
        surrogate_min_energy = surrogate_energy[min_surrogate_index],
    )
end

function main()
    reduced_dir = joinpath(STUDENT_ROOT, "ds_mfg_reduced_flow_objective")
    output_dir = joinpath(STUDENT_ROOT, get(ENV, "DSMFG_JULIQAOA_OUTPUT_DIR", "ds_mfg_qaoa_juliqaoa_angle_search"))
    mkpath(output_dir)

    p_values = parse_int_list(get(ENV, "DSMFG_JULIQAOA_P_VALUES", "1:5"))
    max_p = maximum(p_values)
    niter = parse(Int, get(ENV, "DSMFG_JULIQAOA_NITER", "10"))
    seed = parse(Int, get(ENV, "DSMFG_JULIQAOA_SEED", "91001"))
    targets = split(get(ENV, "DSMFG_JULIQAOA_TARGETS", "energy,top50,top10,global"), ',')

    Random.seed!(seed)

    scalars = read_reduced_scalars(joinpath(reduced_dir, "reduced_scalars.csv"))
    L = read_reduced_vector(joinpath(reduced_dir, "reduced_L_vector.csv"))
    Q = read_reduced_matrix(joinpath(reduced_dir, "reduced_Q_matrix.csv"))
    top_flows = read_exact_top_flows(joinpath(reduced_dir, "reduced_exact_top_flows.csv"))
    masks_by_rank = Dict(hit.rank => mask_from_bits(flow) for (flow, hit) in top_flows)
    surrogate_energy, bitstrings = enumerate_surrogate_energies(L, Q; scale = scalars.scale, offset = scalars.offset)

    energy_shift = mean(surrogate_energy)
    energy_scale = std(surrogate_energy)
    normalized_energy = (surrogate_energy .- energy_shift) ./ energy_scale
    mixer = mixer_x(scalars.n)

    summary_path = joinpath(output_dir, "juliqaoa_angle_summary.csv")
    open(summary_path, "w") do io
        println(
            io,
            csv_join([
                "target", "p", "seed", "basinhopping_niter", "energy_shift", "energy_scale",
                "juliqaoa_expected_surrogate_energy", "juliqaoa_expected_normalized_energy",
                "top50_probability", "top10_probability", "global_probability",
                "best_top50_rank_by_probability", "best_top50_probability",
                "best_top50_exact_repaired_qubo_energy", "best_top50_flow_bits",
                "most_probable_flow_bits", "most_probable_probability", "most_probable_surrogate_energy",
                "surrogate_min_flow_bits", "surrogate_min_probability", "surrogate_min_energy",
                "normalized_angles_beta_then_gamma", "qiskit_angles_beta_then_gamma",
            ]),
        )

        for raw_target in targets
            target = strip(raw_target)
            isempty(target) && continue
            target_config = target_measure(target, normalized_energy, masks_by_rank)
            @info "Finding JuliQAOA angles" target max_p niter seed
            angle_sets, expected_values = find_angles_bh(
                max_p,
                mixer,
                normalized_energy,
                target_config.measure;
                max = target_config.maximize,
                niter = niter,
                verbose = true,
            )

            for p in p_values
                angles = angle_sets[p]
                probabilities = JuliQAOA.probabilities(angles, mixer, normalized_energy)
                expected_surrogate_energy = sum(probabilities .* surrogate_energy)
                expected_normalized_energy = sum(probabilities .* normalized_energy)
                summary = probability_summary(probabilities, surrogate_energy, top_flows, bitstrings)
                println(
                    io,
                    csv_join([
                        target, p, seed, niter, energy_shift, energy_scale,
                        expected_surrogate_energy, expected_normalized_energy,
                        summary.top50_probability, summary.top10_probability, summary.global_probability,
                        summary.best_top50_rank_by_probability, summary.best_top50_probability,
                        summary.best_top50_exact_repaired_qubo_energy, summary.best_top50_flow_bits,
                        summary.most_probable_flow_bits, summary.most_probable_probability,
                        summary.most_probable_surrogate_energy, summary.surrogate_min_flow_bits,
                        summary.surrogate_min_probability, summary.surrogate_min_energy,
                        angles_string(angles), qiskit_angle_string(angles, energy_scale),
                    ]),
                )
                flush(io)
            end
        end
    end

    println("JuliQAOA angle summary written to: ", summary_path)
end

main()
