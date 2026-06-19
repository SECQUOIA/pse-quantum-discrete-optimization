#!/usr/bin/env julia

import Pkg

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(joinpath(@__DIR__, ".."))
    Pkg.instantiate()
end

using Printf
using PythonCall

const STUDY_ROOT = normpath(abspath(joinpath(@__DIR__, "..")))
const ZIP_NAME = "Fw_ DS mfg case qubo information.zip"
const REDUCED_DIR = joinpath(STUDY_ROOT, "ds_mfg_reduced_flow_objective")
const DEFAULT_OUTPUT_DIR = "ds_mfg_direct_full_qubo_audit"
const DIRECT_QAOA_ANGLE_SOURCE = "full-QUBO p=2 final-sampling separation"
const DIRECT_QAOA_ANGLES = [0.555, 0.293, -0.488, -0.898]
const DIRECT_QAOA_P = 2
const DIRECT_QAOA_SOURCE_SUMMARY = joinpath("ds_mfg_final_sampling_sweep_v2", "final_sampling_summary.csv")
const DIRECT_QAOA_SOURCE_DISTRIBUTION = joinpath(
    "ds_mfg_final_sampling_sweep_v2",
    "qaoa_p2_optimizer_reads128_final_reads512_iter25_distribution.csv",
)
const DIRECT_QAOA_PARAMETER_ARTIFACT = "direct_full_qubo_qaoa_hardware_parameters.json"
const DIRECT_FULL_QUBO_ALGORITHM = "QAOA_direct_full_qubo_resource_audit"

include(joinpath(@__DIR__, "hit_rate_stats.jl"))

struct QuboData
    n::Int
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
    if value === nothing || ismissing(value)
        return ""
    end
    text = string(value)
    if occursin(',', text) || occursin('"', text) || occursin('\n', text) || occursin('\r', text)
        return "\"" * replace(text, "\"" => "\"\"") * "\""
    end
    return text
end

csv_row(values) = join(csv_value.(values), ',')

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

function parse_bool_env(name::AbstractString)
    raw = lowercase(strip(get(ENV, name, "false")))
    raw in ("", "false", "0", "no") && return false
    raw in ("true", "1", "yes") && return true
    error("$(name) must be true or false, got $(raw)")
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

function nonempty_lines(text::AbstractString)
    return [strip(line) for line in split(text, '\n') if !isempty(strip(line))]
end

function parse_scalar_csv(text::AbstractString)
    rows = nonempty_lines(text)
    length(rows) >= 2 || error("Expected header and value row in scalars.csv")
    values = split(rows[2], ','; keepempty = true)
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
    matrix = Matrix{Float64}(undef, length(rows), length(split(rows[1], ','; keepempty = true)))
    for (i, row) in enumerate(rows)
        values = split(row, ','; keepempty = true)
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
    return QuboData(scalars.n, scalars.scale, scalars.offset, linear, quadratic)
end

function read_csv_rows(path::AbstractString)
    lines = filter(line -> !isempty(strip(line)), readlines(path))
    length(lines) >= 1 || error("Expected CSV header in $(path)")
    header = split(lines[1], ','; keepempty = true)
    rows = Dict{String,String}[]
    for (line_number, line) in enumerate(lines[2:end])
        values = split(line, ','; keepempty = true)
        length(values) == length(header) ||
            error("CSV row $(line_number + 1) in $(path) has $(length(values)) fields; expected $(length(header))")
        push!(rows, Dict(header[i] => values[i] for i in eachindex(header)))
    end
    return rows
end

function read_aux_components(path::AbstractString, qubo::QuboData; n_flow::Int = 19)
    rows = read_csv_rows(path)
    components = AuxComponent[]
    for row in rows
        aux_indices = [parse(Int, text) for text in split(row["aux_indices"], ';') if !isempty(strip(text))]
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
    top = Dict{String,NamedTuple}()
    for row in read_csv_rows(path)
        top[row["flow_bits"]] = (
            rank = parse(Int, row["rank"]),
            exact_repaired_qubo_energy = parse(Float64, row["exact_repaired_qubo_energy"]),
            repaired_aux_bits = row["repaired_aux_bits"],
            repaired_full_bits = row["repaired_full_bits"],
            match = row["match"],
            matched_ip_obj = isempty(row["matched_ip_obj"]) ? missing : parse(Float64, row["matched_ip_obj"]),
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

function bits_vector(bits::AbstractString)
    return [char == '1' ? 1 : 0 for char in bits]
end

function rank_info(flow_bits::String, top_flows::Dict{String,NamedTuple})
    hit = get(top_flows, flow_bits, nothing)
    isnothing(hit) && return (
        rank = missing,
        match = "not_top50",
        matched_ip_obj = missing,
        repaired_full_bits = "",
        repaired_aux_bits = "",
    )
    return (
        rank = hit.rank,
        match = hit.match,
        matched_ip_obj = hit.matched_ip_obj,
        repaired_full_bits = hit.repaired_full_bits,
        repaired_aux_bits = hit.repaired_aux_bits,
    )
end

function parse_int_field(row::Dict{String,String}, keys::AbstractVector{<:AbstractString})
    for key in keys
        haskey(row, key) && !isempty(row[key]) && return parse(Int, row[key])
    end
    return missing
end

function parse_float_field(row::Dict{String,String}, keys::AbstractVector{<:AbstractString})
    for key in keys
        haskey(row, key) && !isempty(row[key]) && return parse(Float64, row[key])
    end
    return missing
end

function string_field(row::Dict{String,String}, keys::AbstractVector{<:AbstractString}; default = "")
    for key in keys
        haskey(row, key) && return row[key]
    end
    return default
end

function source_specs()
    specs = Vector{Dict{String,Any}}()

    for row in read_csv_rows(joinpath(STUDY_ROOT, "ds_mfg_fixed_angle_sweep", "fixed_angle_qaoa_summary.csv"))
        push!(specs, Dict{String,Any}(
            "algorithm" => "QAOA_full_qubo_fixed_angle",
            "source_artifact" => row["distribution_csv"],
            "seed" => "",
            "p" => row["p"],
            "optimizer_reads" => row["number_of_reads"],
            "final_reads" => row["number_of_reads"],
            "maximum_iterations" => row["maximum_iterations"],
            "solve_time_sec" => row["solve_time_sec"],
            "settings" => "angle_source=$(row["angle_source"]);gamma_sign=$(row["gamma_sign"]);guarantee_3regular_maxcut=$(row["guarantee_3regular_maxcut"])",
        ))
    end

    for row in read_csv_rows(joinpath(STUDY_ROOT, "ds_mfg_iteration_sweep", "sweep_summary.csv"))
        push!(specs, Dict{String,Any}(
            "algorithm" => "QAOA_full_qubo_iteration_sweep",
            "source_artifact" => row["distribution_csv"],
            "seed" => row["seed"],
            "p" => row["number_of_layers"],
            "optimizer_reads" => row["number_of_reads"],
            "final_reads" => row["number_of_reads"],
            "maximum_iterations" => row["maximum_iterations"],
            "solve_time_sec" => row["solve_time_sec"],
            "settings" => "trial=$(row["trial"]);initial_parameters=$(row["initial_parameters"])",
        ))
    end

    for row in read_csv_rows(joinpath(STUDY_ROOT, "ds_mfg_p2_perturbation_sweep_v3", "p2_perturbation_summary.csv"))
        push!(specs, Dict{String,Any}(
            "algorithm" => "QAOA_full_qubo_p2_perturbation",
            "source_artifact" => row["distribution_csv"],
            "seed" => row["seed"],
            "p" => row["p"],
            "optimizer_reads" => row["number_of_reads"],
            "final_reads" => row["number_of_reads"],
            "maximum_iterations" => row["maximum_iterations"],
            "solve_time_sec" => row["solve_time_sec"],
            "settings" => "sigma=$(row["sigma"]);trial=$(row["trial"]);initial_parameters=$(row["initial_parameters"])",
        ))
    end

    for row in read_csv_rows(joinpath(STUDY_ROOT, "ds_mfg_final_sampling_sweep_v2", "final_sampling_summary.csv"))
        push!(specs, Dict{String,Any}(
            "algorithm" => "QAOA_full_qubo_final_sampling",
            "source_artifact" => row["distribution_csv"],
            "seed" => "",
            "p" => row["p"],
            "optimizer_reads" => row["optimizer_reads"],
            "final_reads" => row["final_reads"],
            "maximum_iterations" => row["maximum_iterations"],
            "solve_time_sec" => row["solve_time_sec"],
            "settings" => "separated_optimizer_and_final_sampling=true",
        ))
    end

    for row in read_csv_rows(joinpath(STUDY_ROOT, "ds_mfg_direct_full_qubo_qaoa_highread", "direct_full_qubo_qaoa_highread_summary.csv"))
        push!(specs, Dict{String,Any}(
            "algorithm" => row["algorithm"],
            "source_artifact" => row["distribution_csv"],
            "seed" => row["seed"],
            "p" => row["p"],
            "optimizer_reads" => row["optimizer_reads"],
            "final_reads" => row["final_reads"],
            "maximum_iterations" => row["maximum_iterations"],
            "solve_time_sec" => row["solve_time_sec"],
            "settings" => "parameter_source=$(row["parameter_source"]);optimized_parameter_artifact=$(row["optimized_parameter_artifact"]);parameter_order=$(row["parameter_order"])",
        ))
    end

    for row in read_csv_rows(joinpath(STUDY_ROOT, "ds_mfg_vqe_final_sampling_sweep_v2", "vqe_final_sampling_summary.csv"))
        push!(specs, Dict{String,Any}(
            "algorithm" => "VQE_full_qubo_final_sampling",
            "source_artifact" => row["distribution_csv"],
            "seed" => row["seed"],
            "p" => "",
            "optimizer_reads" => row["optimizer_reads"],
            "final_reads" => row["final_reads"],
            "maximum_iterations" => row["maximum_iterations"],
            "solve_time_sec" => row["solve_time_sec"],
            "settings" => "ansatz=EfficientSU2;initial_parameter_count=$(row["initial_parameter_count"])",
        ))
    end

    return specs
end

function score_distribution(
    spec::Dict{String,Any},
    qubo::QuboData,
    components::Vector{AuxComponent},
    top_flows::Dict{String,NamedTuple},
)
    distribution_path = joinpath(STUDY_ROOT, spec["source_artifact"])
    rows = read_csv_rows(distribution_path)
    total_reads = 0
    top50_hits = 0
    top10_hits = 0
    global_hits = 0
    projected_pool_hits = 0
    encoded_pool_hits = 0
    encoded_global_hits = 0
    best_raw_energy = Inf
    best_raw_flow_bits = ""
    best_raw_projection_match = ""
    best_repaired = nothing
    best_repaired_flow_bits = ""
    repair_cache = Dict{String,NamedTuple}()

    for row in rows
        reads = parse(Int, row["reads"])
        total_reads += reads
        flow_bits = row["flow_bits"]
        full_bits = row["full_bits"]
        info = rank_info(flow_bits, top_flows)
        is_pool = is_gurobi_pool_feasible_match(info.match)
        encoded = is_pool && full_bits == info.repaired_full_bits

        if !ismissing(info.rank)
            top50_hits += reads
            info.rank <= 10 && (top10_hits += reads)
            info.rank == 1 && (global_hits += reads)
        end
        is_pool && (projected_pool_hits += reads)
        encoded && (encoded_pool_hits += reads)
        encoded && !ismissing(info.rank) && info.rank == 1 && (encoded_global_hits += reads)

        raw_energy = parse_float_field(row, ["qubo_energy", "surrogate_qubo_energy", "objective_value"])
        if !ismissing(raw_energy) && raw_energy < best_raw_energy
            best_raw_energy = raw_energy
            best_raw_flow_bits = flow_bits
            best_raw_projection_match = string_field(row, ["projection_match", "match"]; default = info.match)
        end

        repaired = get!(repair_cache, flow_bits) do
            repair_flow(bits_vector(flow_bits), qubo, components)
        end
        if isnothing(best_repaired) ||
            repaired.exact_repaired_qubo_energy < best_repaired.exact_repaired_qubo_energy
            best_repaired = repaired
            best_repaired_flow_bits = flow_bits
        end
    end

    total_reads == sum(parse(Int, row["reads"]) for row in rows; init = 0) ||
        error("Read-count accounting failed for $(spec["source_artifact"])")
    elapsed_text = @sprintf("%.6f", parse(Float64, spec["solve_time_sec"]))
    elapsed_sec = parse(Float64, elapsed_text)
    best_info = rank_info(best_repaired_flow_bits, top_flows)
    hit_stats = hit_rate_stat_values(
        total_reads,
        top50_hits,
        top10_hits,
        global_hits,
        projected_pool_hits,
        elapsed_sec,
    )

    return [
        spec["algorithm"],
        spec["source_artifact"],
        "direct_original_36_variable_qubo",
        "cached_local_aer_emulation",
        spec["settings"],
        spec["seed"],
        spec["p"],
        spec["optimizer_reads"],
        spec["final_reads"],
        spec["maximum_iterations"],
        elapsed_text,
        length(rows),
        total_reads,
        projected_pool_hits,
        encoded_pool_hits,
        encoded_global_hits,
        top50_hits,
        top10_hits,
        global_hits,
        projected_pool_hits,
        hit_stats...,
        isfinite(best_raw_energy) ? best_raw_energy : "",
        best_raw_flow_bits,
        best_raw_projection_match,
        best_info.rank,
        isnothing(best_repaired) ? "" : best_repaired.exact_repaired_qubo_energy,
        best_info.match,
        best_repaired_flow_bits,
        isnothing(best_repaired) ? "" : best_repaired.repaired_aux_bits,
        isnothing(best_repaired) ? "" : best_repaired.repaired_full_bits,
    ]
end

function direct_summary_header()
    return [
        "algorithm", "source_artifact", "problem", "mode", "settings", "seed", "p",
        "optimizer_reads", "final_reads", "maximum_iterations", "solve_time_sec",
        "unique_states", "total_reads", "projected_pool_hits", "encoded_pool_hits",
        "encoded_global_hits", "top50_hits", "top10_hits", "global_hits",
        "gurobi_pool_feasible_hits", hit_rate_stat_headers()..., "best_raw_qubo_energy",
        "best_raw_flow_bits", "best_raw_projection_match", "best_repaired_rank",
        "best_repaired_qubo_energy", "best_repaired_match", "best_repaired_flow_bits",
        "best_repaired_aux_bits", "best_repaired_full_bits",
    ]
end

function effective_qubo_terms(qubo::QuboData)
    linear = [qubo.scale * (qubo.linear[i] + qubo.quadratic[i, i]) for i in 1:qubo.n]
    pairs = Tuple{Int,Int,Float64}[]
    for i in 1:(qubo.n - 1), j in (i + 1):qubo.n
        coefficient = qubo.scale * (qubo.quadratic[i, j] + qubo.quadratic[j, i])
        abs(coefficient) > 1.0e-12 && push!(pairs, (i, j, coefficient))
    end
    return linear, pairs
end

function ising_terms_from_qubo(qubo::QuboData)
    linear, qubo_pairs = effective_qubo_terms(qubo)
    h = zeros(Float64, qubo.n)
    zz_pairs = Tuple{Int,Int,Float64}[]

    for i in eachindex(linear)
        h[i] -= linear[i] / 2.0
    end
    for (i, j, coefficient) in qubo_pairs
        h[i] -= coefficient / 4.0
        h[j] -= coefficient / 4.0
        push!(zz_pairs, (i, j, coefficient / 4.0))
    end
    return h, zz_pairs
end

function validate_direct_qaoa_parameters(p::Integer, angles::AbstractVector{<:Real})
    p > 0 || error("QAOA depth p must be positive, got $(p)")
    length(angles) == 2p || error("Expected $(2p) QAOA angles for p=$(p), got $(length(angles))")
    return Float64.(collect(angles))
end

function build_direct_full_qubo_qaoa_circuit(
    qubo::QuboData,
    p::Integer,
    angles::AbstractVector{<:Real},
)
    qiskit = pyimport("qiskit")
    parameter_values = validate_direct_qaoa_parameters(p, angles)
    betas = parameter_values[1:p]
    gammas = parameter_values[(p + 1):end]
    h, zz_pairs = ising_terms_from_qubo(qubo)

    circuit = qiskit.QuantumCircuit(qubo.n, qubo.n)
    for qubit in 0:(qubo.n - 1)
        circuit.h(qubit)
    end

    for layer in 1:p
        gamma = gammas[layer]
        for i in 1:qubo.n
            theta = 2.0 * gamma * h[i]
            abs(theta) > 1.0e-12 && circuit.rz(theta, i - 1)
        end
        for (i, j, coefficient) in zz_pairs
            theta = 2.0 * gamma * coefficient
            abs(theta) > 1.0e-12 && circuit.rzz(theta, i - 1, j - 1)
        end

        beta = betas[layer]
        for qubit in 0:(qubo.n - 1)
            circuit.rx(2.0 * beta, qubit)
        end
    end

    for qubit in 0:(qubo.n - 1)
        circuit.measure(qubit, qubit)
    end
    return circuit, zz_pairs
end

function build_direct_full_qubo_qaoa_circuit(qubo::QuboData)
    return build_direct_full_qubo_qaoa_circuit(qubo, DIRECT_QAOA_P, DIRECT_QAOA_ANGLES)
end

function direct_qaoa_parameter_artifact(qubo::QuboData)
    return Dict{String,Any}(
        "schema_version" => 1,
        "artifact_role" => "direct_full_qubo_qaoa_hardware_parameter_handoff",
        "algorithm" => "QAOA_direct_full_qubo_fixed_parameter_handoff",
        "problem_name" => "DS-MFG direct original 36-variable QUBO",
        "data_source" => ZIP_NAME,
        "n_qubits" => qubo.n,
        "n_flow_variables" => 19,
        "n_auxiliary_variables" => qubo.n - 19,
        "scale" => qubo.scale,
        "offset" => qubo.offset,
        "angle_source" => DIRECT_QAOA_ANGLE_SOURCE,
        "source_summary_artifact" => DIRECT_QAOA_SOURCE_SUMMARY,
        "source_distribution_artifact" => DIRECT_QAOA_SOURCE_DISTRIBUTION,
        "p" => DIRECT_QAOA_P,
        "parameter_order" => "beta_then_gamma",
        "qiskit_angles_beta_then_gamma" => DIRECT_QAOA_ANGLES,
        "hardware_handoff_policy" => "fixed-parameter sampling only; do not run an optimizer loop on hardware",
        "measurement_bit_order" => "Qiskit count keys must be reversed before scoring as x1..x36.",
        "producer_script" => relpath(@__FILE__, STUDY_ROOT),
    )
end

function py_dict_to_julia_int(py_dict)
    converted = Dict{String,Int}()
    for item in py_dict.items()
        converted[pyconvert(String, item[0])] = pyconvert(Int, item[1])
    end
    return converted
end

function circuit_resource_metadata(qubo::QuboData; transpile::Bool, fake_backend_class::String)
    parameter_artifact = direct_qaoa_parameter_artifact(qubo)
    circuit, zz_pairs = build_direct_full_qubo_qaoa_circuit(
        qubo,
        parameter_artifact["p"],
        parameter_artifact["qiskit_angles_beta_then_gamma"],
    )
    operation_counts = py_dict_to_julia_int(circuit.count_ops())
    rzz_count = get(operation_counts, "rzz", 0)
    cx_count = get(operation_counts, "cx", 0)
    logical_two_qubit_ops = rzz_count + cx_count
    amplitudes = BigInt(2)^qubo.n
    density_matrix_elements = BigInt(2)^(2 * qubo.n)
    metadata = Dict{String,Any}(
        "schema_version" => 1,
        "script" => relpath(@__FILE__, STUDY_ROOT),
        "problem" => Dict{String,Any}(
            "name" => "DS-MFG direct original 36-variable QUBO",
            "data_source" => ZIP_NAME,
            "n_qubits" => qubo.n,
            "n_flow_variables" => 19,
            "n_auxiliary_variables" => qubo.n - 19,
            "scale" => qubo.scale,
            "offset" => qubo.offset,
        ),
        "angle_source" => Dict{String,Any}(
            "description" => parameter_artifact["angle_source"],
            "p" => parameter_artifact["p"],
            "parameter_order" => parameter_artifact["parameter_order"],
            "qiskit_angles_beta_then_gamma" => parameter_artifact["qiskit_angles_beta_then_gamma"],
            "parameter_artifact" => DIRECT_QAOA_PARAMETER_ARTIFACT,
        ),
        "logical_circuit" => Dict{String,Any}(
            "num_qubits" => pyconvert(Int, circuit.num_qubits),
            "num_clbits" => pyconvert(Int, circuit.num_clbits),
            "depth" => pyconvert(Int, circuit.depth()),
            "operation_counts" => operation_counts,
            "nonzero_ising_zz_pairs" => length(zz_pairs),
            "logical_two_qubit_ops" => logical_two_qubit_ops,
            "measurement_bit_order" => "Qiskit count keys must be reversed before scoring as x1..x36.",
        ),
        "statevector_memory_estimate" => Dict{String,Any}(
            "amplitudes" => string(amplitudes),
            "complex64_bytes" => string(amplitudes * 8),
            "complex128_bytes" => string(amplitudes * 16),
        ),
        "classical_noisy_simulation_assessment" => Dict{String,Any}(
            "method_considered" => "qiskit_aer.AerSimulator.from_backend(FakeFez)",
            "dense_density_matrix_elements" => string(density_matrix_elements),
            "dense_density_matrix_complex64_bytes" => string(density_matrix_elements * 8),
            "dense_density_matrix_complex128_bytes" => string(density_matrix_elements * 16),
            "automatic_or_mps_note" => "Aer automatic or MPS-style noisy simulation may avoid dense density-matrix storage for some circuits but the cost is entanglement-dependent and was not run at a useful shot budget for this 36-qubit circuit.",
            "feasibility" => "Exact dense noisy simulation is not feasible on ordinary local hardware; approximate or tensor-network simulation would need a separate bounded feasibility study.",
        ),
        "fake_backend_class" => fake_backend_class,
        "transpile_requested" => transpile,
        "transpile" => Dict{String,Any}("status" => transpile ? "requested" : "not_requested"),
        "simulation_status" => "not_run",
        "blocker" => "No direct 36-qubit noisy samples are cached. Dense noisy density-matrix simulation would require 2^72 complex entries. Aer automatic or MPS simulation is entanglement-dependent and was not run at a useful shot budget.",
    )

    if transpile
        qiskit = pyimport("qiskit")
        fake_provider = pyimport("qiskit_ibm_runtime.fake_provider")
        pyhasattr(fake_provider, fake_backend_class) ||
            error("qiskit_ibm_runtime.fake_provider has no $(fake_backend_class)")
        backend = pygetattr(fake_provider, fake_backend_class)()
        transpiled = qiskit.transpile(circuit; backend = backend, seed_transpiler = 92001, optimization_level = 3)
        transpiled_counts = py_dict_to_julia_int(transpiled.count_ops())
        metadata["transpile"] = Dict{String,Any}(
            "status" => "DONE",
            "fake_backend_class" => fake_backend_class,
            "seed_transpiler" => 92001,
            "optimization_level" => 3,
            "elapsed_sec" => "not_recorded_for_reproducible_artifacts",
            "depth" => pyconvert(Int, transpiled.depth()),
            "operation_counts" => transpiled_counts,
            "cz_count" => get(transpiled_counts, "cz", 0),
            "cx_count" => get(transpiled_counts, "cx", 0),
            "ecr_count" => get(transpiled_counts, "ecr", 0),
        )
    end

    return metadata
end

function resource_summary_header()
    return [
        "problem", "algorithm", "mode", "fake_backend", "angle_source", "p",
        "logical_qubits", "logical_clbits", "logical_depth", "logical_rzz_count",
        "logical_cx_count", "logical_measure_count", "logical_two_qubit_ops",
        "statevector_amplitudes", "statevector_complex64_bytes", "statevector_complex128_bytes",
        "density_matrix_elements", "density_matrix_complex64_bytes", "density_matrix_complex128_bytes",
        "noisy_method_considered", "classical_noisy_feasibility",
        "transpile_requested", "transpile_status", "transpiled_depth", "transpiled_cz_count",
        "transpiled_cx_count", "transpiled_ecr_count", "transpile_elapsed_sec",
        "simulation_status", "blocker",
    ]
end

function write_resource_summary(path::AbstractString, metadata)
    logical = metadata["logical_circuit"]
    operations = logical["operation_counts"]
    memory = metadata["statevector_memory_estimate"]
    noisy = metadata["classical_noisy_simulation_assessment"]
    transpile = metadata["transpile"]
    open(path, "w") do io
        println(io, csv_row(resource_summary_header()))
        println(io, csv_row([
            "direct_original_36_variable_qubo",
            DIRECT_FULL_QUBO_ALGORITHM,
            "resource_audit",
            metadata["fake_backend_class"],
            DIRECT_QAOA_ANGLE_SOURCE,
            DIRECT_QAOA_P,
            logical["num_qubits"],
            logical["num_clbits"],
            logical["depth"],
            get(operations, "rzz", 0),
            get(operations, "cx", 0),
            get(operations, "measure", 0),
            logical["logical_two_qubit_ops"],
            memory["amplitudes"],
            memory["complex64_bytes"],
            memory["complex128_bytes"],
            noisy["dense_density_matrix_elements"],
            noisy["dense_density_matrix_complex64_bytes"],
            noisy["dense_density_matrix_complex128_bytes"],
            noisy["method_considered"],
            noisy["feasibility"],
            metadata["transpile_requested"],
            get(transpile, "status", ""),
            get(transpile, "depth", ""),
            get(transpile, "cz_count", ""),
            get(transpile, "cx_count", ""),
            get(transpile, "ecr_count", ""),
            get(transpile, "elapsed_sec", ""),
            metadata["simulation_status"],
            metadata["blocker"],
        ]))
    end
end

function main()
    output_name = String(strip(get(ENV, "DSMFG_DIRECT_FULL_QUBO_OUTPUT_DIR", DEFAULT_OUTPUT_DIR)))
    isempty(output_name) && error("DSMFG_DIRECT_FULL_QUBO_OUTPUT_DIR must not be empty")
    output_dir = isabspath(output_name) ? output_name : joinpath(STUDY_ROOT, output_name)
    fake_backend_class = String(strip(get(ENV, "DSMFG_DIRECT_FULL_QUBO_FAKE_BACKEND", "FakeFez")))
    isempty(fake_backend_class) && error("DSMFG_DIRECT_FULL_QUBO_FAKE_BACKEND must not be empty")
    transpile = parse_bool_env("DSMFG_DIRECT_FULL_QUBO_TRANSPILE")
    mkpath(output_dir)

    qubo = read_qubo_data(STUDY_ROOT)
    components = read_aux_components(joinpath(REDUCED_DIR, "auxiliary_components.csv"), qubo)
    top_flows = read_exact_top_flows(joinpath(REDUCED_DIR, "reduced_exact_top_flows.csv"))

    summary_path = joinpath(output_dir, "direct_full_qubo_summary.csv")
    open(summary_path, "w") do io
        println(io, csv_row(direct_summary_header()))
        for spec in source_specs()
            println(io, csv_row(score_distribution(spec, qubo, components, top_flows)))
        end
    end

    parameter_path = joinpath(output_dir, DIRECT_QAOA_PARAMETER_ARTIFACT)
    write_json_file(parameter_path, direct_qaoa_parameter_artifact(qubo))

    metadata = circuit_resource_metadata(qubo; transpile = transpile, fake_backend_class = fake_backend_class)
    metadata_path = joinpath(output_dir, "direct_full_qubo_resource_metadata.json")
    resource_summary_path = joinpath(output_dir, "direct_full_qubo_resource_summary.csv")
    write_json_file(metadata_path, metadata)
    write_resource_summary(resource_summary_path, metadata)

    println("Direct full-QUBO summary written to: ", summary_path)
    println("Direct full-QUBO QAOA hardware parameters written to: ", parameter_path)
    println("Direct full-QUBO resource summary written to: ", resource_summary_path)
    println("Direct full-QUBO resource metadata written to: ", metadata_path)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
