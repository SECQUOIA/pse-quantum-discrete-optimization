#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Dates
using Printf
using PythonCall

isdefined(@__MODULE__, :hit_rate_stat_values) || include(joinpath(@__DIR__, "hit_rate_stats.jl"))

const STUDY_ROOT = normpath(abspath(joinpath(@__DIR__, "..")))
const REDUCED_DIR = joinpath(STUDY_ROOT, "ds_mfg_reduced_flow_objective")
const ANGLE_DIR = joinpath(STUDY_ROOT, "ds_mfg_qaoa_juliqaoa_objective_angle_search")
const ANGLE_TARGET = "top10"
const ANGLE_P = 5
const ALGORITHM = "QAOA_reduced_surrogate_JuliQAOA_IBM_pilot"
const SCHEMA_VERSION = 1
const DEFAULT_IBM_RUNTIME_CHANNEL = "ibm_quantum_platform"

struct PilotConfig
    channel::String
    backend_name::String
    instance::Union{Nothing,String}
    final_reads::Int
    repeats::Int
    transpile_seeds::Vector{Int}
    output_dir::String
    run_hardware::Bool
end

function require_env(name::AbstractString)
    value = strip(get(ENV, name, ""))
    isempty(value) && error("Missing required environment variable $(name)")
    return value
end

function parse_positive_int_env(name::AbstractString, default::Integer)
    raw = strip(get(ENV, name, string(default)))
    value = parse(Int, raw)
    value > 0 || error("$(name) must be a positive integer, got $(raw)")
    return value
end

function parse_bool_env(name::AbstractString)
    raw = lowercase(strip(get(ENV, name, "false")))
    raw in ("", "false", "0", "no") && return false
    raw in ("true", "1", "yes") && return true
    error("$(name) must be true or false, got $(raw)")
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
    isempty(items) && error("At least one transpile seed is required")
    return items
end

function read_config()
    backend_name = require_env("QISKIT_IBM_BACKEND")
    channel = strip(get(ENV, "QISKIT_IBM_CHANNEL", DEFAULT_IBM_RUNTIME_CHANNEL))
    isempty(channel) && error("QISKIT_IBM_CHANNEL must not be empty")
    instance_text = strip(get(ENV, "QISKIT_IBM_INSTANCE", ""))
    output_name = strip(get(ENV, "DSMFG_HARDWARE_OUTPUT_DIR", "ds_mfg_ibm_qaoa_pilot"))
    isempty(output_name) && error("DSMFG_HARDWARE_OUTPUT_DIR must not be empty")
    output_dir = isabspath(output_name) ? output_name : joinpath(STUDY_ROOT, output_name)

    return PilotConfig(
        channel,
        backend_name,
        isempty(instance_text) ? nothing : instance_text,
        parse_positive_int_env("DSMFG_HARDWARE_FINAL_READS", 4096),
        parse_positive_int_env("DSMFG_HARDWARE_REPEATS", 1),
        parse_int_list(get(ENV, "DSMFG_HARDWARE_TRANSPILE_SEEDS", "92001")),
        normpath(abspath(output_dir)),
        parse_bool_env("DSMFG_RUN_IBM_HARDWARE"),
    )
end

function read_reduced_scalars(path::AbstractString)
    rows = readlines(path)
    length(rows) >= 2 || error("Missing scalar row in $(path)")
    values = split(rows[2], ','; keepempty = true)
    return (
        n = round(Int, parse(Float64, values[1])),
        scale = parse(Float64, values[2]),
        offset = parse(Float64, values[3]),
    )
end

function read_reduced_vector(path::AbstractString)
    values = Float64[]
    for row in readlines(path)[2:end]
        parts = split(row, ','; keepempty = true)
        push!(values, parse(Float64, parts[2]))
    end
    return values
end

function read_reduced_matrix(path::AbstractString)
    rows = readlines(path)[2:end]
    matrix = Matrix{Float64}(undef, length(rows), length(rows))
    for (i, row) in enumerate(rows)
        parts = split(row, ','; keepempty = true)
        for j in 1:length(rows)
            matrix[i, j] = parse(Float64, parts[j + 1])
        end
    end
    return matrix
end

function read_exact_top_flows(path::AbstractString)
    rows = readlines(path)[2:end]
    top = Dict{String,NamedTuple}()
    for row in rows
        isempty(strip(row)) && continue
        parts = split(row, ','; keepempty = true)
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

function parse_angle_summary(path::AbstractString)
    rows = readlines(path)
    length(rows) >= 2 || error("Missing angle rows in $(path)")
    header = split(rows[1], ','; keepempty = true)
    parsed = Dict{Tuple{String,Int},Dict{String,String}}()
    for (line_number, row) in enumerate(rows[2:end])
        isempty(strip(row)) && continue
        values = split(row, ','; keepempty = true)
        length(values) == length(header) ||
            error("Angle row $(line_number + 1) in $(path) has $(length(values)) fields; expected $(length(header))")
        record = Dict(header[i] => values[i] for i in eachindex(header))
        parsed[(record["target"], parse(Int, record["p"]))] = record
    end
    return parsed
end

function parse_angles(value::AbstractString)
    angles = [parse(Float64, strip(part)) for part in split(value, ';') if !isempty(strip(part))]
    iseven(length(angles)) || error("Expected beta-then-gamma QAOA angles, got an odd count")
    return angles
end

function load_problem_data()
    scalars = read_reduced_scalars(joinpath(REDUCED_DIR, "reduced_scalars.csv"))
    L = read_reduced_vector(joinpath(REDUCED_DIR, "reduced_L_vector.csv"))
    Q = read_reduced_matrix(joinpath(REDUCED_DIR, "reduced_Q_matrix.csv"))
    top_flows = read_exact_top_flows(joinpath(REDUCED_DIR, "reduced_exact_top_flows.csv"))
    angle_path = joinpath(ANGLE_DIR, "juliqaoa_angle_summary.csv")
    angle_rows = parse_angle_summary(angle_path)
    angle_record = get(angle_rows, (ANGLE_TARGET, ANGLE_P), nothing)
    isnothing(angle_record) && error("Missing $(ANGLE_TARGET) p=$(ANGLE_P) angle row in $(angle_path)")
    angles = parse_angles(angle_record["qiskit_angles_beta_then_gamma"])
    length(angles) == 2 * ANGLE_P ||
        error("Expected $(2 * ANGLE_P) angles for p=$(ANGLE_P), got $(length(angles))")

    scalars.n == length(L) == size(Q, 1) == size(Q, 2) ||
        error("Reduced surrogate dimensions are inconsistent")

    return (
        scalars = scalars,
        L = L,
        Q = Q,
        top_flows = top_flows,
        angle_path = angle_path,
        angle_record = angle_record,
        angles = angles,
    )
end

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

utc_timestamp() = string(Dates.now(Dates.UTC), "Z")

function effective_qubo_terms(L::Vector{Float64}, Q::Matrix{Float64}; scale::Float64)
    n = length(L)
    linear = [scale * (L[i] + Q[i, i]) for i in 1:n]
    pairs = Tuple{Int,Int,Float64}[]
    for i in 1:(n - 1), j in (i + 1):n
        coefficient = scale * (Q[i, j] + Q[j, i])
        abs(coefficient) > 1.0e-12 && push!(pairs, (i, j, coefficient))
    end
    return linear, pairs
end

function ising_terms_from_qubo(L::Vector{Float64}, Q::Matrix{Float64}; scale::Float64)
    linear, qubo_pairs = effective_qubo_terms(L, Q; scale = scale)
    h = zeros(Float64, length(L))
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

function build_qaoa_circuit(data)
    qiskit = pyimport("qiskit")
    n = data.scalars.n
    p = length(data.angles) ÷ 2
    betas = data.angles[1:p]
    gammas = data.angles[(p + 1):end]
    h, zz_pairs = ising_terms_from_qubo(data.L, data.Q; scale = data.scalars.scale)

    circuit = qiskit.QuantumCircuit(n, n)
    for qubit in 0:(n - 1)
        circuit.h(qubit)
    end

    for layer in 1:p
        gamma = gammas[layer]
        for i in 1:n
            theta = 2.0 * gamma * h[i]
            abs(theta) > 1.0e-12 && circuit.rz(theta, i - 1)
        end
        for (i, j, coefficient) in zz_pairs
            theta = 2.0 * gamma * coefficient
            abs(theta) > 1.0e-12 && circuit.rzz(theta, i - 1, j - 1)
        end

        beta = betas[layer]
        for qubit in 0:(n - 1)
            circuit.rx(2.0 * beta, qubit)
        end
    end

    for qubit in 0:(n - 1)
        circuit.measure(qubit, qubit)
    end
    return circuit
end

function py_dict_to_julia(py_dict)
    try
        return Dict{String,Int}(pyconvert(Dict{String,Int}, py_dict))
    catch
        converted = Dict{String,Int}()
        for item in py_dict.items()
            converted[pyconvert(String, item[0])] = pyconvert(Int, item[1])
        end
        return converted
    end
end

function circuit_metadata(circuit)
    return Dict{String,Any}(
        "num_qubits" => pyconvert(Int, circuit.num_qubits),
        "num_clbits" => pyconvert(Int, circuit.num_clbits),
        "depth" => pyconvert(Int, circuit.depth()),
        "operation_counts" => py_dict_to_julia(circuit.count_ops()),
        "cost_unitary" => "exp(-i gamma C(x)) for the reduced surrogate QUBO",
        "mixer_unitary" => "exp(-i beta sum X)",
        "measurement_bit_order" => "Qiskit count keys are reversed before scoring as x1..x19 flow bits",
    )
end

function normalize_qiskit_key(key::AbstractString, n::Integer)
    compact = replace(strip(key), " " => "")
    length(compact) == n || error("Expected $(n)-bit Qiskit count key, got $(key)")
    all(char -> char == '0' || char == '1', compact) ||
        error("Qiskit count key contains non-binary characters: $(key)")
    return compact
end

function qiskit_key_to_flow_bits(key::AbstractString, n::Integer)
    return reverse(normalize_qiskit_key(key, n))
end

function bits_vector(bits::AbstractString)
    return [char == '1' ? 1 : 0 for char in bits]
end

function surrogate_energy(flow_bits::AbstractString, data)
    bits = bits_vector(flow_bits)
    energy = data.scalars.offset
    for i in eachindex(bits)
        energy += data.L[i] * bits[i]
    end
    for i in eachindex(bits), j in eachindex(bits)
        energy += data.Q[i, j] * bits[i] * bits[j]
    end
    return data.scalars.scale * energy
end

function planned_jobs(config::PilotConfig)
    jobs = Vector{Dict{String,Any}}()
    job_index = 1
    for repeat in 1:config.repeats, seed in config.transpile_seeds
        push!(
            jobs,
            Dict{String,Any}(
                "job_index" => job_index,
                "repeat" => repeat,
                "transpile_seed" => seed,
                "shots" => config.final_reads,
                "submitted" => false,
                "job_id" => nothing,
                "status" => "planned",
            ),
        )
        job_index += 1
    end
    return jobs
end

function output_paths(config::PilotConfig)
    return Dict{String,String}(
        "job_manifest" => joinpath(config.output_dir, "job_manifest.json"),
        "backend_metadata" => joinpath(config.output_dir, "backend_metadata.json"),
        "raw_counts" => joinpath(config.output_dir, "raw_counts.csv"),
        "scored_counts" => joinpath(config.output_dir, "scored_counts.csv"),
        "summary" => joinpath(config.output_dir, "summary.csv"),
    )
end

function angle_metadata(data)
    record = data.angle_record
    return Dict{String,Any}(
        "path" => relpath(data.angle_path, STUDY_ROOT),
        "target" => ANGLE_TARGET,
        "p" => ANGLE_P,
        "seed" => parse(Int, record["seed"]),
        "basinhopping_niter" => parse(Int, record["basinhopping_niter"]),
        "parameter_order" => "beta_then_gamma",
        "qiskit_angles_beta_then_gamma" => data.angles,
        "top50_probability" => parse(Float64, record["top50_probability"]),
        "top10_probability" => parse(Float64, record["top10_probability"]),
        "global_probability" => parse(Float64, record["global_probability"]),
    )
end

function manifest(config::PilotConfig, data, circuit_info, jobs, paths)
    return Dict{String,Any}(
        "schema_version" => SCHEMA_VERSION,
        "created_at_utc" => utc_timestamp(),
        "script" => relpath(@__FILE__, STUDY_ROOT),
        "mode" => config.run_hardware ? "hardware" : "dry_run",
        "configuration" => Dict{String,Any}(
            "channel" => config.channel,
            "backend_name" => config.backend_name,
            "qiskit_ibm_instance_configured" => !isnothing(config.instance),
            "final_reads" => config.final_reads,
            "repeats" => config.repeats,
            "transpile_seeds" => config.transpile_seeds,
            "output_dir" => config.output_dir,
            "run_hardware" => config.run_hardware,
        ),
        "problem" => Dict{String,Any}(
            "name" => "DS-MFG reduced 19-flow surrogate",
            "n_qubits" => data.scalars.n,
            "scale" => data.scalars.scale,
            "offset" => data.scalars.offset,
            "reduced_dir" => relpath(REDUCED_DIR, STUDY_ROOT),
        ),
        "angle_source" => angle_metadata(data),
        "circuit" => circuit_info,
        "jobs" => jobs,
        "output_files" => Dict(key => relpath(path, config.output_dir) for (key, path) in paths),
        "credential_policy" => "No IBM tokens, account files, credential paths, or backend secrets are written by this script.",
    )
end

function persist_job_manifest!(config::PilotConfig, data, circuit_info, jobs, paths)
    write_json_file(paths["job_manifest"], manifest(config, data, circuit_info, jobs, paths))
    return nothing
end

function dry_run_backend_metadata(config::PilotConfig)
    return Dict{String,Any}(
        "schema_version" => SCHEMA_VERSION,
        "created_at_utc" => utc_timestamp(),
        "mode" => "dry_run",
        "channel" => config.channel,
        "backend_name" => config.backend_name,
        "queried" => false,
        "reason" => "Dry run does not contact IBM Quantum Runtime.",
        "qiskit_ibm_instance_configured" => !isnothing(config.instance),
        "credential_fields_written" => String[],
    )
end

function maybe_py_attr(obj, attr::AbstractString)
    pyhasattr(obj, attr) || return nothing
    return pygetattr(obj, attr)
end

function real_backend_metadata(config::PilotConfig, backend)
    metadata = Dict{String,Any}(
        "schema_version" => SCHEMA_VERSION,
        "created_at_utc" => utc_timestamp(),
        "mode" => "hardware",
        "channel" => config.channel,
        "backend_name_requested" => config.backend_name,
        "queried" => true,
        "qiskit_ibm_instance_configured" => !isnothing(config.instance),
        "credential_fields_written" => String[],
    )

    name_attr = maybe_py_attr(backend, "name")
    isnothing(name_attr) || (metadata["backend_name_resolved"] = pyconvert(String, name_attr))

    num_qubits = maybe_py_attr(backend, "num_qubits")
    isnothing(num_qubits) || (metadata["num_qubits"] = pyconvert(Int, num_qubits))

    target = maybe_py_attr(backend, "target")
    if !isnothing(target) && pyhasattr(target, "operation_names")
        metadata["operation_names"] = sort(pyconvert(Vector{String}, target.operation_names))
    end

    if pyhasattr(backend, "configuration")
        try
            configuration = backend.configuration()
            if pyhasattr(configuration, "basis_gates")
                metadata["basis_gates"] = sort(pyconvert(Vector{String}, configuration.basis_gates))
            end
            if pyhasattr(configuration, "simulator")
                metadata["simulator"] = pyconvert(Bool, configuration.simulator)
            end
        catch
            metadata["configuration_warning"] = "backend.configuration() metadata unavailable"
        end
    end

    if pyhasattr(backend, "status")
        try
            status = backend.status()
            pyhasattr(status, "operational") && (metadata["operational"] = pyconvert(Bool, status.operational))
            pyhasattr(status, "pending_jobs") && (metadata["pending_jobs"] = pyconvert(Int, status.pending_jobs))
        catch
            metadata["status_warning"] = "backend.status() metadata unavailable"
        end
    end

    return metadata
end

function write_empty_count_files(paths)
    open(paths["raw_counts"], "w") do io
        println(io, csv_row([
            "algorithm", "mode", "backend", "repeat", "transpile_seed", "job_id",
            "qiskit_bitstring", "flow_bits", "reads",
        ]))
    end
    open(paths["scored_counts"], "w") do io
        println(io, csv_row([
            "algorithm", "mode", "backend", "repeat", "transpile_seed", "job_id",
            "flow_bits", "reads", "probability", "surrogate_qubo_energy",
            "top50_rank", "exact_repaired_qubo_energy", "match", "matched_ip_obj",
            "repaired_aux_bits", "repaired_full_bits",
        ]))
    end
end

function summary_header()
    return [
        "algorithm", "mode", "backend", "angle_target", "p", "final_reads",
        "repeats", "transpile_seeds", "submitted_jobs", "elapsed_sec",
        "unique_flows", "total_reads", "top50_hits", "top10_hits", "global_hits",
        "gurobi_pool_feasible_hits", hit_rate_stat_headers()..., "best_top50_rank",
        "best_top50_exact_repaired_qubo_energy", "best_top50_match",
        "best_top50_flow_bits", "manifest_path", "backend_metadata_path",
        "raw_counts_path", "scored_counts_path",
    ]
end

function summarize_counts(aggregate::Dict{String,Int}, data; elapsed_sec::Real)
    total_reads = sum(values(aggregate); init = 0)
    top50_hits = 0
    top10_hits = 0
    global_hits = 0
    gurobi_pool_feasible_hits = 0
    best_hit = nothing

    for (flow_bits, reads) in aggregate
        hit = get(data.top_flows, flow_bits, nothing)
        isnothing(hit) && continue
        top50_hits += reads
        hit.rank <= 10 && (top10_hits += reads)
        hit.rank == 1 && (global_hits += reads)
        is_gurobi_pool_feasible_match(hit.match) && (gurobi_pool_feasible_hits += reads)
        if isnothing(best_hit) || hit.rank < best_hit.rank
            best_hit = merge(hit, (flow_bits = flow_bits,))
        end
    end

    hit_stats = hit_rate_stat_values(
        total_reads,
        top50_hits,
        top10_hits,
        global_hits,
        gurobi_pool_feasible_hits,
        elapsed_sec,
    )

    return (
        unique_flows = length(aggregate),
        total_reads = total_reads,
        top50_hits = top50_hits,
        top10_hits = top10_hits,
        global_hits = global_hits,
        gurobi_pool_feasible_hits = gurobi_pool_feasible_hits,
        hit_stats = hit_stats,
        best_hit = best_hit,
    )
end

function write_summary(path, config::PilotConfig, paths, aggregate::Dict{String,Int}; submitted_jobs::Int, elapsed_sec::Real, data)
    summary = summarize_counts(aggregate, data; elapsed_sec = elapsed_sec)
    best_rank = isnothing(summary.best_hit) ? "" : summary.best_hit.rank
    best_energy = isnothing(summary.best_hit) ? "" : summary.best_hit.exact_repaired_qubo_energy
    best_match = isnothing(summary.best_hit) ? "" : summary.best_hit.match
    best_flow = isnothing(summary.best_hit) ? "" : summary.best_hit.flow_bits

    open(path, "w") do io
        println(io, csv_row(summary_header()))
        println(
            io,
            csv_row([
                ALGORITHM,
                config.run_hardware ? "hardware" : "dry_run",
                config.backend_name,
                ANGLE_TARGET,
                ANGLE_P,
                config.final_reads,
                config.repeats,
                join(config.transpile_seeds, ';'),
                submitted_jobs,
                @sprintf("%.6f", elapsed_sec),
                summary.unique_flows,
                summary.total_reads,
                summary.top50_hits,
                summary.top10_hits,
                summary.global_hits,
                summary.gurobi_pool_feasible_hits,
                summary.hit_stats...,
                best_rank,
                best_energy,
                best_match,
                best_flow,
                relpath(paths["job_manifest"], config.output_dir),
                relpath(paths["backend_metadata"], config.output_dir),
                relpath(paths["raw_counts"], config.output_dir),
                relpath(paths["scored_counts"], config.output_dir),
            ]),
        )
    end
end

function sampler_counts(result)
    pub_result = result[0]
    data = pub_result.data
    if pyhasattr(data, "c")
        return py_dict_to_julia(data.c.get_counts())
    elseif pyhasattr(data, "meas")
        return py_dict_to_julia(data.meas.get_counts())
    elseif pyhasattr(pub_result, "join_data")
        return py_dict_to_julia(pub_result.join_data().get_counts())
    end
    error("Could not find sampler count data in Runtime SamplerV2 result")
end

function write_count_rows!(raw_io, scored_io, aggregate::Dict{String,Int}, config::PilotConfig, data, job, counts)
    total_reads = sum(values(counts); init = 0)
    n = data.scalars.n

    for (qiskit_bitstring, reads) in sort(collect(counts); by = first)
        flow_bits = qiskit_key_to_flow_bits(qiskit_bitstring, n)
        aggregate[flow_bits] = get(aggregate, flow_bits, 0) + reads
        probability = total_reads == 0 ? 0.0 : reads / total_reads
        repaired = get(data.top_flows, flow_bits, nothing)

        println(
            raw_io,
            csv_row([
                ALGORITHM,
                "hardware",
                config.backend_name,
                job["repeat"],
                job["transpile_seed"],
                job["job_id"],
                qiskit_bitstring,
                flow_bits,
                reads,
            ]),
        )

        println(
            scored_io,
            csv_row([
                ALGORITHM,
                "hardware",
                config.backend_name,
                job["repeat"],
                job["transpile_seed"],
                job["job_id"],
                flow_bits,
                reads,
                probability,
                surrogate_energy(flow_bits, data),
                isnothing(repaired) ? "" : repaired.rank,
                isnothing(repaired) ? "" : repaired.exact_repaired_qubo_energy,
                isnothing(repaired) ? "" : repaired.match,
                isnothing(repaired) ? "" : repaired.matched_ip_obj,
                isnothing(repaired) ? "" : repaired.repaired_aux_bits,
                isnothing(repaired) ? "" : repaired.repaired_full_bits,
            ]),
        )
    end
end

function job_status_text(job)
    if pyhasattr(job, "status")
        try
            return string(job.status())
        catch
            return "unknown"
        end
    end
    return "unknown"
end

function record_job_failure!(job_plan, runtime_job, submitted::Bool, err)
    if runtime_job !== nothing
        job_plan["status"] = job_status_text(runtime_job)
    elseif submitted
        job_plan["status"] = get(job_plan, "status", "submitted")
    else
        job_plan["status"] = "failed_before_submission"
    end
    job_plan["error"] = submitted ?
        "Runtime job result retrieval or scoring failed after submission." :
        "Runtime job transpile or submission failed before a job ID was assigned."
    job_plan["error_type"] = string(typeof(err))
    job_plan["error_recorded_at_utc"] = utc_timestamp()
    return nothing
end

function runtime_service(runtime, config::PilotConfig)
    token = strip(get(ENV, "QISKIT_IBM_TOKEN", ""))
    if isempty(token)
        if isnothing(config.instance)
            return runtime.QiskitRuntimeService(; channel = config.channel)
        else
            return runtime.QiskitRuntimeService(; channel = config.channel, instance = config.instance)
        end
    else
        if isnothing(config.instance)
            return runtime.QiskitRuntimeService(; channel = config.channel, token = token)
        else
            return runtime.QiskitRuntimeService(; channel = config.channel, token = token, instance = config.instance)
        end
    end
end

function runtime_service_or_error(runtime, config::PilotConfig, paths)
    original_error_type = nothing
    try
        return runtime_service(runtime, config)
    catch err
        original_error_type = typeof(err)
    end
    write_json_file(
        paths["backend_metadata"],
        Dict{String,Any}(
            "schema_version" => SCHEMA_VERSION,
            "created_at_utc" => utc_timestamp(),
            "mode" => "hardware",
            "channel" => config.channel,
            "backend_name_requested" => config.backend_name,
            "queried" => false,
            "qiskit_ibm_instance_configured" => !isnothing(config.instance),
            "credential_fields_written" => String[],
            "error" => "IBM Runtime service could not be initialized. Set QISKIT_IBM_INSTANCE to a valid Runtime instance CRN when the account cannot be auto-resolved.",
        ),
    )
    error(
        "IBM Runtime service could not be initialized. ",
        "Set QISKIT_IBM_INSTANCE to a valid Runtime instance CRN when the account cannot be auto-resolved. ",
        "Backend metadata written to $(paths["backend_metadata"]). Original error type: $(original_error_type).",
    )
end

function run_hardware_jobs!(config::PilotConfig, data, circuit, circuit_info, jobs, paths)
    qiskit = pyimport("qiskit")
    runtime = pyimport("qiskit_ibm_runtime")
    service = runtime_service_or_error(runtime, config, paths)
    backend = service.backend(config.backend_name)
    write_json_file(paths["backend_metadata"], real_backend_metadata(config, backend))

    aggregate = Dict{String,Int}()
    submitted_jobs = 0

    open(paths["raw_counts"], "w") do raw_io
        println(raw_io, csv_row([
            "algorithm", "mode", "backend", "repeat", "transpile_seed", "job_id",
            "qiskit_bitstring", "flow_bits", "reads",
        ]))
        open(paths["scored_counts"], "w") do scored_io
            println(scored_io, csv_row([
                "algorithm", "mode", "backend", "repeat", "transpile_seed", "job_id",
                "flow_bits", "reads", "probability", "surrogate_qubo_energy",
                "top50_rank", "exact_repaired_qubo_energy", "match", "matched_ip_obj",
                "repaired_aux_bits", "repaired_full_bits",
            ]))

            for job_plan in jobs
                runtime_job = nothing
                submitted = false
                try
                    @info "Submitting fixed-parameter QAOA sampler job" repeat=job_plan["repeat"] transpile_seed=job_plan["transpile_seed"] shots=config.final_reads
                    transpiled = qiskit.transpile(
                        circuit;
                        backend = backend,
                        seed_transpiler = job_plan["transpile_seed"],
                        optimization_level = 3,
                    )
                    sampler = runtime.SamplerV2(; mode = backend)
                    runtime_job = sampler.run([transpiled]; shots = config.final_reads)
                    submitted = true
                    job_plan["submitted"] = true
                    job_plan["job_id"] = pyconvert(String, runtime_job.job_id())
                    job_plan["status"] = job_status_text(runtime_job)
                    submitted_jobs += 1
                    persist_job_manifest!(config, data, circuit_info, jobs, paths)

                    result = runtime_job.result()
                    counts = sampler_counts(result)
                    job_plan["status"] = job_status_text(runtime_job)
                    job_plan["unique_raw_bitstrings"] = length(counts)
                    job_plan["total_reads"] = sum(values(counts); init = 0)
                    persist_job_manifest!(config, data, circuit_info, jobs, paths)
                    write_count_rows!(raw_io, scored_io, aggregate, config, data, job_plan, counts)
                    flush(raw_io)
                    flush(scored_io)
                catch err
                    record_job_failure!(job_plan, runtime_job, submitted || get(job_plan, "submitted", false), err)
                    persist_job_manifest!(config, data, circuit_info, jobs, paths)
                    rethrow()
                end
            end
        end
    end

    return aggregate, submitted_jobs
end

function main()
    config = read_config()
    data = load_problem_data()
    mkpath(config.output_dir)

    circuit = build_qaoa_circuit(data)
    circuit_info = circuit_metadata(circuit)
    paths = output_paths(config)
    jobs = planned_jobs(config)
    aggregate = Dict{String,Int}()
    submitted_jobs = 0
    start_time = time()

    persist_job_manifest!(config, data, circuit_info, jobs, paths)

    if config.run_hardware
        aggregate, submitted_jobs = run_hardware_jobs!(config, data, circuit, circuit_info, jobs, paths)
    else
        write_json_file(paths["backend_metadata"], dry_run_backend_metadata(config))
        write_empty_count_files(paths)
    end

    elapsed_sec = time() - start_time
    write_summary(
        paths["summary"],
        config,
        paths,
        aggregate;
        submitted_jobs = submitted_jobs,
        elapsed_sec = elapsed_sec,
        data = data,
    )
    persist_job_manifest!(config, data, circuit_info, jobs, paths)

    if config.run_hardware
        println("IBM QAOA pilot hardware run complete.")
    else
        println("IBM QAOA pilot dry run complete. No hardware jobs were submitted.")
    end
    println("Output directory: ", config.output_dir)
    println("Manifest: ", paths["job_manifest"])
    println("Summary: ", paths["summary"])
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
