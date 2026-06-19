#!/usr/bin/env julia

import Pkg

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(joinpath(@__DIR__, ".."))
    Pkg.instantiate()
end

using Dates
using Printf
using PythonCall

include(joinpath(@__DIR__, "run_direct_full_qubo_audit.jl"))

const DIRECT_HARDWARE_ALGORITHM = "QAOA_direct_full_qubo_IBM_handoff"
const DIRECT_HARDWARE_SCHEMA_VERSION = 1
const DIRECT_HARDWARE_DEFAULT_OUTPUT_DIR = "ds_mfg_direct_full_qubo_hardware_pilot"
const DIRECT_HARDWARE_DEFAULT_CHANNEL = "ibm_quantum_platform"
const DIRECT_FLOW_VARIABLES = 19

struct DirectHardwareConfig
    channel::String
    backend_name::String
    instance::Union{Nothing,String}
    final_reads::Int
    repeats::Int
    transpile_seeds::Vector{Int}
    output_dir::String
    run_hardware::Bool
    parameter_path::String
end

struct DirectQaoaParameters
    path::String
    p::Int
    angles::Vector{Float64}
    angle_source::String
    source_summary_artifact::String
    source_distribution_artifact::String
    hardware_handoff_policy::String
end

direct_utc_timestamp() = string(Dates.now(Dates.UTC), "Z")

function direct_require_env(name::AbstractString)
    value = strip(get(ENV, name, ""))
    isempty(value) && error("Missing required environment variable $(name)")
    return value
end

function direct_parse_positive_int_env(name::AbstractString, default::Integer)
    raw = strip(get(ENV, name, string(default)))
    value = parse(Int, raw)
    value > 0 || error("$(name) must be a positive integer, got $(raw)")
    return value
end

function direct_parse_int_list(value::AbstractString)
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

function direct_default_parameter_path()
    return joinpath(STUDY_ROOT, DEFAULT_OUTPUT_DIR, DIRECT_QAOA_PARAMETER_ARTIFACT)
end

function direct_abs_path(raw_path::AbstractString)
    path = strip(raw_path)
    isempty(path) && error("Path value must not be empty")
    return normpath(abspath(isabspath(path) ? path : joinpath(STUDY_ROOT, path)))
end

function read_direct_hardware_config()
    backend_name = direct_require_env("QISKIT_IBM_BACKEND")
    channel = strip(get(ENV, "QISKIT_IBM_CHANNEL", DIRECT_HARDWARE_DEFAULT_CHANNEL))
    isempty(channel) && error("QISKIT_IBM_CHANNEL must not be empty")
    instance_text = strip(get(ENV, "QISKIT_IBM_INSTANCE", ""))
    output_name = strip(get(ENV, "DSMFG_DIRECT_HARDWARE_OUTPUT_DIR", DIRECT_HARDWARE_DEFAULT_OUTPUT_DIR))
    isempty(output_name) && error("DSMFG_DIRECT_HARDWARE_OUTPUT_DIR must not be empty")
    output_dir = isabspath(output_name) ? output_name : joinpath(STUDY_ROOT, output_name)
    parameter_path = direct_abs_path(get(
        ENV,
        "DSMFG_DIRECT_QAOA_PARAMETER_PATH",
        direct_default_parameter_path(),
    ))

    return DirectHardwareConfig(
        channel,
        backend_name,
        isempty(instance_text) ? nothing : instance_text,
        direct_parse_positive_int_env("DSMFG_DIRECT_HARDWARE_FINAL_READS", 4096),
        direct_parse_positive_int_env("DSMFG_DIRECT_HARDWARE_REPEATS", 1),
        direct_parse_int_list(get(ENV, "DSMFG_DIRECT_HARDWARE_TRANSPILE_SEEDS", "92001")),
        normpath(abspath(output_dir)),
        parse_bool_env("DSMFG_RUN_DIRECT_FULL_QUBO_HARDWARE"),
        parameter_path,
    )
end

function direct_py_string(py_dict, key::AbstractString; default = "")
    try
        return pyconvert(String, py_dict[key])
    catch
        return default
    end
end

function read_direct_qaoa_parameters(path::AbstractString, qubo::QuboData)
    isfile(path) || error(
        "Missing direct QAOA parameter artifact at $(path). ",
        "Run scripts/run_direct_full_qubo_audit.jl first or set DSMFG_DIRECT_QAOA_PARAMETER_PATH.",
    )

    parsed = pyimport("json").loads(read(path, String))
    schema_version = pyconvert(Int, parsed["schema_version"])
    schema_version == 1 || error("Unsupported direct QAOA parameter schema version: $(schema_version)")
    role = pyconvert(String, parsed["artifact_role"])
    role == "direct_full_qubo_qaoa_hardware_parameter_handoff" ||
        error("Unexpected direct QAOA parameter artifact role: $(role)")
    n_qubits = pyconvert(Int, parsed["n_qubits"])
    n_qubits == qubo.n || error("Parameter artifact uses $(n_qubits) qubits; expected $(qubo.n)")
    n_flow = pyconvert(Int, parsed["n_flow_variables"])
    n_flow == DIRECT_FLOW_VARIABLES ||
        error("Parameter artifact uses $(n_flow) flow variables; expected $(DIRECT_FLOW_VARIABLES)")
    order = pyconvert(String, parsed["parameter_order"])
    order == "beta_then_gamma" || error("Expected beta_then_gamma QAOA parameter order, got $(order)")

    p = pyconvert(Int, parsed["p"])
    angles = pyconvert(Vector{Float64}, parsed["qiskit_angles_beta_then_gamma"])
    validate_direct_qaoa_parameters(p, angles)

    return DirectQaoaParameters(
        normpath(abspath(path)),
        p,
        angles,
        direct_py_string(parsed, "angle_source"; default = DIRECT_QAOA_ANGLE_SOURCE),
        direct_py_string(parsed, "source_summary_artifact"; default = ""),
        direct_py_string(parsed, "source_distribution_artifact"; default = ""),
        direct_py_string(parsed, "hardware_handoff_policy"; default = ""),
    )
end

function direct_py_dict_to_julia(py_dict)
    converted = Dict{String,Int}()
    for item in py_dict.items()
        converted[pyconvert(String, item[0])] = pyconvert(Int, item[1])
    end
    return converted
end

function direct_circuit_metadata(circuit, zz_pairs)
    operation_counts = direct_py_dict_to_julia(circuit.count_ops())
    rzz_count = get(operation_counts, "rzz", 0)
    cx_count = get(operation_counts, "cx", 0)
    return Dict{String,Any}(
        "num_qubits" => pyconvert(Int, circuit.num_qubits),
        "num_clbits" => pyconvert(Int, circuit.num_clbits),
        "depth" => pyconvert(Int, circuit.depth()),
        "operation_counts" => operation_counts,
        "nonzero_ising_zz_pairs" => length(zz_pairs),
        "logical_two_qubit_ops" => rzz_count + cx_count,
        "measurement_bit_order" => "Qiskit count keys are reversed before scoring as x1..x36 full QUBO bits",
    )
end

function direct_planned_jobs(config::DirectHardwareConfig)
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

function direct_output_paths(config::DirectHardwareConfig)
    return Dict{String,String}(
        "job_manifest" => joinpath(config.output_dir, "job_manifest.json"),
        "backend_metadata" => joinpath(config.output_dir, "backend_metadata.json"),
        "raw_counts" => joinpath(config.output_dir, "raw_counts.csv"),
        "scored_counts" => joinpath(config.output_dir, "scored_counts.csv"),
        "summary" => joinpath(config.output_dir, "summary.csv"),
    )
end

function direct_parameter_metadata(params::DirectQaoaParameters)
    return Dict{String,Any}(
        "path" => relpath(params.path, STUDY_ROOT),
        "angle_source" => params.angle_source,
        "source_summary_artifact" => params.source_summary_artifact,
        "source_distribution_artifact" => params.source_distribution_artifact,
        "p" => params.p,
        "parameter_order" => "beta_then_gamma",
        "qiskit_angles_beta_then_gamma" => params.angles,
        "hardware_handoff_policy" => params.hardware_handoff_policy,
    )
end

function direct_manifest(
    config::DirectHardwareConfig,
    qubo::QuboData,
    params::DirectQaoaParameters,
    circuit_info,
    jobs,
    paths,
)
    return Dict{String,Any}(
        "schema_version" => DIRECT_HARDWARE_SCHEMA_VERSION,
        "created_at_utc" => direct_utc_timestamp(),
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
            "parameter_path" => relpath(config.parameter_path, STUDY_ROOT),
        ),
        "problem" => Dict{String,Any}(
            "name" => "DS-MFG direct original 36-variable QUBO",
            "data_source" => ZIP_NAME,
            "n_qubits" => qubo.n,
            "n_flow_variables" => DIRECT_FLOW_VARIABLES,
            "n_auxiliary_variables" => qubo.n - DIRECT_FLOW_VARIABLES,
            "scale" => qubo.scale,
            "offset" => qubo.offset,
        ),
        "qaoa_parameters" => direct_parameter_metadata(params),
        "circuit" => circuit_info,
        "jobs" => jobs,
        "output_files" => Dict(key => relpath(path, config.output_dir) for (key, path) in paths),
        "credential_policy" => "No IBM tokens, account files, credential paths, or backend secrets are written by this script.",
    )
end

function direct_persist_job_manifest!(
    config::DirectHardwareConfig,
    qubo::QuboData,
    params::DirectQaoaParameters,
    circuit_info,
    jobs,
    paths,
)
    write_json_file(paths["job_manifest"], direct_manifest(config, qubo, params, circuit_info, jobs, paths))
    return nothing
end

function direct_dry_run_backend_metadata(config::DirectHardwareConfig)
    return Dict{String,Any}(
        "schema_version" => DIRECT_HARDWARE_SCHEMA_VERSION,
        "created_at_utc" => direct_utc_timestamp(),
        "mode" => "dry_run",
        "channel" => config.channel,
        "backend_name" => config.backend_name,
        "queried" => false,
        "reason" => "Dry run does not contact IBM Quantum Runtime.",
        "qiskit_ibm_instance_configured" => !isnothing(config.instance),
        "credential_fields_written" => String[],
    )
end

function direct_maybe_py_attr(obj, attr::AbstractString)
    pyhasattr(obj, attr) || return nothing
    return pygetattr(obj, attr)
end

function direct_real_backend_metadata(config::DirectHardwareConfig, backend)
    metadata = Dict{String,Any}(
        "schema_version" => DIRECT_HARDWARE_SCHEMA_VERSION,
        "created_at_utc" => direct_utc_timestamp(),
        "mode" => "hardware",
        "channel" => config.channel,
        "backend_name_requested" => config.backend_name,
        "queried" => true,
        "qiskit_ibm_instance_configured" => !isnothing(config.instance),
        "credential_fields_written" => String[],
    )

    name_attr = direct_maybe_py_attr(backend, "name")
    isnothing(name_attr) || (metadata["backend_name_resolved"] = pyconvert(String, name_attr))

    num_qubits = direct_maybe_py_attr(backend, "num_qubits")
    isnothing(num_qubits) || (metadata["num_qubits"] = pyconvert(Int, num_qubits))

    target = direct_maybe_py_attr(backend, "target")
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

function direct_write_empty_count_files(paths)
    open(paths["raw_counts"], "w") do io
        println(io, csv_row([
            "algorithm", "mode", "backend", "repeat", "transpile_seed", "job_id",
            "qiskit_bitstring", "full_bits", "flow_bits", "reads",
        ]))
    end
    open(paths["scored_counts"], "w") do io
        println(io, csv_row([
            "algorithm", "mode", "backend", "repeat", "transpile_seed", "job_id",
            "full_bits", "flow_bits", "reads", "probability", "raw_qubo_energy",
            "exact_repaired_qubo_energy", "top50_rank", "match", "matched_ip_obj",
            "repaired_aux_bits", "repaired_full_bits", "encoded_repaired_match",
        ]))
    end
end

function direct_summary_header()
    return [
        "algorithm", "mode", "backend", "parameter_artifact", "p", "final_reads",
        "repeats", "transpile_seeds", "submitted_jobs", "elapsed_sec",
        "unique_full_states", "unique_flows", "total_reads", "projected_pool_hits",
        "encoded_pool_hits", "encoded_global_hits", "top50_hits", "top10_hits",
        "global_hits", "gurobi_pool_feasible_hits", hit_rate_stat_headers()...,
        "best_repaired_rank", "best_repaired_qubo_energy", "best_repaired_match",
        "best_repaired_flow_bits", "best_repaired_full_bits", "manifest_path",
        "backend_metadata_path", "raw_counts_path", "scored_counts_path",
    ]
end

function normalize_direct_qiskit_key(key::AbstractString, n::Integer)
    compact = replace(strip(key), " " => "")
    length(compact) == n || error("Expected $(n)-bit Qiskit count key, got $(key)")
    all(char -> char == '0' || char == '1', compact) ||
        error("Qiskit count key contains non-binary characters: $(key)")
    return compact
end

function qiskit_key_to_full_bits(key::AbstractString, n::Integer)
    return reverse(normalize_direct_qiskit_key(key, n))
end

function flow_bits_from_full_bits(full_bits::AbstractString)
    length(full_bits) >= DIRECT_FLOW_VARIABLES ||
        error("Expected at least $(DIRECT_FLOW_VARIABLES) full bits, got $(full_bits)")
    return full_bits[1:DIRECT_FLOW_VARIABLES]
end

function direct_full_qubo_energy(full_bits::AbstractString, qubo::QuboData)
    bits = bits_vector(full_bits)
    length(bits) == qubo.n || error("Expected $(qubo.n) bits, got $(length(bits))")
    energy = qubo.offset
    for i in eachindex(bits)
        bits[i] == 1 && (energy += qubo.linear[i])
    end
    for i in eachindex(bits), j in eachindex(bits)
        bits[i] == 1 && bits[j] == 1 && (energy += qubo.quadratic[i, j])
    end
    return qubo.scale * energy
end

function direct_summarize_counts(
    aggregate_full::Dict{String,Int},
    qubo::QuboData,
    components::Vector{AuxComponent},
    top_flows::Dict{String,NamedTuple};
    elapsed_sec::Real,
)
    total_reads = sum(values(aggregate_full); init = 0)
    flow_counts = Dict{String,Int}()
    top50_hits = 0
    top10_hits = 0
    global_hits = 0
    projected_pool_hits = 0
    encoded_pool_hits = 0
    encoded_global_hits = 0
    repair_cache = Dict{String,NamedTuple}()
    best_repaired = nothing
    best_repaired_flow_bits = ""

    for (full_bits, reads) in aggregate_full
        flow_bits = flow_bits_from_full_bits(full_bits)
        flow_counts[flow_bits] = get(flow_counts, flow_bits, 0) + reads
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

        repaired = get!(repair_cache, flow_bits) do
            repair_flow(bits_vector(flow_bits), qubo, components)
        end
        if isnothing(best_repaired) ||
            repaired.exact_repaired_qubo_energy < best_repaired.exact_repaired_qubo_energy
            best_repaired = repaired
            best_repaired_flow_bits = flow_bits
        end
    end

    hit_stats = hit_rate_stat_values(
        total_reads,
        top50_hits,
        top10_hits,
        global_hits,
        projected_pool_hits,
        elapsed_sec,
    )
    best_info = isempty(best_repaired_flow_bits) ? nothing : rank_info(best_repaired_flow_bits, top_flows)

    return (
        unique_full_states = length(aggregate_full),
        unique_flows = length(flow_counts),
        total_reads = total_reads,
        projected_pool_hits = projected_pool_hits,
        encoded_pool_hits = encoded_pool_hits,
        encoded_global_hits = encoded_global_hits,
        top50_hits = top50_hits,
        top10_hits = top10_hits,
        global_hits = global_hits,
        gurobi_pool_feasible_hits = projected_pool_hits,
        hit_stats = hit_stats,
        best_repaired = best_repaired,
        best_repaired_flow_bits = best_repaired_flow_bits,
        best_info = best_info,
    )
end

function direct_write_summary(
    path,
    config::DirectHardwareConfig,
    params::DirectQaoaParameters,
    paths,
    aggregate_full::Dict{String,Int},
    qubo::QuboData,
    components::Vector{AuxComponent},
    top_flows::Dict{String,NamedTuple};
    submitted_jobs::Int,
    elapsed_sec::Real,
)
    summary = direct_summarize_counts(
        aggregate_full,
        qubo,
        components,
        top_flows;
        elapsed_sec = elapsed_sec,
    )
    best_info = summary.best_info
    best_rank = isnothing(best_info) || ismissing(best_info.rank) ? "" : best_info.rank
    best_match = isnothing(best_info) ? "" : best_info.match
    best_energy = isnothing(summary.best_repaired) ? "" : summary.best_repaired.exact_repaired_qubo_energy
    best_full = isnothing(summary.best_repaired) ? "" : summary.best_repaired.repaired_full_bits

    open(path, "w") do io
        println(io, csv_row(direct_summary_header()))
        println(
            io,
            csv_row([
                DIRECT_HARDWARE_ALGORITHM,
                config.run_hardware ? "hardware" : "dry_run",
                config.backend_name,
                relpath(params.path, STUDY_ROOT),
                params.p,
                config.final_reads,
                config.repeats,
                join(config.transpile_seeds, ';'),
                submitted_jobs,
                @sprintf("%.6f", elapsed_sec),
                summary.unique_full_states,
                summary.unique_flows,
                summary.total_reads,
                summary.projected_pool_hits,
                summary.encoded_pool_hits,
                summary.encoded_global_hits,
                summary.top50_hits,
                summary.top10_hits,
                summary.global_hits,
                summary.gurobi_pool_feasible_hits,
                summary.hit_stats...,
                best_rank,
                best_energy,
                best_match,
                summary.best_repaired_flow_bits,
                best_full,
                relpath(paths["job_manifest"], config.output_dir),
                relpath(paths["backend_metadata"], config.output_dir),
                relpath(paths["raw_counts"], config.output_dir),
                relpath(paths["scored_counts"], config.output_dir),
            ]),
        )
    end
end

function direct_sampler_counts(result)
    pub_result = result[0]
    data = pub_result.data
    if pyhasattr(data, "c")
        return direct_py_dict_to_julia(data.c.get_counts())
    elseif pyhasattr(data, "meas")
        return direct_py_dict_to_julia(data.meas.get_counts())
    elseif pyhasattr(pub_result, "join_data")
        return direct_py_dict_to_julia(pub_result.join_data().get_counts())
    end
    error("Could not find sampler count data in Runtime SamplerV2 result")
end

function direct_write_count_rows!(
    raw_io,
    scored_io,
    aggregate_full::Dict{String,Int},
    config::DirectHardwareConfig,
    qubo::QuboData,
    components::Vector{AuxComponent},
    top_flows::Dict{String,NamedTuple},
    job,
    counts,
)
    total_reads = sum(values(counts); init = 0)

    for (qiskit_bitstring, reads) in sort(collect(counts); by = first)
        full_bits = qiskit_key_to_full_bits(qiskit_bitstring, qubo.n)
        flow_bits = flow_bits_from_full_bits(full_bits)
        aggregate_full[full_bits] = get(aggregate_full, full_bits, 0) + reads
        probability = total_reads == 0 ? 0.0 : reads / total_reads
        repaired = repair_flow(bits_vector(flow_bits), qubo, components)
        info = rank_info(flow_bits, top_flows)

        println(
            raw_io,
            csv_row([
                DIRECT_HARDWARE_ALGORITHM,
                "hardware",
                config.backend_name,
                job["repeat"],
                job["transpile_seed"],
                job["job_id"],
                qiskit_bitstring,
                full_bits,
                flow_bits,
                reads,
            ]),
        )

        println(
            scored_io,
            csv_row([
                DIRECT_HARDWARE_ALGORITHM,
                "hardware",
                config.backend_name,
                job["repeat"],
                job["transpile_seed"],
                job["job_id"],
                full_bits,
                flow_bits,
                reads,
                probability,
                direct_full_qubo_energy(full_bits, qubo),
                repaired.exact_repaired_qubo_energy,
                info.rank,
                info.match,
                info.matched_ip_obj,
                repaired.repaired_aux_bits,
                repaired.repaired_full_bits,
                full_bits == repaired.repaired_full_bits,
            ]),
        )
    end
end

function direct_job_status_text(job)
    if pyhasattr(job, "status")
        try
            return string(job.status())
        catch
            return "unknown"
        end
    end
    return "unknown"
end

function direct_record_job_failure!(job_plan, runtime_job, submitted::Bool, err)
    if runtime_job !== nothing
        job_plan["status"] = direct_job_status_text(runtime_job)
    elseif submitted
        job_plan["status"] = get(job_plan, "status", "submitted")
    else
        job_plan["status"] = "failed_before_submission"
    end
    job_plan["error"] = submitted ?
        "Runtime job result retrieval or scoring failed after submission." :
        "Runtime job transpile or submission failed before a job ID was assigned."
    job_plan["error_type"] = string(typeof(err))
    job_plan["error_recorded_at_utc"] = direct_utc_timestamp()
    return nothing
end

function direct_runtime_service(runtime, config::DirectHardwareConfig)
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

function direct_runtime_service_or_error(runtime, config::DirectHardwareConfig, paths)
    original_error_type = nothing
    try
        return direct_runtime_service(runtime, config)
    catch err
        original_error_type = typeof(err)
    end
    write_json_file(
        paths["backend_metadata"],
        Dict{String,Any}(
            "schema_version" => DIRECT_HARDWARE_SCHEMA_VERSION,
            "created_at_utc" => direct_utc_timestamp(),
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

function direct_run_hardware_jobs!(
    config::DirectHardwareConfig,
    qubo::QuboData,
    components::Vector{AuxComponent},
    top_flows::Dict{String,NamedTuple},
    params::DirectQaoaParameters,
    circuit,
    circuit_info,
    jobs,
    paths,
)
    qiskit = pyimport("qiskit")
    runtime = pyimport("qiskit_ibm_runtime")
    service = direct_runtime_service_or_error(runtime, config, paths)
    backend = service.backend(config.backend_name)
    write_json_file(paths["backend_metadata"], direct_real_backend_metadata(config, backend))

    aggregate_full = Dict{String,Int}()
    submitted_jobs = 0

    open(paths["raw_counts"], "w") do raw_io
        println(raw_io, csv_row([
            "algorithm", "mode", "backend", "repeat", "transpile_seed", "job_id",
            "qiskit_bitstring", "full_bits", "flow_bits", "reads",
        ]))
        open(paths["scored_counts"], "w") do scored_io
            println(scored_io, csv_row([
                "algorithm", "mode", "backend", "repeat", "transpile_seed", "job_id",
                "full_bits", "flow_bits", "reads", "probability", "raw_qubo_energy",
                "exact_repaired_qubo_energy", "top50_rank", "match", "matched_ip_obj",
                "repaired_aux_bits", "repaired_full_bits", "encoded_repaired_match",
            ]))

            for job_plan in jobs
                runtime_job = nothing
                submitted = false
                try
                    @info "Submitting direct full-QUBO fixed-parameter QAOA sampler job" repeat=job_plan["repeat"] transpile_seed=job_plan["transpile_seed"] shots=config.final_reads parameter_artifact=params.path
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
                    job_plan["status"] = direct_job_status_text(runtime_job)
                    submitted_jobs += 1
                    direct_persist_job_manifest!(config, qubo, params, circuit_info, jobs, paths)

                    result = runtime_job.result()
                    counts = direct_sampler_counts(result)
                    job_plan["status"] = direct_job_status_text(runtime_job)
                    job_plan["unique_raw_bitstrings"] = length(counts)
                    job_plan["total_reads"] = sum(values(counts); init = 0)
                    direct_persist_job_manifest!(config, qubo, params, circuit_info, jobs, paths)
                    direct_write_count_rows!(
                        raw_io,
                        scored_io,
                        aggregate_full,
                        config,
                        qubo,
                        components,
                        top_flows,
                        job_plan,
                        counts,
                    )
                    flush(raw_io)
                    flush(scored_io)
                catch err
                    direct_record_job_failure!(job_plan, runtime_job, submitted || get(job_plan, "submitted", false), err)
                    direct_persist_job_manifest!(config, qubo, params, circuit_info, jobs, paths)
                    rethrow()
                end
            end
        end
    end

    return aggregate_full, submitted_jobs
end

function main()
    config = read_direct_hardware_config()
    mkpath(config.output_dir)

    qubo = read_qubo_data(STUDY_ROOT)
    components = read_aux_components(joinpath(REDUCED_DIR, "auxiliary_components.csv"), qubo)
    top_flows = read_exact_top_flows(joinpath(REDUCED_DIR, "reduced_exact_top_flows.csv"))
    params = read_direct_qaoa_parameters(config.parameter_path, qubo)
    circuit, zz_pairs = build_direct_full_qubo_qaoa_circuit(qubo, params.p, params.angles)
    circuit_info = direct_circuit_metadata(circuit, zz_pairs)
    paths = direct_output_paths(config)
    jobs = direct_planned_jobs(config)
    aggregate_full = Dict{String,Int}()
    submitted_jobs = 0
    start_time = time()

    direct_persist_job_manifest!(config, qubo, params, circuit_info, jobs, paths)

    if config.run_hardware
        aggregate_full, submitted_jobs = direct_run_hardware_jobs!(
            config,
            qubo,
            components,
            top_flows,
            params,
            circuit,
            circuit_info,
            jobs,
            paths,
        )
    else
        write_json_file(paths["backend_metadata"], direct_dry_run_backend_metadata(config))
        direct_write_empty_count_files(paths)
    end

    elapsed_sec = time() - start_time
    direct_write_summary(
        paths["summary"],
        config,
        params,
        paths,
        aggregate_full,
        qubo,
        components,
        top_flows;
        submitted_jobs = submitted_jobs,
        elapsed_sec = elapsed_sec,
    )
    direct_persist_job_manifest!(config, qubo, params, circuit_info, jobs, paths)

    if config.run_hardware
        println("Direct full-QUBO QAOA hardware run complete.")
    else
        println("Direct full-QUBO QAOA hardware dry run complete. No hardware jobs were submitted.")
    end
    println("Output directory: ", config.output_dir)
    println("Parameter artifact: ", params.path)
    println("Manifest: ", paths["job_manifest"])
    println("Summary: ", paths["summary"])
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
