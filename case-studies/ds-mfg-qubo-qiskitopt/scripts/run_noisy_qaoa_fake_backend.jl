#!/usr/bin/env julia

import Pkg

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    Pkg.activate(joinpath(@__DIR__, ".."))
    Pkg.instantiate()
end

using Printf
using PythonCall

include(joinpath(@__DIR__, "run_ibm_qaoa_pilot.jl"))

const NOISY_ALGORITHM = "QAOA_reduced_surrogate_JuliQAOA_FakeFez_noisy_simulation"
const NOISY_MODE = "model_based_noisy_simulation"
const NOISY_SCHEMA_VERSION = 1

struct NoisyConfig
    fake_backend_class::String
    final_reads::Int
    repeats::Int
    transpile_seeds::Vector{Int}
    simulator_seed_start::Int
    optimization_level::Int
    output_dir::String
    run_simulation::Bool
    max_jobs::Int
end

function parse_nonnegative_int_env(name::AbstractString, default::Integer)
    raw = strip(get(ENV, name, string(default)))
    value = parse(Int, raw)
    value >= 0 || error("$(name) must be a nonnegative integer, got $(raw)")
    return value
end

function read_noisy_config()
    output_name = strip(get(ENV, "DSMFG_NOISY_OUTPUT_DIR", "ds_mfg_fake_fez_qaoa_noisy_4096x3"))
    isempty(output_name) && error("DSMFG_NOISY_OUTPUT_DIR must not be empty")
    output_dir = isabspath(output_name) ? output_name : joinpath(STUDY_ROOT, output_name)

    backend_class = strip(get(ENV, "DSMFG_NOISY_FAKE_BACKEND", "FakeFez"))
    isempty(backend_class) && error("DSMFG_NOISY_FAKE_BACKEND must not be empty")

    optimization_level = parse_nonnegative_int_env("DSMFG_NOISY_OPTIMIZATION_LEVEL", 3)
    0 <= optimization_level <= 3 ||
        error("DSMFG_NOISY_OPTIMIZATION_LEVEL must be between 0 and 3, got $(optimization_level)")

    return NoisyConfig(
        backend_class,
        parse_positive_int_env("DSMFG_NOISY_FINAL_READS", 4096),
        parse_positive_int_env("DSMFG_NOISY_REPEATS", 3),
        parse_int_list(get(ENV, "DSMFG_NOISY_TRANSPILE_SEEDS", "92001:92003")),
        parse_positive_int_env("DSMFG_NOISY_SIMULATOR_SEED", 93001),
        optimization_level,
        normpath(abspath(output_dir)),
        parse_bool_env("DSMFG_RUN_NOISY_SIMULATION"),
        parse_nonnegative_int_env("DSMFG_NOISY_MAX_JOBS", 0),
    )
end

function planned_noisy_jobs(config::NoisyConfig)
    jobs = Vector{Dict{String,Any}}()
    job_index = 1
    for repeat in 1:config.repeats, seed in config.transpile_seeds
        push!(
            jobs,
            Dict{String,Any}(
                "job_index" => job_index,
                "repeat" => repeat,
                "transpile_seed" => seed,
                "simulator_seed" => config.simulator_seed_start + job_index - 1,
                "shots" => config.final_reads,
                "status" => "planned",
            ),
        )
        job_index += 1
    end
    if config.max_jobs > 0
        return jobs[1:min(config.max_jobs, length(jobs))]
    end
    return jobs
end

function noisy_output_paths(config::NoisyConfig)
    return Dict{String,String}(
        "manifest" => joinpath(config.output_dir, "job_manifest.json"),
        "backend_metadata" => joinpath(config.output_dir, "backend_metadata.json"),
        "raw_counts" => joinpath(config.output_dir, "raw_counts.csv"),
        "scored_counts" => joinpath(config.output_dir, "scored_counts.csv"),
        "summary" => joinpath(config.output_dir, "summary.csv"),
    )
end

function fake_backend_instance(config::NoisyConfig)
    fake_provider = pyimport("qiskit_ibm_runtime.fake_provider")
    pyhasattr(fake_provider, config.fake_backend_class) ||
        error("qiskit_ibm_runtime.fake_provider has no $(config.fake_backend_class)")
    return pygetattr(fake_provider, config.fake_backend_class)()
end

function noisy_backend_metadata(config::NoisyConfig, backend)
    metadata = Dict{String,Any}(
        "schema_version" => NOISY_SCHEMA_VERSION,
        "created_at_utc" => utc_timestamp(),
        "mode" => config.run_simulation ? NOISY_MODE : "dry_run",
        "fake_backend_class" => config.fake_backend_class,
        "model_source" => "qiskit_ibm_runtime.fake_provider",
        "simulator_source" => "qiskit_aer.AerSimulator.from_backend",
        "calibration_scope" => "packaged fake backend model; not calibrated to the IBM hardware pilot execution window",
    )

    for attr in ("name", "num_qubits", "backend_version")
        value = maybe_py_attr(backend, attr)
        isnothing(value) && continue
        if attr == "num_qubits"
            metadata[attr] = pyconvert(Int, value)
        else
            metadata[attr] = pyconvert(String, value)
        end
    end

    target = maybe_py_attr(backend, "target")
    if !isnothing(target) && pyhasattr(target, "operation_names")
        metadata["operation_names"] = sort(pyconvert(Vector{String}, target.operation_names))
    end

    if pyhasattr(backend, "configuration")
        try
            configuration = backend.configuration()
            pyhasattr(configuration, "basis_gates") &&
                (metadata["basis_gates"] = sort(pyconvert(Vector{String}, configuration.basis_gates)))
            pyhasattr(configuration, "simulator") &&
                (metadata["backend_reports_simulator"] = pyconvert(Bool, configuration.simulator))
        catch
            metadata["configuration_warning"] = "backend.configuration() metadata unavailable"
        end
    end

    return metadata
end

function noisy_manifest(config::NoisyConfig, data, circuit_info, jobs, paths)
    return Dict{String,Any}(
        "schema_version" => NOISY_SCHEMA_VERSION,
        "created_at_utc" => utc_timestamp(),
        "script" => relpath(@__FILE__, STUDY_ROOT),
        "mode" => config.run_simulation ? NOISY_MODE : "dry_run",
        "configuration" => Dict{String,Any}(
            "fake_backend_class" => config.fake_backend_class,
            "final_reads" => config.final_reads,
            "repeats" => config.repeats,
            "transpile_seeds" => config.transpile_seeds,
            "simulator_seed_start" => config.simulator_seed_start,
            "optimization_level" => config.optimization_level,
            "output_dir" => public_artifact_path(config.output_dir, STUDY_ROOT),
            "run_simulation" => config.run_simulation,
            "max_jobs" => config.max_jobs,
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
        "interpretation_warning" => "This is a model-based fake-backend simulation. It is not a calibrated prediction of the IBM hardware pilot.",
    )
end

function persist_noisy_manifest!(config::NoisyConfig, data, circuit_info, jobs, paths)
    write_json_file(paths["manifest"], noisy_manifest(config, data, circuit_info, jobs, paths))
    return nothing
end

function write_noisy_empty_count_files(paths)
    open(paths["raw_counts"], "w") do io
        println(io, csv_row([
            "algorithm", "mode", "fake_backend", "repeat", "transpile_seed", "simulator_seed",
            "qiskit_bitstring", "flow_bits", "reads",
        ]))
    end
    open(paths["scored_counts"], "w") do io
        println(io, csv_row([
            "algorithm", "mode", "fake_backend", "repeat", "transpile_seed", "simulator_seed",
            "flow_bits", "reads", "probability", "surrogate_qubo_energy",
            "top50_rank", "exact_repaired_qubo_energy", "match", "matched_ip_obj",
            "repaired_aux_bits", "repaired_full_bits",
        ]))
    end
end

function write_noisy_count_rows!(raw_io, scored_io, aggregate::Dict{String,Int}, config::NoisyConfig, data, job, counts)
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
                NOISY_ALGORITHM,
                NOISY_MODE,
                config.fake_backend_class,
                job["repeat"],
                job["transpile_seed"],
                job["simulator_seed"],
                qiskit_bitstring,
                flow_bits,
                reads,
            ]),
        )

        println(
            scored_io,
            csv_row([
                NOISY_ALGORITHM,
                NOISY_MODE,
                config.fake_backend_class,
                job["repeat"],
                job["transpile_seed"],
                job["simulator_seed"],
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

function noisy_summary_header()
    return [
        "algorithm", "mode", "fake_backend", "angle_target", "p", "final_reads",
        "repeats", "transpile_seeds", "simulator_seed_start", "optimization_level",
        "completed_jobs", "elapsed_sec", "unique_flows", "total_reads", "top50_hits",
        "top10_hits", "global_hits", "gurobi_pool_feasible_hits", hit_rate_stat_headers()...,
        "best_top50_rank", "best_top50_exact_repaired_qubo_energy", "best_top50_match",
        "best_top50_flow_bits", "transpiled_depth_min", "transpiled_depth_max",
        "transpiled_cz_min", "transpiled_cz_max", "manifest_path", "backend_metadata_path",
        "raw_counts_path", "scored_counts_path",
    ]
end

function completed_job_values(jobs, key::AbstractString)
    values = Int[]
    for job in jobs
        get(job, "status", "") == "DONE" || continue
        haskey(job, key) && push!(values, Int(job[key]))
    end
    return values
end

function min_or_blank(values)
    isempty(values) ? "" : minimum(values)
end

function max_or_blank(values)
    isempty(values) ? "" : maximum(values)
end

function write_noisy_summary(path, config::NoisyConfig, paths, jobs, aggregate::Dict{String,Int}; elapsed_sec::Real, data)
    summary = summarize_counts(aggregate, data; elapsed_sec = elapsed_sec)
    best_rank = isnothing(summary.best_hit) ? "" : summary.best_hit.rank
    best_energy = isnothing(summary.best_hit) ? "" : summary.best_hit.exact_repaired_qubo_energy
    best_match = isnothing(summary.best_hit) ? "" : summary.best_hit.match
    best_flow = isnothing(summary.best_hit) ? "" : summary.best_hit.flow_bits
    depth_values = completed_job_values(jobs, "transpiled_depth")
    cz_values = completed_job_values(jobs, "transpiled_cz_count")
    completed_jobs = count(job -> get(job, "status", "") == "DONE", jobs)

    open(path, "w") do io
        println(io, csv_row(noisy_summary_header()))
        println(
            io,
            csv_row([
                NOISY_ALGORITHM,
                config.run_simulation ? NOISY_MODE : "dry_run",
                config.fake_backend_class,
                ANGLE_TARGET,
                ANGLE_P,
                config.final_reads,
                config.repeats,
                join(config.transpile_seeds, ';'),
                config.simulator_seed_start,
                config.optimization_level,
                completed_jobs,
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
                min_or_blank(depth_values),
                max_or_blank(depth_values),
                min_or_blank(cz_values),
                max_or_blank(cz_values),
                relpath(paths["manifest"], config.output_dir),
                relpath(paths["backend_metadata"], config.output_dir),
                relpath(paths["raw_counts"], config.output_dir),
                relpath(paths["scored_counts"], config.output_dir),
            ]),
        )
    end
end

function run_noisy_jobs!(config::NoisyConfig, data, circuit, circuit_info, jobs, paths)
    qiskit = pyimport("qiskit")
    aer = pyimport("qiskit_aer")
    backend = fake_backend_instance(config)
    write_json_file(paths["backend_metadata"], noisy_backend_metadata(config, backend))
    simulator = aer.AerSimulator.from_backend(backend)

    aggregate = Dict{String,Int}()
    open(paths["raw_counts"], "w") do raw_io
        println(raw_io, csv_row([
            "algorithm", "mode", "fake_backend", "repeat", "transpile_seed", "simulator_seed",
            "qiskit_bitstring", "flow_bits", "reads",
        ]))
        open(paths["scored_counts"], "w") do scored_io
            println(scored_io, csv_row([
                "algorithm", "mode", "fake_backend", "repeat", "transpile_seed", "simulator_seed",
                "flow_bits", "reads", "probability", "surrogate_qubo_energy",
                "top50_rank", "exact_repaired_qubo_energy", "match", "matched_ip_obj",
                "repaired_aux_bits", "repaired_full_bits",
            ]))

            for job in jobs
                try
                    @info "Running local fake-backend noisy QAOA simulation" repeat=job["repeat"] transpile_seed=job["transpile_seed"] shots=config.final_reads simulator_seed=job["simulator_seed"]
                    job["status"] = "transpiling"
                    persist_noisy_manifest!(config, data, circuit_info, jobs, paths)
                    transpiled = qiskit.transpile(
                        circuit;
                        backend = backend,
                        seed_transpiler = job["transpile_seed"],
                        optimization_level = config.optimization_level,
                    )
                    ops = py_dict_to_julia(transpiled.count_ops())
                    job["transpiled_depth"] = pyconvert(Int, transpiled.depth())
                    job["transpiled_cz_count"] = get(ops, "cz", 0)
                    job["status"] = "simulating"
                    persist_noisy_manifest!(config, data, circuit_info, jobs, paths)

                    result = simulator.run(
                        transpiled;
                        shots = config.final_reads,
                        seed_simulator = job["simulator_seed"],
                    ).result()
                    counts = py_dict_to_julia(result.get_counts())
                    job["unique_raw_bitstrings"] = length(counts)
                    job["total_reads"] = sum(values(counts); init = 0)
                    write_noisy_count_rows!(raw_io, scored_io, aggregate, config, data, job, counts)
                    flush(raw_io)
                    flush(scored_io)
                    job["status"] = "DONE"
                    persist_noisy_manifest!(config, data, circuit_info, jobs, paths)
                catch err
                    job["status"] = "FAILED"
                    job["error"] = "Local fake-backend simulation failed."
                    job["error_type"] = string(typeof(err))
                    job["error_recorded_at_utc"] = utc_timestamp()
                    persist_noisy_manifest!(config, data, circuit_info, jobs, paths)
                    rethrow()
                end
            end
        end
    end
    return aggregate
end

function main()
    config = read_noisy_config()
    data = load_problem_data()
    mkpath(config.output_dir)

    circuit = build_qaoa_circuit(data)
    circuit_info = circuit_metadata(circuit)
    paths = noisy_output_paths(config)
    jobs = planned_noisy_jobs(config)
    start_time = time()
    aggregate = Dict{String,Int}()

    persist_noisy_manifest!(config, data, circuit_info, jobs, paths)
    backend = fake_backend_instance(config)
    write_json_file(paths["backend_metadata"], noisy_backend_metadata(config, backend))

    if config.run_simulation
        aggregate = run_noisy_jobs!(config, data, circuit, circuit_info, jobs, paths)
    else
        write_noisy_empty_count_files(paths)
    end

    elapsed_sec = time() - start_time
    write_noisy_summary(paths["summary"], config, paths, jobs, aggregate; elapsed_sec = elapsed_sec, data = data)
    persist_noisy_manifest!(config, data, circuit_info, jobs, paths)

    if config.run_simulation
        println("Fake-backend noisy QAOA simulation complete.")
    else
        println("Fake-backend noisy QAOA dry run complete. Set DSMFG_RUN_NOISY_SIMULATION=true to run overnight.")
    end
    println("Output directory: ", config.output_dir)
    println("Manifest: ", paths["manifest"])
    println("Summary: ", paths["summary"])
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
