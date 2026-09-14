#!/usr/bin/env julia
# Run in the JuliQAOA project documented for find_reduced_qaoa_angles_juliqaoa.jl.
const ENTRY_NS = time_ns()
module RevisionQAOA
using Random, Statistics, LinearAlgebra, TOML, Dates, SHA
module Angles
include(joinpath(@__DIR__, "find_reduced_qaoa_angles_juliqaoa.jl"))
end
const A = Angles

function main()
    BLAS.set_num_threads(1)
    Threads.nthreads() == 1 || error("Run with one Julia thread")
    output = joinpath(A.STUDENT_ROOT, "ds_mfg_revision_audit")
    mkpath(output)
    environment_seconds = (time_ns()-Main.ENTRY_NS)/1e9
    start = time_ns()
    reduced = joinpath(A.STUDENT_ROOT, "ds_mfg_reduced_flow_objective")
    scalars = A.read_reduced_scalars(joinpath(reduced, "reduced_scalars.csv"))
    L = A.read_reduced_vector(joinpath(reduced, "reduced_L_vector.csv"))
    Q = A.read_reduced_matrix(joinpath(reduced, "reduced_Q_matrix.csv"))
    energy, bitstrings = A.enumerate_surrogate_energies(L, Q; scale=scalars.scale, offset=scalars.offset)
    top = A.read_exact_top_flows(joinpath(reduced, "reduced_exact_top_flows.csv"))
    normalized = (energy .- mean(energy)) ./ std(energy)
    mixer = A.mixer_x(scalars.n)
    preprocessing_seconds = (time_ns()-start)/1e9
    # Compile the optimizer on an unrelated three-qubit cost before timed runs.
    start = time_ns()
    Random.seed!(1)
    A.find_angles_bh(1, A.mixer_x(3), collect(0.:7.); max=false, niter=1, verbose=false)
    warmup_seconds = (time_ns()-start)/1e9
    rows = NamedTuple[]
    runs = NamedTuple[]
    for (seed, niter) in ((91001,1), (91002,1), (91003,1), (91001,5))
        Random.seed!(seed)
        start = time_ns()
        sets, expectations = A.find_angles_bh(5, mixer, normalized;
            max=false, niter=niter, verbose=true)
        optimization_seconds = (time_ns()-start)/1e9
        start = time_ns()
        for p in 1:5
            probabilities = A.JuliQAOA.probabilities(sets[p], mixer, normalized)
            s = A.probability_summary(probabilities, energy, top, bitstrings)
            push!(rows, (seed=seed, basinhopping_iterations=niter, p=p,
                global_probability=s.global_probability, top50_probability=s.top50_probability,
                expected_surrogate_energy=sum(probabilities .* energy),
                normalized_angles_beta_then_gamma=A.angles_string(sets[p])))
        end
        scoring_seconds = (time_ns()-start)/1e9
        push!(runs, (seed=seed, basinhopping_iterations=niter, maximum_p=5,
            optimization_seconds=optimization_seconds, probability_scoring_seconds=scoring_seconds))
        for (filename, data) in (("qaoa_sensitivity.csv",rows),("qaoa_timing.csv",runs))
            open(joinpath(output,filename), "w") do io
                println(io, join(string.(propertynames(first(data))), ','))
                for row in data
                    println(io, join(string.(values(row)), ','))
                end
            end
        end
        println(last(runs))
        flush(stdout)
    end
    metadata = Dict("julia_version"=>string(VERSION), "julia_threads"=>Threads.nthreads(),
        "cpu"=>Sys.cpu_info()[1].model, "created_utc"=>string(now(UTC)),
        "environment_load_seconds"=>environment_seconds,
        "surrogate_enumeration_and_setup_seconds"=>preprocessing_seconds,
        "optimizer_warmup_seconds"=>warmup_seconds,
        "wall_seconds_since_script_entry"=>(time_ns()-Main.ENTRY_NS)/1e9,
        "script_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "angle_helper_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"find_reduced_qaoa_angles_juliqaoa.jl")))),
        "juliqaoa_manifest_sha256"=>bytes2hex(sha256(read(joinpath(dirname(Base.active_project()),"Manifest.toml")))),
        "timing_policy"=>"new warm CPU rerun; each optimization timer includes the complete p=1..5 ladder; not a reconstruction of historical times")
    open(joinpath(output,"qaoa_metadata.toml"), "w") do io
        TOML.print(io, metadata)
    end
end
end
RevisionQAOA.main()
