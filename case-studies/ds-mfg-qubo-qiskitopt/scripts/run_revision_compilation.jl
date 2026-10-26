#!/usr/bin/env julia
# Offline timing proxy using the actual pilot circuit builder and FakeFez.
# No Runtime service is constructed and no jobs are submitted.
const ENTRY_NS = time_ns()
module Pilot
include(joinpath(@__DIR__, "run_ibm_qaoa_pilot.jl"))
end
using TOML, Dates, SHA

function main()
    output = joinpath(Pilot.STUDY_ROOT,"ds_mfg_revision_audit")
    mkpath(output)
    environment_seconds = (time_ns()-ENTRY_NS)/1e9
    start = time_ns()
    data = Pilot.load_problem_data()
    qiskit = Pilot.pyimport("qiskit")
    fake = Pilot.pyimport("qiskit_ibm_runtime.fake_provider")
    backend = fake.FakeFez()
    setup_seconds = (time_ns()-start)/1e9
    start = time_ns()
    circuit = Pilot.build_qaoa_circuit(data)
    build_seconds = (time_ns()-start)/1e9
    metadata = Pilot.circuit_metadata(circuit)
    @assert metadata["depth"] == 86
    @assert metadata["operation_counts"] == Dict("h"=>19,"rx"=>95,"rz"=>95,"rzz"=>255,"measure"=>19)
    start = time_ns()
    qiskit.transpile(circuit; backend=backend,seed_transpiler=91999,
        optimization_level=Pilot.TRANSPILER_OPTIMIZATION_LEVEL)
    warmup_seconds = (time_ns()-start)/1e9
    open(joinpath(output,"compilation_proxy.csv"),"w") do io
        println(io,"backend_model,seed,optimization_level,compilation_seconds,depth,two_qubit_gates")
        for seed in 92001:92003
            start = time_ns()
            compiled = qiskit.transpile(circuit; backend=backend,seed_transpiler=seed,
                optimization_level=Pilot.TRANSPILER_OPTIMIZATION_LEVEL)
            elapsed = (time_ns()-start)/1e9
            summary = Pilot.circuit_summary(compiled)
            println(io,Pilot.csv_row(["FakeFez",seed,Pilot.TRANSPILER_OPTIMIZATION_LEVEL,
                elapsed,summary["depth"],summary["two_qubit_gate_count"]]))
        end
    end
    summary = Dict("created_utc"=>string(now(UTC)),"julia_version"=>string(VERSION),
        "qiskit_version"=>Pilot.pyconvert(String,qiskit.__version__),
        "runtime_version"=>Pilot.pyconvert(String,Pilot.pyimport("qiskit_ibm_runtime").__version__),
        "environment_load_seconds"=>environment_seconds,"load_and_fake_backend_seconds"=>setup_seconds,
        "logical_circuit_construction_seconds"=>build_seconds,
        "transpilation_warmup_seconds"=>warmup_seconds,
        "script_sha256"=>bytes2hex(sha256(read(@__FILE__))),
        "pilot_builder_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"run_ibm_qaoa_pilot.jl")))),
        "scope"=>"new offline FakeFez proxy; not historical ibm_fez compilation or queue time")
    open(joinpath(output,"compilation_metadata.toml"),"w") do io
        TOML.print(io,summary)
    end
end
main()
