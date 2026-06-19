#!/usr/bin/env julia

using Printf

const STUDY_ROOT = abspath(joinpath(@__DIR__, ".."))
const OUTPUT_DIR = joinpath(STUDY_ROOT, "ds_mfg_simulator_hardware_comparison")
const IDEAL_REL = "ds_mfg_qaoa_juliqaoa_objective_transfer_highread/qaoa_juliqaoa_transfer_summary.csv"
const NOISY_REL = "ds_mfg_fake_fez_qaoa_noisy_4096x3/summary.csv"
const NOISY_SCRIPT_REL = "scripts/run_noisy_qaoa_fake_backend.jl"
const HARDWARE_REL = "ds_mfg_ibm_qaoa_pilot_fez_4096x3/summary.csv"
const COMPARISON_REL = "ds_mfg_simulator_hardware_comparison/simulator_hardware_comparison_summary.csv"
const DEGRADATION_REL = "ds_mfg_simulator_hardware_comparison/simulator_hardware_degradation_summary.csv"

const EVENTS = (
    ("top50", "top-50 repaired-flow"),
    ("top10", "top-10 repaired-flow"),
    ("global", "global optimum"),
    ("gurobi_pool_feasible", "Gurobi-pool feasible"),
)

function read_csv_rows(relpath::AbstractString)
    path = joinpath(STUDY_ROOT, relpath)
    lines = filter(line -> !isempty(strip(line)), readlines(path))
    length(lines) >= 2 || error("Expected header and at least one row in $(relpath)")
    header = split(lines[1], ","; keepempty = true)
    rows = Dict{String,String}[]
    for (line_number, line) in enumerate(lines[2:end])
        fields = split(line, ","; keepempty = true)
        length(fields) == length(header) ||
            error("Row $(line_number + 1) in $(relpath) has $(length(fields)) fields; expected $(length(header))")
        push!(rows, Dict(header[index] => fields[index] for index in eachindex(header)))
    end
    return rows
end

function maybe_read_csv_rows(relpath::AbstractString)
    path = joinpath(STUDY_ROOT, relpath)
    isfile(path) || return nothing
    return read_csv_rows(relpath)
end

function csv_value(value)
    value === nothing && return ""
    text = string(value)
    if occursin(',', text) || occursin('"', text) || occursin('\n', text) || occursin('\r', text)
        return "\"" * replace(text, "\"" => "\"\"") * "\""
    end
    return text
end

csv_row(values) = join(csv_value.(values), ',')

function write_csv(path::AbstractString, header, rows)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, csv_row(header))
        for row in rows
            println(io, csv_row(row))
        end
    end
end

function require_row(rows, predicate, label::AbstractString)
    matches = filter(predicate, rows)
    length(matches) == 1 || error("Expected exactly one $(label) row, got $(length(matches))")
    return only(matches)
end

function required_int(row, key::AbstractString)
    haskey(row, key) || error("Missing required column $(key)")
    return parse(Int, row[key])
end

function source_hit_rate(row, event::AbstractString)
    key = "$(event)_hit_rate"
    haskey(row, key) && !isempty(row[key]) && return row[key]
    total = required_int(row, "total_reads")
    hits = required_int(row, "$(event)_hits")
    return format_number(hits / total)
end

function format_number(value)
    x = Float64(value)
    isnan(x) && return ""
    isinf(x) && return x > 0.0 ? "Inf" : "-Inf"
    return @sprintf("%.12g", x)
end

function comparison_row(;
    evidence_tier,
    row_status,
    label,
    backend_or_model,
    source_artifact,
    row = nothing,
    notes = "",
)
    if row === nothing
        return [
            evidence_tier,
            row_status,
            label,
            backend_or_model,
            source_artifact,
            "top10",
            "5",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            notes,
        ]
    end

    return [
        evidence_tier,
        row_status,
        label,
        backend_or_model,
        source_artifact,
        get(row, "angle_target", "top10"),
        get(row, "p", "5"),
        row["total_reads"],
        row["top50_hits"],
        source_hit_rate(row, "top50"),
        row["top10_hits"],
        source_hit_rate(row, "top10"),
        row["global_hits"],
        source_hit_rate(row, "global"),
        row["gurobi_pool_feasible_hits"],
        source_hit_rate(row, "gurobi_pool_feasible"),
        row["best_top50_rank"],
        row["best_top50_exact_repaired_qubo_energy"],
        row["best_top50_match"],
        row["best_top50_flow_bits"],
        notes,
    ]
end

function event_stats(row, event::AbstractString)
    hits = required_int(row, "$(event)_hits")
    total = required_int(row, "total_reads")
    return hits, total, hits / Float64(total)
end

function optional_event_stats(row, event::AbstractString)
    row === nothing && return nothing
    required_int(row, "total_reads") > 0 || return nothing
    return event_stats(row, event)
end

function degradation_rows(ideal, noisy, hardware)
    rows = Vector{Vector{String}}()
    for (event, label) in EVENTS
        ideal_hits, ideal_total, ideal_rate = event_stats(ideal, event)
        noisy_stats = optional_event_stats(noisy, event)
        hardware_hits, hardware_total, hardware_rate = event_stats(hardware, event)
        noisy_hits = ""
        noisy_total = ""
        noisy_rate = ""
        noisy_ratio = ""
        if noisy_stats !== nothing
            noisy_hits_value, noisy_total_value, noisy_rate_value = noisy_stats
            noisy_hits = string(noisy_hits_value)
            noisy_total = string(noisy_total_value)
            noisy_rate = format_number(noisy_rate_value)
            noisy_ratio = format_number(noisy_rate_value / ideal_rate)
        end
        ratio = ideal_rate == 0.0 ? NaN : hardware_rate / ideal_rate
        absolute_drop = ideal_rate - hardware_rate
        relative_drop = ideal_rate == 0.0 ? NaN : absolute_drop / ideal_rate
        expected_hardware_budget_hits = ideal_rate * hardware_total
        hardware_minus_expected = hardware_hits - expected_hardware_budget_hits
        push!(
            rows,
            [
                event,
                label,
                string(ideal_hits),
                string(ideal_total),
                format_number(ideal_rate),
                noisy_hits,
                noisy_total,
                noisy_rate,
                noisy_ratio,
                string(hardware_hits),
                string(hardware_total),
                format_number(hardware_rate),
                format_number(absolute_drop),
                format_number(relative_drop),
                format_number(ratio),
                format_number(expected_hardware_budget_hits),
                format_number(hardware_minus_expected),
                noisy_stats === nothing ?
                    "No cached model-based noisy result is present yet; run scripts/run_noisy_qaoa_fake_backend.jl overnight to populate the noisy columns." :
                    "Rates use each row's own read budget; expected hardware-budget hits apply the ideal Aer rate to $(hardware_total) reads.",
            ],
        )
    end
    return rows
end

function main()
    ideal_rows = read_csv_rows(IDEAL_REL)
    ideal = require_row(
        ideal_rows,
        row -> row["angle_target"] == "top10" && row["p"] == "5",
        "top-10 targeted ideal Aer QAOA",
    )
    hardware = only(read_csv_rows(HARDWARE_REL))
    noisy_rows = maybe_read_csv_rows(NOISY_REL)
    noisy = nothing
    if noisy_rows !== nothing
        candidate = only(noisy_rows)
        if get(candidate, "mode", "") == "model_based_noisy_simulation" && required_int(candidate, "total_reads") > 0
            noisy = candidate
        end
    end

    pending_noisy_note = join(
        [
            "Pinned stack exposes qiskit_ibm_runtime.fake_provider.FakeFez and qiskit_aer.AerSimulator.from_backend.",
            "The overnight runner scripts/run_noisy_qaoa_fake_backend.jl uses hardware-comparable defaults but has not produced a cached result yet.",
            "Local timing probes of the depth-1557 and 1141-CZ transpiled circuit took 17.5 s for 64 shots and 116.6 s for 512 shots.",
            "Any produced result must be labeled model-based simulation and not treated as a calibrated hardware predictor.",
        ],
        " ",
    )

    comparison_header = [
        "evidence_tier",
        "row_status",
        "label",
        "backend_or_model",
        "source_artifact",
        "angle_target",
        "p",
        "total_reads",
        "top50_hits",
        "top50_hit_rate",
        "top10_hits",
        "top10_hit_rate",
        "global_hits",
        "global_hit_rate",
        "gurobi_pool_feasible_hits",
        "gurobi_pool_feasible_hit_rate",
        "best_top50_rank",
        "best_top50_exact_repaired_qubo_energy",
        "best_top50_match",
        "best_top50_flow_bits",
        "notes",
    ]
    comparison_rows = [
        comparison_row(
            evidence_tier = "ideal_simulator",
            row_status = "cached_result",
            label = "Ideal Aer transferred-angle QAOA",
            backend_or_model = "Qiskit Aer MPS local simulator",
            source_artifact = IDEAL_REL,
            row = ideal,
            notes = "Top-10-targeted p=5 JuliQAOA angles transferred into QiskitOpt.QAOA with MaximumIterations() = 0.",
        ),
        noisy === nothing ?
        comparison_row(
            evidence_tier = "model_based_noisy_simulator",
            row_status = "script_available_not_run",
            label = "FakeFez/Aer noise-model simulation",
            backend_or_model = "qiskit_ibm_runtime.fake_provider.FakeFez + qiskit_aer.AerSimulator.from_backend",
            source_artifact = NOISY_SCRIPT_REL,
            notes = pending_noisy_note,
        ) :
        comparison_row(
            evidence_tier = "model_based_noisy_simulator",
            row_status = "cached_result",
            label = "FakeFez/Aer noise-model simulation",
            backend_or_model = get(noisy, "fake_backend", "FakeFez"),
            source_artifact = NOISY_REL,
            row = noisy,
            notes = "Model-based fake-backend simulation; not a calibrated hardware predictor.",
        ),
        comparison_row(
            evidence_tier = "hardware",
            row_status = "cached_result",
            label = "IBM hardware pilot",
            backend_or_model = "ibm_fez",
            source_artifact = HARDWARE_REL,
            row = hardware,
            notes = "Fixed-parameter hardware pilot; descriptive feasibility evidence only.",
        ),
    ]
    write_csv(joinpath(STUDY_ROOT, COMPARISON_REL), comparison_header, comparison_rows)

    degradation_header = [
        "event",
        "event_label",
        "ideal_aer_hits",
        "ideal_aer_total_reads",
        "ideal_aer_hit_rate",
        "noisy_model_hits",
        "noisy_model_total_reads",
        "noisy_model_hit_rate",
        "noisy_to_ideal_hit_rate_ratio",
        "hardware_hits",
        "hardware_total_reads",
        "hardware_hit_rate",
        "absolute_hit_rate_drop",
        "relative_hit_rate_drop",
        "hardware_to_ideal_hit_rate_ratio",
        "ideal_expected_hits_at_hardware_reads",
        "hardware_hits_minus_ideal_expected",
        "notes",
    ]
    write_csv(joinpath(STUDY_ROOT, DEGRADATION_REL), degradation_header, degradation_rows(ideal, noisy, hardware))

    println("Wrote $(COMPARISON_REL)")
    println("Wrote $(DEGRADATION_REL)")
end

main()
