#!/usr/bin/env julia

const STUDY_ROOT = abspath(joinpath(@__DIR__, ".."))
const REPORT_DIR = joinpath(STUDY_ROOT, "ds_mfg_hit_rate_reports")

include(joinpath(@__DIR__, "hit_rate_stats.jl"))

const STAT_COLUMNS = hit_rate_stat_headers()
const STAT_COLUMN_SET = Set(STAT_COLUMNS)
const SUMMARY_TARGETS = [
    (
        relpath = "ds_mfg_qaoa_juliqaoa_transfer/qaoa_juliqaoa_transfer_summary.csv",
        time_column = "solve_time_sec",
    ),
    (
        relpath = "ds_mfg_qaoa_juliqaoa_transfer_highread/qaoa_juliqaoa_transfer_summary.csv",
        time_column = "solve_time_sec",
    ),
    (
        relpath = "ds_mfg_vqe_reduced_flow_objective_v3/vqe_reduced_top50_sampling_summary.csv",
        time_column = "solve_time_sec",
    ),
    (
        relpath = "ds_mfg_vqe_reduced_flow_objective_v4/vqe_reduced_top50_sampling_summary.csv",
        time_column = "solve_time_sec",
    ),
    (
        relpath = "ds_mfg_vqe_reduced_flow_objective_final/vqe_reduced_top50_sampling_summary.csv",
        time_column = "solve_time_sec",
    ),
    (
        relpath = "ds_mfg_classical_baselines/classical_baseline_summary.csv",
        time_column = "wall_time_sec",
    ),
]

function read_csv_table(path::AbstractString)
    lines = filter(line -> !isempty(strip(line)), readlines(path))
    length(lines) >= 2 || error("Expected header and at least one row in $(path)")
    header = String.(split(lines[1], ","; keepempty = true))
    rows = [String.(split(line, ","; keepempty = true)) for line in lines[2:end]]
    for (row_index, row) in enumerate(rows)
        length(row) == length(header) ||
            error("Row $(row_index + 1) in $(path) has $(length(row)) fields; expected $(length(header))")
    end
    return header, rows
end

function write_csv_table(path::AbstractString, header::Vector{String}, rows::Vector{Vector{String}})
    open(path, "w") do io
        println(io, join(header, ','))
        for row in rows
            println(io, join(row, ','))
        end
    end
end

function row_dict(header::Vector{String}, row::Vector{String})
    return Dict(header[index] => row[index] for index in eachindex(header))
end

function remove_stat_columns(header::Vector{String}, row::Vector{String})
    keep_indices = [index for (index, name) in enumerate(header) if !(name in STAT_COLUMN_SET)]
    return header[keep_indices], row[keep_indices]
end

function required_int(row::Dict{String,String}, key::String)
    haskey(row, key) || error("Missing required column $(key)")
    return parse(Int, row[key])
end

function required_float(row::Dict{String,String}, key::String)
    haskey(row, key) || error("Missing required column $(key)")
    return parse(Float64, row[key])
end

function required_total_trials(row::Dict{String,String})
    haskey(row, "total_reads") && return parse(Int, row["total_reads"])
    haskey(row, "total_samples") && return parse(Int, row["total_samples"])
    error("Missing required column total_reads or total_samples")
end

function run_label(row::Dict{String,String})
    keys = [
        "angle_target",
        "p",
        "seed",
        "optimizer_reads",
        "final_reads",
        "sample_budget",
        "maximum_iterations",
    ]
    parts = String[]
    for key in keys
        haskey(row, key) && !isempty(row[key]) && push!(parts, "$(key)=$(row[key])")
    end
    return join(parts, ';')
end

function report_rows_for_summary(relpath::String, time_column::String, header::Vector{String}, rows)
    report_rows = Vector{String}[]
    for row_values in rows
        row = row_dict(header, row_values)
        total_reads = required_total_trials(row)
        elapsed_sec = required_float(row, time_column)
        label = run_label(row)
        algorithm = get(row, "algorithm", "")
        for (event, hit_column) in HIT_RATE_EVENTS
            hits = required_int(row, hit_column)
            stats = hit_rate_event_stats(hits, total_reads, elapsed_sec)
            push!(
                report_rows,
                [
                    relpath,
                    algorithm,
                    label,
                    time_column,
                    format_stat_value(elapsed_sec),
                    event,
                    string(hits),
                    string(total_reads),
                    format_stat_value(stats.hit_rate),
                    format_stat_value(stats.hit_rate_wilson95_low),
                    format_stat_value(stats.hit_rate_wilson95_high),
                    format_stat_value(stats.tts95_sec),
                ],
            )
        end
    end
    return report_rows
end

function update_summary!(target)
    path = joinpath(STUDY_ROOT, target.relpath)
    header, rows = read_csv_table(path)
    updated_rows = Vector{String}[]
    updated_header = nothing

    for row in rows
        base_header, base_row = remove_stat_columns(header, row)
        row = row_dict(base_header, base_row)
        total_reads = required_total_trials(row)
        top50_hits = required_int(row, "top50_hits")
        top10_hits = required_int(row, "top10_hits")
        global_hits = required_int(row, "global_hits")
        elapsed_sec = required_float(row, target.time_column)
        stat_values = hit_rate_stat_values(total_reads, top50_hits, top10_hits, global_hits, elapsed_sec)

        insert_after = findfirst(==("global_hits"), base_header)
        isnothing(insert_after) && error("Missing global_hits column in $(target.relpath)")

        if isnothing(updated_header)
            updated_header = vcat(base_header[1:insert_after], STAT_COLUMNS, base_header[(insert_after + 1):end])
        end
        push!(updated_rows, vcat(base_row[1:insert_after], stat_values, base_row[(insert_after + 1):end]))
    end

    write_csv_table(path, updated_header, updated_rows)
    return report_rows_for_summary(target.relpath, target.time_column, updated_header, updated_rows)
end

function main()
    mkpath(REPORT_DIR)
    report_rows = Vector{String}[]
    for target in SUMMARY_TARGETS
        append!(report_rows, update_summary!(target))
    end

    report_header = [
        "source_summary",
        "algorithm",
        "run_label",
        "time_basis_column",
        "time_basis_sec",
        "event",
        "hits",
        "total_reads",
        "hit_rate",
        "hit_rate_wilson95_low",
        "hit_rate_wilson95_high",
        "tts95_sec",
    ]
    report_path = joinpath(REPORT_DIR, "time_to_solution_report.csv")
    write_csv_table(report_path, report_header, report_rows)
    println("Updated cached hit-rate summaries.")
    println("Time-to-solution report written to: ", report_path)
end

main()
