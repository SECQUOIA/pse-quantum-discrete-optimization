using Printf

const HIT_RATE_CONFIDENCE_LEVEL = 0.95
const HIT_RATE_WILSON_Z = 1.959963984540054
const TIME_TO_SOLUTION_CONFIDENCE = 0.95
const HIT_RATE_EVENTS = (
    ("top50", "top50_hits"),
    ("top10", "top10_hits"),
    ("global", "global_hits"),
)

function validate_hit_count(hits::Integer, total::Integer)
    total >= 0 || error("Total trials must be nonnegative, got $(total)")
    0 <= hits <= total || error("Hit count $(hits) is outside [0, $(total)]")
end

function wilson_interval(hits::Integer, total::Integer; z::Real = HIT_RATE_WILSON_Z)
    validate_hit_count(hits, total)
    total == 0 && return (low = NaN, high = NaN)

    n = Float64(total)
    phat = hits / n
    z2 = z^2
    denominator = 1.0 + z2 / n
    center = phat + z2 / (2.0 * n)
    margin = z * sqrt(phat * (1.0 - phat) / n + z2 / (4.0 * n^2))
    return (
        low = clamp((center - margin) / denominator, 0.0, 1.0),
        high = clamp((center + margin) / denominator, 0.0, 1.0),
    )
end

function time_to_solution_seconds(
    hits::Integer,
    total::Integer,
    elapsed_sec::Real;
    target_confidence::Real = TIME_TO_SOLUTION_CONFIDENCE,
)
    validate_hit_count(hits, total)
    0.0 < target_confidence < 1.0 ||
        error("Target confidence must be in (0, 1), got $(target_confidence)")
    elapsed = Float64(elapsed_sec)
    total == 0 && return NaN
    elapsed >= 0.0 || return NaN
    hits == 0 && return Inf
    hits == total && return elapsed / total

    hit_rate = hits / Float64(total)
    reads_needed = ceil(log1p(-target_confidence) / log1p(-hit_rate))
    return (elapsed / total) * reads_needed
end

function hit_rate_event_stats(hits::Integer, total::Integer, elapsed_sec::Real)
    validate_hit_count(hits, total)
    interval = wilson_interval(hits, total)
    return (
        hit_rate = total == 0 ? NaN : hits / Float64(total),
        hit_rate_wilson95_low = interval.low,
        hit_rate_wilson95_high = interval.high,
        tts95_sec = time_to_solution_seconds(hits, total, elapsed_sec),
    )
end

function format_stat_value(value)
    x = Float64(value)
    isnan(x) && return ""
    isinf(x) && return x > 0.0 ? "Inf" : "-Inf"
    return @sprintf("%.12g", x)
end

function hit_rate_stat_headers()
    headers = String[]
    for (event, _) in HIT_RATE_EVENTS
        append!(
            headers,
            [
                "$(event)_hit_rate",
                "$(event)_hit_rate_wilson95_low",
                "$(event)_hit_rate_wilson95_high",
                "$(event)_tts95_sec",
            ],
        )
    end
    return headers
end

function hit_rate_stat_values(
    total_reads::Integer,
    top50_hits::Integer,
    top10_hits::Integer,
    global_hits::Integer,
    elapsed_sec::Real,
)
    values = String[]
    for hits in (top50_hits, top10_hits, global_hits)
        stats = hit_rate_event_stats(hits, total_reads, elapsed_sec)
        append!(
            values,
            format_stat_value.([
                stats.hit_rate,
                stats.hit_rate_wilson95_low,
                stats.hit_rate_wilson95_high,
                stats.tts95_sec,
            ]),
        )
    end
    return values
end
