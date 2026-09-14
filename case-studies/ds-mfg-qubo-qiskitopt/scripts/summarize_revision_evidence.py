#!/usr/bin/env python3
"""Inventory every cached VQE seed campaign and preserve missing cost data.

Standard library only. Timings in different cohorts/implementations are never
combined into a synthetic historical end-to-end cost.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUTPUT = ROOT / "ds_mfg_revision_audit"


def read_rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_rows(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def wilson(hits, reads):
    z = 1.959963984540054
    p = hits / reads
    denominator = 1 + z*z/reads
    center = (p + z*z/(2*reads))/denominator
    radius = z*math.sqrt(p*(1-p)/reads+z*z/(4*reads*reads))/denominator
    return max(0., center-radius), min(1., center+radius)


def tts(hits, reads, seconds):
    if hits == 0:
        return math.inf
    return seconds/reads * math.ceil(math.log(.01)/math.log1p(-hits/reads))


def vqe_inventory():
    rows, cohorts, sources = [], [], []
    for path in sorted(ROOT.glob("ds_mfg_vqe*/*summary.csv")):
        sources.append(path)
        campaign = path.parent.name
        data = read_rows(path)
        for row in data:
            rows.append(dict(campaign=campaign, source=path.relative_to(ROOT).as_posix(),
                seed=int(row["seed"]), optimizer_reads=int(row["optimizer_reads"]),
                maximum_iterations=int(row["maximum_iterations"]),
                final_reads=int(row["total_reads"]),
                global_hits=int(row.get("global_hits",row.get("global_reads",0))),
                feasible_hits=int(row.get("gurobi_pool_feasible_hits",row.get("pool_reads",0))),
                top50_hits=row.get("top50_hits",""), solve_seconds=float(row["solve_time_sec"])))
        group = [row for row in rows if row["campaign"] == campaign]
        reads = sum(row["final_reads"] for row in group)
        hits = sum(row["global_hits"] for row in group)
        feasible = sum(row["feasible_hits"] for row in group)
        # Explicit accurate summation keeps artifacts stable across Python's
        # change to the built-in float sum algorithm in Python 3.12.
        seconds = math.fsum(row["solve_seconds"] for row in group)
        lo,hi = wilson(hits,reads)
        flo,fhi = wilson(feasible,reads)
        cohorts.append(dict(campaign=campaign, runs=len(group), reads=reads, global_hits=hits,
            seeds_with_global_hit=sum(row["global_hits"]>0 for row in group),
            minimum_global_hits=min(row["global_hits"] for row in group),
            maximum_global_hits=max(row["global_hits"] for row in group),
            feasible_hits=feasible, solve_seconds=seconds,
            global_rate=hits/reads, global_wilson_low=lo, global_wilson_high=hi,
            feasible_rate=feasible/reads, feasible_wilson_low=flo, feasible_wilson_high=fhi,
            descriptive_global_tts99=tts(hits,reads,seconds),
            descriptive_feasible_tts99=tts(feasible,reads,seconds)))
    assert len(rows) == 42, "Re-inventory and describe any new VQE campaigns"
    # The original 128-read baseline has no summary table; retain its seed and
    # explicitly missing time rather than silently excluding that campaign.
    initial_path = ROOT/"ds_mfg_saved_distributions/vqe_distribution.csv"
    metadata_path = ROOT/"ds_mfg_saved_distributions/run_metadata.csv"
    top_path = ROOT/"ds_mfg_reduced_flow_objective/reduced_exact_top_flows.csv"
    sources.extend([initial_path,metadata_path,top_path])
    initial = read_rows(initial_path)
    metadata = {r["key"]:r["value"] for r in read_rows(metadata_path)}
    top = {r["flow_bits"]:r for r in read_rows(top_path)}
    reads = sum(int(r["reads"]) for r in initial)
    hits = sum(int(r["reads"]) for r in initial if r["flow_bits"]=="1001110100111100011")
    feasible = sum(int(r["reads"]) for r in initial
                   if top.get(r["flow_bits"],{}).get("match") in ("global_optimum","gurobi_pool"))
    rows.append(dict(campaign="ds_mfg_saved_distributions",source=initial_path.relative_to(ROOT).as_posix(),
        seed=int(metadata["VQE.seed"]),optimizer_reads=int(metadata["number_of_reads"]),
        maximum_iterations=int(metadata["maximum_iterations"]),final_reads=reads,
        global_hits=hits,feasible_hits=feasible,
        top50_hits=sum(int(r["reads"]) for r in initial if r["flow_bits"] in top),solve_seconds=""))
    lo,hi = wilson(hits,reads)
    flo,fhi = wilson(feasible,reads)
    cohorts.append(dict(campaign="ds_mfg_saved_distributions",runs=1,reads=reads,global_hits=hits,
        seeds_with_global_hit=int(hits>0),minimum_global_hits=hits,maximum_global_hits=hits,
        feasible_hits=feasible,solve_seconds="",global_rate=hits/reads,
        global_wilson_low=lo,global_wilson_high=hi,feasible_rate=feasible/reads,
        feasible_wilson_low=flo,feasible_wilson_high=fhi,
        descriptive_global_tts99="",descriptive_feasible_tts99=""))
    assert len(rows) == 43
    screen = [row for row in rows if row["campaign"]=="ds_mfg_vqe_reduced_flow_objective_v3"]
    assert {r["seed"] for r in screen} == set(range(74001,74021))
    write_rows(OUTPUT/"all_vqe_seed_runs.csv",rows)
    write_rows(OUTPUT/"vqe_cohorts.csv",cohorts)
    return sources


def historical_costs():
    source = ROOT/"ds_mfg_ibm_qaoa_pilot_fez_4096x3/summary.csv"
    hardware = read_rows(source)[0]
    rows = []
    def add(stage, value, status, evidence, scope):
        rows.append(dict(stage=stage, seconds=value, status=status, source=evidence, scope=scope))
    for stage in ("environment_initialization", "IP_to_QUBO_reformulation", "exhaustive_auxiliary_repair",
                  "surrogate_fit", "offline_angle_optimization", "logical_circuit_construction",
                  "local_sample_repair_and_scoring"):
        add(stage,"","not_recorded","historical artifact bundle",
            "unknown; do not substitute zero or a later rerun")
    for stage in ("circuit_compilation", "hardware_queue", "hardware_execution",
                  "hardware_result_retrieval_and_scoring"):
        add(stage,"","not_separately_recorded",source.relative_to(ROOT).as_posix(),
            "included in the runner aggregate where applicable; no individual timestamps")
    add("hardware_runner_aggregate",float(hardware["elapsed_sec"]),"recorded",
        source.relative_to(ROOT).as_posix(),
        "73.5-second runner interval; does not identify queue share or include offline training")
    sources = [source]
    for name,path in (
        ("energy_QAOA_local_solve","ds_mfg_qaoa_juliqaoa_transfer_highread/qaoa_juliqaoa_transfer_summary.csv"),
        ("VQE_20_seed_local_solves","ds_mfg_vqe_reduced_flow_objective_v3/vqe_reduced_top50_sampling_summary.csv")):
        sources.append(ROOT/path)
        add(name,math.fsum(float(r["solve_time_sec"]) for r in read_rows(ROOT/path)),"recorded",path,
            "solver call boundary; excludes preceding preprocessing and subsequent scoring")
    add("historical_end_to_end_total","","not_reconstructible","historical artifact bundle",
        "unknown stages prevent a numerical total; sampling TTS is not an end-to-end comparison")
    write_rows(OUTPUT/"historical_cost_ledger.csv",rows)
    return sources


def cached_qaoa():
    rows, sources = [], []
    experiments = [
        ("full_QUBO_depth","ds_mfg_fixed_angle_sweep/fixed_angle_qaoa_summary.csv"),
        ("full_QUBO_iteration_cap","ds_mfg_iteration_sweep/sweep_summary.csv"),
        ("full_QUBO_angle_perturbation","ds_mfg_p2_perturbation_sweep_v3/p2_perturbation_summary.csv"),
        ("historical_energy_transfer_depth","ds_mfg_qaoa_juliqaoa_angle_search/juliqaoa_angle_summary.csv")]
    for campaign,path in experiments:
        sources.append(ROOT/path)
        for row in read_rows(ROOT/path):
            rows.append(dict(campaign=campaign, source=path,
                p=row.get("p",row.get("number_of_layers","")), seed=row.get("seed",""),
                iteration_cap=row.get("maximum_iterations",""), sigma=row.get("sigma",""),
                reads=row.get("total_reads",""),
                optimum_hits=row.get("global_reads",row.get("projected_global_reads","")),
                optimum_probability=row.get("global_probability",""),
                seconds=row.get("solve_time_sec","")))
    write_rows(OUTPUT/"cached_qaoa_sweeps.csv",rows)
    return sources


def generate():
    OUTPUT.mkdir(exist_ok=True)
    sources = vqe_inventory()+historical_costs()+cached_qaoa()
    sources.append(Path(__file__))
    manifest = {p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
                for p in sorted(set(sources))}
    (OUTPUT/"source_checksums.json").write_text(json.dumps(manifest,indent=2)+"\n")


def main():
    global OUTPUT
    if sys.argv[1:] == ["--check"]:
        committed = OUTPUT
        with tempfile.TemporaryDirectory(prefix="ds-mfg-evidence-") as tmp:
            OUTPUT = Path(tmp)
            generate()
            for path in OUTPUT.iterdir():
                assert (committed/path.name).read_bytes() == path.read_bytes(), path.name
        print("PASS: cached seed inventory, cohort statistics, QAOA sweeps, and cost ledger reproduce.")
    elif not sys.argv[1:]:
        generate()
        print("Inventoried 43 VQE seed runs, four cached QAOA sweeps, and historical cost gaps.")
    else:
        raise SystemExit("Usage: summarize_revision_evidence.py [--check]")


if __name__ == "__main__":
    main()
