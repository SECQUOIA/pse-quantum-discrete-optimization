#!/usr/bin/env python3
"""Enumerate the DS-MFG 19-flow IP for provenance recovery.

This is a deterministic fallback for environments without Gurobi. It enumerates
all 2^19 binary flow assignments using the flow-conservation and logic
constraints implemented by the original Pyomo model in SECQUOIA/pd_ising.
"""

from __future__ import annotations

import csv
import itertools
import json
import zipfile
from pathlib import Path


STUDY_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = STUDY_ROOT / "ds_mfg_gurobi_provenance"
POOL_ARCHIVE = STUDY_ROOT / "Fw_ DS mfg case qubo information.zip"
POOL_MEMBER = "gurobi_pool_solutions_n36_20260126_105706.json"

FLOW_ROWS = [
    ("f00", "n00", "n01", 0.0, "none"),
    ("f01", "n01", "n02", 0.0289, "PFR"),
    ("f02", "n01", "n02", 0.0399, "CSTR"),
    ("f03", "n01", "n02", 0.1225, "batch"),
    ("f04", "n02", "n03", 0.0, "none"),
    ("f05", "n03", "n04", -0.0, "noholdT"),
    ("f06", "n03", "n04", 4.959, "holdT"),
    ("f07", "n04", "n05", 0.0, "none"),
    ("f08", "n05", "n06", 6.193, "PFR"),
    ("f09", "n05", "n06", 9.717, "CSTR"),
    ("f10", "n05", "n06", 3.352, "batch"),
    ("f11", "n06", "n07", 6.163, "batch"),
    ("f12", "n07", "n08", 0.0, "none"),
    ("f13", "n08", "n09", 0.739, "cont1"),
    ("f14", "n08", "n09", 0.84, "cont2"),
    ("f15", "n08", "n09", 1.271, "cont3"),
    ("f16", "n07", "n10", 6.99, "batchU"),
    ("f17", "n09", "n10", 1.333, "holdT"),
    ("f18", "n10", "n11", 0.0, "none"),
]

FLOW_ORDER = [row[0] for row in FLOW_ROWS]
NETWORK = {flow: (source, target) for flow, source, target, _, _ in FLOW_ROWS}
COST = {flow: cost for flow, _, _, cost, _ in FLOW_ROWS}
UO = {flow: uo for flow, _, _, _, uo in FLOW_ROWS}
NODES = sorted({node for _, source, target, _, _ in FLOW_ROWS for node in (source, target)})
SOURCE_NODE = NODES[0]
SINK_NODE = NODES[-1]


def flow_bits(assign):
    return "".join(str(assign[flow]) for flow in FLOW_ORDER)


def objective(assign):
    return sum(COST[flow] * assign[flow] for flow in FLOW_ORDER)


def operating_modes(assign):
    return [
        UO[flow]
        for flow in FLOW_ORDER
        if assign[flow] == 1 and UO[flow] not in ("none", "holdT", "noholdT")
    ]


def is_feasible(assign):
    for node in NODES:
        inflow = sum(assign[flow] for flow in FLOW_ORDER if NETWORK[flow][1] == node)
        outflow = sum(assign[flow] for flow in FLOW_ORDER if NETWORK[flow][0] == node)
        if node == SOURCE_NODE:
            if outflow != 1:
                return False
        elif node == SINK_NODE:
            if inflow != 1:
                return False
        elif inflow != outflow:
            return False

    logic_constraints = (
        (1 - assign["f01"]) + (1 - assign["f10"]) + assign["f06"] >= 1,
        (1 - assign["f02"]) + (1 - assign["f10"]) + assign["f06"] >= 1,
        (1 - assign["f03"]) + (1 - assign["f10"]) + (1 - assign["f06"]) >= 1,
        assign["f01"] + assign["f02"] + assign["f05"] >= 1,
        assign["f10"] + assign["f05"] >= 1,
        (1 - assign["f13"]) + assign["f17"] >= 1,
        (1 - assign["f14"]) + assign["f17"] >= 1,
        (1 - assign["f15"]) + assign["f17"] >= 1,
    )
    return all(logic_constraints)


def retained_pool_records():
    if not POOL_ARCHIVE.exists():
        return {}
    with zipfile.ZipFile(POOL_ARCHIVE) as archive:
        with archive.open(POOL_MEMBER) as handle:
            data = json.load(handle)

    records = {}
    for solution in data["solutions"]:
        values = solution["values"]
        bits = "".join(str(round(values[f"f({flow})"])) for flow in FLOW_ORDER)
        records[bits] = {
            "solution_id": solution["solution_id"],
            "ip_obj_value": solution["ip_obj_value"],
        }
    return records


def enumerate_feasible():
    feasible = []
    for bits in itertools.product((0, 1), repeat=len(FLOW_ORDER)):
        assign = dict(zip(FLOW_ORDER, bits))
        if is_feasible(assign):
            feasible.append(
                {
                    "objective": objective(assign),
                    "flow_bits": flow_bits(assign),
                    "operating_modes": operating_modes(assign),
                }
            )
    feasible.sort(key=lambda row: (row["objective"], row["flow_bits"]))
    return feasible


def write_csv(path, feasible, pool):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "rank",
                "ip_obj_value",
                "flow_bits",
                "retained_pool_solution_id",
                "retained_pool_obj_value",
                "operating_modes",
            ]
        )
        for rank, row in enumerate(feasible, start=1):
            pool_row = pool.get(row["flow_bits"], {})
            writer.writerow(
                [
                    rank,
                    format(row["objective"], ".12g"),
                    row["flow_bits"],
                    pool_row.get("solution_id", ""),
                    pool_row.get("ip_obj_value", ""),
                    ";".join(row["operating_modes"]),
                ]
            )


def write_summary(path, feasible, pool):
    feasible_bits = {row["flow_bits"] for row in feasible}
    pool_bits = set(pool)
    summary = {
        "schema_version": 1,
        "artifact_role": "ds_mfg_exact_ip_enumeration_provenance",
        "script": "scripts/enumerate_ip_provenance.py",
        "source_model": "Flow-conservation and logic constraints from SECQUOIA/pd_ising ds-mfg/discrete_ip/flst_opti_IP.py",
        "flow_order": FLOW_ORDER,
        "enumeration": {
            "number_of_binary_flow_variables": len(FLOW_ORDER),
            "total_binary_assignments": 2 ** len(FLOW_ORDER),
            "feasible_assignment_count": len(feasible),
            "objective_sense": "Min",
            "optimum_objective": feasible[0]["objective"],
            "optimum_flow_bits": feasible[0]["flow_bits"],
        },
        "retained_pool_comparison": {
            "source_archive": POOL_ARCHIVE.name,
            "source_member": POOL_MEMBER,
            "retained_pool_solution_count": len(pool),
            "feasible_set_matches_retained_pool": feasible_bits == pool_bits,
            "missing_from_retained_pool": sorted(feasible_bits - pool_bits),
            "extra_in_retained_pool": sorted(pool_bits - feasible_bits),
        },
        "claim_boundary": "This deterministic enumeration certifies the 19-flow IP feasible set and optimum for the encoded constraints. It does not recover Gurobi version, status code, MIPGap, or PoolGap.",
        "credential_policy": "No solver licenses, credentials, account files, or machine-local private paths are stored.",
    }
    path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    feasible = enumerate_feasible()
    pool = retained_pool_records()
    write_csv(OUTPUT_DIR / "ip_exact_feasible_flows.csv", feasible, pool)
    write_summary(OUTPUT_DIR / "ip_exact_enumeration_summary.json", feasible, pool)
    print(f"Enumerated {len(feasible)} feasible assignments.")
    print(f"Optimum: {feasible[0]['objective']:.12g} {feasible[0]['flow_bits']}")


if __name__ == "__main__":
    main()
