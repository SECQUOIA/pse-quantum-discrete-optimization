#!/usr/bin/env python3
"""Rerun the DS-MFG 19-flow IP with Gurobi for provenance.

This script is intentionally independent of the external source repository. It
builds the same 19-binary flow IP used by the retained Gurobi export, runs
Gurobi solution-pool search, and writes sanitized metadata without license
paths, machine-local paths, or credentials.
"""

from __future__ import annotations

import csv
import json
import math
import os
import zipfile
from pathlib import Path

import gurobipy as gp
from gurobipy import GRB

from enumerate_ip_provenance import FLOW_ORDER, FLOW_ROWS, NETWORK, POOL_ARCHIVE, POOL_MEMBER, UO


STUDY_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = STUDY_ROOT / "ds_mfg_gurobi_provenance"
SUMMARY_PATH = OUTPUT_DIR / "gurobi_local_pool_rerun_summary.json"
SOLUTIONS_PATH = OUTPUT_DIR / "gurobi_local_pool_rerun_solutions.csv"


def json_number(value):
    if value is None:
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return value
    return numeric if math.isfinite(numeric) else None


def status_text(status_code):
    for name in dir(GRB):
        if name.isupper() and getattr(GRB, name) == status_code:
            if name in {"OPTIMAL", "INFEASIBLE", "INF_OR_UNBD", "UNBOUNDED", "TIME_LIMIT"}:
                return name
    return str(status_code)


def retained_pool_records():
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


def operating_modes(bits):
    return [
        UO[flow]
        for flow, bit in zip(FLOW_ORDER, bits)
        if bit == "1" and UO[flow] not in ("none", "holdT", "noholdT")
    ]


def add_flow_conservation(model, variables):
    nodes = sorted({node for flow in FLOW_ORDER for node in NETWORK[flow]})
    source_node = nodes[0]
    sink_node = nodes[-1]
    for node in nodes:
        inflow = gp.quicksum(variables[flow] for flow in FLOW_ORDER if NETWORK[flow][1] == node)
        outflow = gp.quicksum(variables[flow] for flow in FLOW_ORDER if NETWORK[flow][0] == node)
        if node == source_node:
            model.addConstr(outflow == 1, name=f"flow_conservation[{node}]")
        elif node == sink_node:
            model.addConstr(inflow == 1, name=f"flow_conservation[{node}]")
        else:
            model.addConstr(inflow == outflow, name=f"flow_conservation[{node}]")


def build_model():
    model = gp.Model("capex")
    variables = {
        flow: model.addVar(vtype=GRB.BINARY, name=f"f({flow})")
        for flow in FLOW_ORDER
    }
    model.setObjective(
        gp.quicksum(cost * variables[flow] for flow, _, _, cost, _ in FLOW_ROWS),
        GRB.MINIMIZE,
    )
    add_flow_conservation(model, variables)

    model.addConstr((1 - variables["f01"]) + (1 - variables["f10"]) + variables["f06"] >= 1)
    model.addConstr((1 - variables["f02"]) + (1 - variables["f10"]) + variables["f06"] >= 1)
    model.addConstr((1 - variables["f03"]) + (1 - variables["f10"]) + (1 - variables["f06"]) >= 1)
    model.addConstr(variables["f01"] + variables["f02"] + variables["f05"] >= 1)
    model.addConstr(variables["f10"] + variables["f05"] >= 1)
    for flow in ("f13", "f14", "f15"):
        model.addConstr((1 - variables[flow]) + variables["f17"] >= 1)

    return model, variables


def collect_solutions(model, variables):
    rows = []
    for solution_number in range(model.SolCount):
        model.Params.SolutionNumber = solution_number
        bits = "".join(str(round(variables[flow].Xn)) for flow in FLOW_ORDER)
        rows.append(
            {
                "rank": solution_number + 1,
                "objective": model.PoolObjVal,
                "flow_bits": bits,
            }
        )
    rows.sort(key=lambda row: (row["objective"], row["flow_bits"]))
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank
    return rows


def write_solutions(rows, retained_pool):
    with SOLUTIONS_PATH.open("w", newline="") as handle:
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
        for row in rows:
            retained = retained_pool.get(row["flow_bits"], {})
            writer.writerow(
                [
                    row["rank"],
                    format(row["objective"], ".12g"),
                    row["flow_bits"],
                    retained.get("solution_id", ""),
                    retained.get("ip_obj_value", ""),
                    ";".join(operating_modes(row["flow_bits"])),
                ]
            )


def write_summary(model, rows, retained_pool, pool_solutions):
    solution_bits = {row["flow_bits"] for row in rows}
    retained_bits = set(retained_pool)
    summary = {
        "schema_version": 1,
        "artifact_role": "ds_mfg_local_gurobi_pool_rerun",
        "script": "scripts/run_gurobi_pool_provenance.py",
        "credential_policy": "No solver license paths, account files, credentials, or machine-local private paths are stored.",
        "solver": {
            "gurobi_version": ".".join(str(part) for part in gp.gurobi.version()),
            "python_api": "gurobipy",
        },
        "configuration": {
            "pool_search_mode": 2,
            "pool_solutions": pool_solutions,
            "output_flag": 0,
        },
        "model_info": {
            "model_name": model.ModelName,
            "number_of_variables": model.NumVars,
            "number_of_constraints": model.NumConstrs,
            "number_of_objectives": model.NumObj,
            "objective_sense": "Min" if model.ModelSense == 1 else "Max",
            "is_mip": bool(model.IsMIP),
            "number_of_binary_variables": model.NumBinVars,
            "number_of_integer_variables": model.NumIntVars,
            "number_of_nonzeros": model.NumNZs,
        },
        "solve_metadata": {
            "status_code": model.Status,
            "status_text": status_text(model.Status),
            "objective_value": json_number(model.ObjVal if model.SolCount else None),
            "objective_bound": json_number(model.ObjBound),
            "mip_gap": json_number(model.MIPGap if model.SolCount else None),
            "pool_obj_bound": json_number(model.PoolObjBound),
            "solution_count": model.SolCount,
            "solve_time_sec": json_number(model.Runtime),
        },
        "retained_pool_comparison": {
            "source_archive": POOL_ARCHIVE.name,
            "source_member": POOL_MEMBER,
            "retained_pool_solution_count": len(retained_pool),
            "local_rerun_solution_count": len(rows),
            "local_rerun_matches_retained_pool": solution_bits == retained_bits,
            "missing_from_local_rerun": sorted(retained_bits - solution_bits),
            "extra_in_local_rerun": sorted(solution_bits - retained_bits),
        },
        "artifacts": {
            "solutions_csv": "ds_mfg_gurobi_provenance/gurobi_local_pool_rerun_solutions.csv",
        },
        "claim_boundary": "This is a local Gurobi rerun under the recorded Gurobi version. It validates the model/pool behavior for the committed constraints but does not recover the original run's unrecorded status fields.",
    }
    SUMMARY_PATH.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    pool_solutions = int(os.environ.get("DSMFG_GUROBI_POOL_SOLUTIONS", "100"))
    retained_pool = retained_pool_records()
    model, variables = build_model()
    model.Params.OutputFlag = 0
    model.Params.PoolSearchMode = 2
    model.Params.PoolSolutions = pool_solutions
    model.optimize()

    rows = collect_solutions(model, variables)
    write_solutions(rows, retained_pool)
    write_summary(model, rows, retained_pool, pool_solutions)
    print(f"Gurobi status: {status_text(model.Status)}")
    print(f"Solutions found: {len(rows)}")
    print(f"Optimum: {rows[0]['objective']:.12g} {rows[0]['flow_bits']}")
    print(f"Summary written to: {SUMMARY_PATH.relative_to(STUDY_ROOT)}")


if __name__ == "__main__":
    main()
