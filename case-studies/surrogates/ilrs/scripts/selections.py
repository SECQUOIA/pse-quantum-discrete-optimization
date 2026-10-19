"""Unit/ionic-pair identity, independent of routing assignments."""
import csv
from collections import defaultdict
import numpy as np
from ilrs import bits, write_json

SELECTION_ORDER = ['y(r1)', 'y(r2)', 'y(s1)', 'y(s2)', 'y(s3)',
                   'z_cat(1)', 'z_cat(2)', 'z_an(1)', 'z_an(2)']


def selection_catalog(feasible, objectives):
    groups = defaultdict(list)
    for index in np.flatnonzero(feasible):
        groups[int(index) >> 11].append(int(index))
    rows = []
    for selection, assignments in sorted(groups.items()):
        costs = objectives[assignments]
        np.testing.assert_allclose(costs, costs[0], atol=1e-9, rtol=0)
        y = bits(selection, 9)
        rows.append(dict(
            selection_bits=''.join(map(str, y)),
            units=';'.join(name for name, active in zip(['r1','r2','s1','s2','s3'], y[:5]) if active),
            cation=int(np.flatnonzero(y[5:7])[0])+1,
            anion=int(np.flatnonzero(y[7:9])[0])+1,
            objective=float(costs[0]),
            feasible_flow_assignments=len(assignments)))
    if len(rows) != 84:
        raise AssertionError('Expected 3 reactor subsets × 7 separator subsets × 4 ionic pairs')
    return rows


def save_catalog(path, rows):
    with (path/'feasible_selections.csv').open('w', newline='') as f:
        writer=csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    return dict(
        selection_order=SELECTION_ORDER,
        feasible_selection_count=len(rows),
        feasible_flow_assignment_count=sum(r['feasible_flow_assignments'] for r in rows),
        interpretation='Unit selections and ionic pair define an option; feasible routes are grouped within it.',
        coverage_policy='Count a selection as covered only when an originally feasible sampled assignment represents it.')
