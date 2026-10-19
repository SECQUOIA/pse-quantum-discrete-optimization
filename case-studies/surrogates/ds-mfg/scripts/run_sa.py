#!/usr/bin/env python3
"""DS-MFG local simulated annealing with paired input/result artifacts."""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import importlib.util
import inspect
import io
import json
import platform
import subprocess
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import dimod
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'inputs'
REDUCED = ROOT / 'model'
ARCHIVE = SOURCE / 'full_qubo.zip'
PROTOCOL_REVISION = 'a2cb3124e07df04d984454034f913e33ef5adcfd'
import sys
sys.path.insert(0, str(ROOT.parent / 'scripts'))
from run_io import (config, model_dict, create_run, read_json, result_template,
                    save_result, add_samples, check_result, timing, utc_now, fail_run, implementation)

OPTIMUM_BITS = '1001110100111100011'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def assignments(n):
    return ((np.arange(2**n, dtype=np.uint64)[:, None] >>
             np.arange(n - 1, -1, -1, dtype=np.uint64)) & 1).astype(np.int8)


class Model:
    def __init__(self, raw):
        self.hashes = {k: digest(v) for k, v in raw.items()}
        self.n, self.scale, self.offset = np.loadtxt(
            io.BytesIO(raw['scalars.csv']), delimiter=',', skiprows=1)
        self.n = int(self.n)
        if raw['L_vector.csv'].startswith(b'index,coefficient'):
            table = np.loadtxt(io.BytesIO(raw['L_vector.csv']), delimiter=',', skiprows=1)
            np.testing.assert_array_equal(table[:, 0], np.arange(1, self.n + 1))
            self.linear = table[:, 1]
            rows = list(csv.reader(io.StringIO(raw['Q_matrix.csv'].decode())))
            labels = [f'x{i}' for i in range(1, self.n + 1)]
            if rows[0][1:] != labels or [r[0] for r in rows[1:]] != labels:
                raise ValueError('Unexpected surrogate matrix variable order')
            self.q = np.array([r[1:] for r in rows[1:]], dtype=float)
        else:
            self.linear = np.loadtxt(io.BytesIO(raw['L_vector.csv']), delimiter=',')
            self.q = np.loadtxt(io.BytesIO(raw['Q_matrix.csv']), delimiter=',')
        if self.linear.shape != (self.n,) or self.q.shape != (self.n, self.n):
            raise ValueError('Coefficient dimensions do not match scalars')
        if self.scale <= 0 or not all(np.isfinite(v).all() for v in
                                     [self.linear, self.q, np.array([self.scale, self.offset])]):
            raise ValueError('Expected finite coefficients and positive scale')
        # Binary x_i^2=x_i; combine both matrix triangles exactly once.
        linear = {i: float(self.scale * (self.linear[i] + self.q[i, i]))
                  for i in range(self.n)}
        quadratic = {(i, j): float(self.scale * (self.q[i, j] + self.q[j, i]))
                     for i in range(self.n) for j in range(i + 1, self.n)
                     if self.q[i, j] + self.q[j, i] != 0}
        self.bqm = dimod.BinaryQuadraticModel(linear, quadratic,
                                             self.scale * self.offset, dimod.BINARY)

    def energies(self, x):
        return self.scale * (self.offset + x @ self.linear +
                            np.einsum('bi,ij,bj->b', x, self.q, x, optimize=False))


def load_models():
    provenance = json.loads((ROOT/'inputs/provenance.json').read_text())
    for path, expected in provenance['sha256'].items():
        if digest((ROOT/path).read_bytes()) != expected:
            raise ValueError(f'Frozen input changed: {path}')
    names = ['scalars.csv', 'L_vector.csv', 'Q_matrix.csv']
    with zipfile.ZipFile(ARCHIVE) as z:
        full = Model({k: z.read(k) for k in names})
    surrogate = Model({k: (REDUCED / ('reduced_' + k)).read_bytes() for k in names})
    if (full.n, surrogate.n) != (36, 19):
        raise ValueError('Expected 36 full and 19 surrogate variables')
    return full, surrogate


def load_ip():
    path = SOURCE / 'evaluate_design.py'
    spec = importlib.util.spec_from_file_location('ds_mfg_ip', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Scorer:
    def __init__(self, full, surrogate):
        self.full, self.surrogate, self.ip = full, surrogate, load_ip()
        with (REDUCED / 'auxiliary_components.csv').open() as f:
            self.components = [np.array([int(i) - 1 for i in r['aux_indices'].split(';')])
                               for r in csv.DictReader(f)]
        if sorted(np.concatenate(self.components).tolist()) != list(range(19, 36)):
            raise ValueError('Auxiliary partition invalid')
        for a, indices in enumerate(self.components):
            for other in self.components[a + 1:]:
                if np.any(full.q[np.ix_(indices, other)] + full.q[np.ix_(other, indices)].T):
                    raise ValueError('Auxiliary components have cross coupling')
        self.cache = {}

    def repair(self, y):
        x = np.zeros(36, dtype=np.int8)
        x[:19] = y
        for indices in self.components:
            candidates = np.repeat(x[None, :], 2**len(indices), axis=0)
            candidates[:, indices] = assignments(len(indices))
            x = candidates[np.argmin(self.full.energies(candidates))]
        return x, float(self.full.energies(x[None, :])[0])

    def score(self, bits):
        key = ''.join(map(str, bits[:19]))
        if key not in self.cache:
            y = np.array(list(map(int, key)), dtype=np.int8)
            repaired, energy = self.repair(y)
            assign = dict(zip(self.ip.FLOW_ORDER, map(int, y)))
            feasible = self.ip.is_feasible(assign)
            self.cache[key] = dict(
                flow_bits=key, surrogate_energy=float(self.surrogate.energies(y[None, :])[0]),
                repaired_energy=energy, repaired_full_bits=''.join(map(str, repaired)),
                feasible=bool(feasible), feasible_cost=self.ip.objective(assign) if feasible else None,
                optimum_hit=key == OPTIMUM_BITS)
        return self.cache[key]


def validate(full, surrogate, scorer):
    """Exhaustive conversion, independent brute-force repair, and cached truth checks."""
    max_error = 0.0
    best = (float('inf'), None)
    optimum_energy = float(surrogate.energies(np.array([list(map(int, OPTIMUM_BITS))]))[0])
    rank = 1
    for start in range(0, 2**19, 8192):
        x = ((np.arange(start, start + 8192, dtype=np.uint64)[:, None] >>
              np.arange(18, -1, -1, dtype=np.uint64)) & 1).astype(np.int8)
        expected = surrogate.energies(x)
        actual = surrogate.bqm.energies((x, list(range(19))))
        max_error = max(max_error, float(np.max(np.abs(actual - expected))))
        rank += int(np.sum(expected < optimum_energy - 1e-8))
        i = int(np.argmin(expected))
        if expected[i] < best[0]:
            best = (float(expected[i]), x[i].copy())
    if max_error > 1e-8 or rank != 33:
        raise AssertionError(f'Surrogate conversion/rank mismatch: {max_error}, {rank}')
    rng = np.random.default_rng(74018)
    probes = rng.integers(0, 2, size=(100, 36), dtype=np.int8)
    np.testing.assert_allclose(full.bqm.energies((probes, list(range(36)))),
                               full.energies(probes), atol=1e-8, rtol=0)
    # Enumerate all 2^17 auxiliary choices independently for three flows.
    for y in [np.zeros(19, dtype=np.int8), np.array(list(map(int, OPTIMUM_BITS))), probes[0, :19]]:
        x = np.column_stack((np.tile(y, (2**17, 1)), assignments(17)))
        np.testing.assert_allclose(scorer.repair(y)[1], np.min(full.energies(x)), atol=1e-8, rtol=0)
    with (REDUCED / 'reduced_exact_top_flows.csv').open() as f:
        rows = list(csv.DictReader(f))
    feasible_count = 0
    for row in rows:
        score = scorer.score(list(map(int, row['flow_bits'])))
        np.testing.assert_allclose(score['repaired_energy'], float(row['exact_repaired_qubo_energy']), atol=1e-8, rtol=0)
        if score['feasible']:
            feasible_count += 1
            np.testing.assert_allclose(score['feasible_cost'], score['repaired_energy'], atol=1e-8, rtol=0)
        if score['feasible'] != (row['match'] in ['global_optimum', 'gurobi_pool']):
            raise AssertionError('Original constraints disagree with cached feasible pool')
    if feasible_count != 36:
        raise AssertionError('Expected 36 cached feasible designs')
    minimum = scorer.score(best[1])
    np.testing.assert_allclose(minimum['repaired_energy'], 14.5505, atol=1e-8, rtol=0)
    return dict(surrogate_states_checked=2**19, maximum_conversion_error=max_error,
                full_conversion_probes=100, brute_force_repair_flows=3,
                cached_top_flows_checked=len(rows), feasible_pool_count=feasible_count,
                true_optimum_surrogate_rank=rank, surrogate_minimum=minimum)


def source_audit(path, full):
    if path is None:
        return {'status': 'not_requested'}
    directory = path / 'ds-mfg/discrete_ip'
    other = Model({k: (directory / 'julia_exports' / k).read_bytes() for k in full.hashes})
    # Compare canonical objective, allowing formatting and floating-point roundoff.
    np.testing.assert_allclose(full.linear, other.linear, atol=1e-10, rtol=0)
    np.testing.assert_allclose(full.q, other.q, atol=1e-10, rtol=0)
    np.testing.assert_allclose([full.scale, full.offset], [other.scale, other.offset], atol=1e-10, rtol=0)
    metadata = []
    for file in sorted((directory / 'result_raw').glob('*results*.json')):
        if not file.name.startswith(('SA_', 'Dwave_')):
            continue
        data = json.loads(file.read_text())
        metadata.append(dict(file=file.name, sha256=digest(file.read_bytes()),
                             solver=data.get('solver'), num_samples=data.get('num_samples'),
                             has_assignments=bool(data.get('solutions') or data.get('samples'))))
    revision = subprocess.check_output(['git', '-C', str(path), 'rev-parse', 'HEAD'], text=True).strip()
    dirty = subprocess.check_output(['git', '-C', str(path), 'status', '--porcelain'], text=True)
    return dict(path=str(path), revision=revision, worktree_dirty=bool(dirty),
                coefficient_match=True, coefficient_hashes=other.hashes, results=metadata)


def run(model, name, scorer, output, seed, validation, provenance, shared_setup):
    import neal
    kwargs = {'num_reads': 1000}
    if seed is not None:
        kwargs['seed'] = seed
    started = time.perf_counter()
    inp = config('ds-mfg', 'full_qubo' if name == 'full' else 'surrogate', 'sa',
                 model_dict(model.bqm), kwargs, provenance=provenance)
    folder = create_run(output, inp)
    result = result_template(read_json(folder / 'input.json'), 'RUNNING')
    result['timestamps']['started_utc'] = utc_now()
    save_result(folder, result)
    try:
        call_started = time.perf_counter()
        samples = neal.SimulatedAnnealingSampler().sample(model.bqm, **kwargs)
        sampling_seconds = time.perf_counter() - call_started
        post_started = time.perf_counter()
        rows = []
        for datum in samples.data(['sample', 'energy', 'num_occurrences'], sorted_by=None):
            bits = [int(datum.sample[i]) for i in range(model.n)]
            raw = float(model.energies(np.array([bits]))[0])
            np.testing.assert_allclose(datum.energy, raw, atol=1e-8, rtol=0)
            rows.append(dict(sample_bits=''.join(map(str, bits)), occurrences=int(datum.num_occurrences),
                             sampled_model_energy=raw, **scorer.score(bits),
                             auxiliary_gap=raw-scorer.score(bits)['repaired_energy'] if name=='full' else None))
        reads = sum(r['occurrences'] for r in rows)
        if reads != 1000:
            raise AssertionError('Read count mismatch')
        feasible = [r for r in rows if r['feasible']]
        summary = dict(formulation=name, num_reads=reads,
                       optimum_hits=sum(r['occurrences'] for r in rows if r['optimum_hit']),
                       feasible_hits=sum(r['occurrences'] for r in feasible),
                       unique_feasible_designs=len({r['flow_bits'] for r in feasible}),
                       best_feasible_cost=min((r['feasible_cost'] for r in feasible), default=None))
        summary.update(optimum_rate=summary['optimum_hits']/reads, feasible_rate=summary['feasible_hits']/reads)
        add_samples(result, rows, summary)
        result.update(status='COMPLETED', validation=validation, info={'sampler_info': samples.info},
                      execution_time=sampling_seconds, execution_time_basis='local sampler-call wallclock')
        result['timing'] = timing(sampling_seconds=sampling_seconds,
                                 postprocessing_seconds=time.perf_counter()-post_started,
                                 wallclock_seconds=time.perf_counter()-started,
                                 shared_setup_seconds=shared_setup['seconds'], shared_setup_id=shared_setup['id'],
                                 measurement_notes='Per-run wallclock excludes shared batch setup and final JSON write; shared setup is counted once per batch.')
        result['timestamps']['completed_utc'] = utc_now()
        save_result(folder, check_result(result))
    except Exception as error:
        fail_run(folder, started, error)
        raise
    print(json.dumps({'run': str(folder), 'metrics': result['metrics']}))


def main():
    import neal
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'results', help='Results root; solver/formulation/timestamp generated automatically')
    parser.add_argument('--pd-ising', type=Path, help='Optional read-only baseline provenance check')
    parser.add_argument('--validate-only', action='store_true', help='Validate locally without saving a solver run')
    parser.add_argument('--repeats', type=int, default=1)
    parser.add_argument('--seed', type=int)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error('--repeats must be positive')
    if args.seed is not None and not 0 <= args.seed <= 2**31 - args.repeats:
        parser.error('Invalid seed/repeat range')
    started = time.perf_counter()
    full, surrogate = load_models()
    scorer = Scorer(full, surrogate)
    audit = source_audit(args.pd_ising, full)
    validation = validate(full, surrogate, scorer)
    if args.validate_only:
        print(json.dumps(validation, indent=2))
        return
    provenance = dict(packages={p: importlib.metadata.version(p) for p in ['numpy','dimod','dwave-neal','dwave-samplers']},
                      implementation=implementation(__file__),
                      sampler_source=inspect.getsource(neal.SimulatedAnnealingSampler.sample),
                      input_provenance=read_json(ROOT/'inputs/provenance.json'),
                      full_coefficient_hashes=full.hashes, surrogate_coefficient_hashes=surrogate.hashes,
                      baseline_audit=audit, seed_policy='automatic' if args.seed is None else 'explicit override')
    shared = {'seconds': time.perf_counter()-started, 'id': utc_now()}
    for repeat in range(args.repeats):
        for name, model in [('full', full), ('surrogate', surrogate)]:
            run(model, name, scorer, args.output, None if args.seed is None else args.seed+repeat,
                validation, provenance, shared)


if __name__ == '__main__':
    main()
