"""Common original-design metrics for ILRS samplers."""
import hashlib
import json
import numpy as np
from ilrs import ROOT, Problem, index_of
from selections import selection_catalog


class Scorer:
    def __init__(self,model=ROOT/'model'):
        self.path=model;self.problem=Problem();self.model_metadata=json.loads((model/'model.json').read_text());self.summary=self.model_metadata['summary']
        for name, expected in self.model_metadata.get('metadata_sections_sha256', {}).items():
            if hashlib.sha256(json.dumps(self.model_metadata[name],sort_keys=True).encode()).hexdigest()!=expected:raise ValueError('Frozen model metadata changed')
        for name,expected in (self.summary['coefficients_sha256'] | self.summary['model_artifacts_sha256']).items():
            if hashlib.sha256((model/name).read_bytes()).hexdigest()!=expected:raise ValueError('Frozen surrogate coefficients changed')
        self.repaired=np.load(model/'repaired_energies.npy',mmap_mode='r')
        self.surrogate=np.load(model/'surrogate_energies.npy',mmap_mode='r')
        self.feasible=np.load(model/'feasible.npy',mmap_mode='r')
        self.objectives=np.load(model/'objectives.npy',mmap_mode='r')
        self.optima=set(self.summary['optimum_indices'])
        self.selections=selection_catalog(self.feasible,self.objectives)
        self.linear=np.loadtxt(model/'linear.csv',delimiter=',');self.q=np.loadtxt(model/'quadratic.csv',delimiter=',')
        self.offset=self.summary['offset']

    def score_counts(self,counts,algorithm,seconds,metadata,full=False):
        rows=[]
        for key,count in sorted(counts.items()):
            x=np.array(list(map(int,key)),dtype=np.int8)
            if len(x)!=(40 if full else 20):raise ValueError('Invalid bitstring length')
            i=index_of(x[:20]);repaired,completion=self.problem.repair(x[:20],return_bits=True)
            np.testing.assert_allclose(repaired,self.repaired[i],atol=1e-8,rtol=0)
            feasible=bool(self.problem.feasible(x[:20]));objective=float(self.problem.objective(x[:20]))
            if feasible!=bool(self.feasible[i]):raise AssertionError('Feasibility table mismatch')
            raw=float(self.problem.energy(x)) if full else float(self.surrogate[i])
            row=dict(sample_bits=key,occurrences=int(count),design_bits=key[:20],selection_bits=key[:9],sampled_model_energy=raw,
                     surrogate_energy=float(self.surrogate[i]),repaired_energy=float(repaired),
                     repaired_full_bits=''.join(map(str,completion)),feasible=feasible,
                     feasible_cost=objective if feasible else None,optimum_hit=i in self.optima,
                     auxiliary_gap=raw-float(repaired) if full else None)
            if full and row['auxiliary_gap'] < -1e-8:raise AssertionError('Repair increases original energy')
            rows.append(row)
        reads=sum(counts.values());feasible_rows=[r for r in rows if r['feasible']]
        groups={}
        for row in feasible_rows:
            key=row['selection_bits']
            if key not in groups:
                groups[key]=dict(selection_bits=key,occurrences=0,distinct_flow_assignments=set(),feasible_cost=row['feasible_cost'],optimum_hit=row['optimum_hit'])
            groups[key]['occurrences']+=row['occurrences']
            groups[key]['distinct_flow_assignments'].add(row['design_bits'])
        grouped=[]
        for key,group in sorted(groups.items()):
            grouped.append(group | {'distinct_flow_assignments':len(group['distinct_flow_assignments'])})
        summary=dict(algorithm=algorithm,reads=reads,optimum_hits=sum(r['occurrences'] for r in rows if r['optimum_hit']),
                     feasible_hits=sum(r['occurrences'] for r in feasible_rows),
                     distinct_feasible_designs=len({r['design_bits'] for r in feasible_rows}),
                     distinct_feasible_flow_assignments=len({r['design_bits'] for r in feasible_rows}),
                     distinct_feasible_selections=len(groups),
                     total_feasible_selections=len(self.selections),
                     selection_coverage=len(groups)/len(self.selections),
                     best_feasible_cost=min((r['feasible_cost'] for r in feasible_rows),default=None),
                     sampling_seconds=seconds,metadata=metadata)
        return rows, grouped, summary
