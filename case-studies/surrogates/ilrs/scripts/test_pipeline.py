"""Independent model, fitting, gate-convention, and artifact regression checks."""
import csv
import json
import unittest
from collections import Counter

import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit.library import DiagonalGate
from qiskit.quantum_info import Statevector

from ilrs import ROOT,Problem,bits,index_of
from build_surrogate import quadratic_projection
from run_gates import qaoa_state,design_statevector
from scoring import Scorer
from selections import selection_catalog


class PipelineTests(unittest.TestCase):
    def test_quadratic_projection_preserves_exact_quadratic(self):
        y=bits(np.arange(32),5)
        energy=7+3*y[:,0]-2*y[:,3]+11*y[:,0]*y[:,4]-5*y[:,1]*y[:,2]
        b,l,q=quadratic_projection(energy,5)
        np.testing.assert_allclose(b+y@l+np.einsum('bi,ij,bj->b',y,q,y),energy,atol=1e-12)

    def test_qaoa_against_qiskit_diagonal_circuit(self):
        energy=np.array([2.,-1.,4.,0.,3.,5.,-2.,1.]);angles=np.array([.23,-.41,.32,.17])
        circuit=QuantumCircuit(3);circuit.h(range(3))
        for beta,gamma in zip(angles[:2],angles[2:]):
            circuit.append(DiagonalGate(np.exp(-1j*gamma*energy)),range(3))
            for i in range(3):circuit.rx(2*beta,i)
        np.testing.assert_allclose(qaoa_state(energy,angles),Statevector.from_instruction(circuit).data,atol=1e-12)

    def test_vqe_design_bit_order(self):
        circuit=QuantumCircuit(3);circuit.x(0)
        state=design_statevector(circuit,[])
        self.assertEqual(int(np.argmax(np.abs(state))),4)
        np.testing.assert_array_equal(bits(4,3),[1,0,0])

    def test_feasible_pool_and_derived_product_elimination(self):
        p=Problem();y=bits(677380)
        self.assertTrue(p.feasible(y));self.assertAlmostEqual(p.objective(y),19.08)
        repaired,full=p.repair(y,return_bits=True)
        self.assertAlmostEqual(repaired,19.08);self.assertAlmostEqual(p.energy(full),repaired)
        np.testing.assert_array_equal(full[20:24],[y[5]*y[7],y[6]*y[7],y[5]*y[8],y[6]*y[8]])
        invalid=y.copy();invalid[0]=1-invalid[0]
        self.assertFalse(p.feasible(invalid))
        repaired,completion=p.repair(invalid,return_bits=True)
        np.testing.assert_array_equal(completion[:20],invalid)
        summary=json.loads((ROOT/'model/model.json').read_text())['summary']
        self.assertEqual(summary['feasible_count'],252);self.assertEqual(summary['missing_from_recorded_pool'],24)

    def test_unit_pair_selections_group_flow_variants(self):
        scorer=Scorer();catalog=selection_catalog(scorer.feasible,scorer.objectives)
        self.assertEqual(len(catalog),84)
        self.assertEqual(sum(r['feasible_flow_assignments'] for r in catalog),252)
        self.assertEqual(dict(Counter(r['feasible_flow_assignments'] for r in catalog)),{1:68,7:12,25:4})
        # Distinct ionic selections with equal costs must remain distinct options.
        self.assertLess(len({r['objective'] for r in catalog}),len(catalog))
        self.assertEqual(len({r['selection_bits'] for r in catalog}),84)
        summary=json.loads(next((ROOT/'results/historical/sa').glob('*/result.json')).read_text())['metrics']
        self.assertEqual(summary['distinct_feasible_selections'],84)
        self.assertEqual(summary['distinct_feasible_flow_assignments'],231)

    def test_future_gate_runners_save_complete_pairs(self):
        import tempfile,sys
        from pathlib import Path
        from types import SimpleNamespace
        from unittest.mock import patch
        import run_gates
        class FixedRng:
            def __init__(self,seed):self.original=np.random.Generator(np.random.PCG64(seed));self.bit_generator=self.original.bit_generator
            def uniform(self,*a,**kw):return self.original.uniform(*a,**kw)
        rng=FixedRng(74018)
        # Replace stochastic sampling and optimization with deterministic fixtures.
        rng.choice=lambda population,*,size,p:np.full(size,677380,dtype=int)
        probability=np.zeros(2**20,dtype=complex);probability[677380]=1
        def minimize(objective,initial,**kw):
            objective(initial)
            return SimpleNamespace(x=initial,success=False,nfev=1,message='offline fixture')
        with tempfile.TemporaryDirectory() as directory,patch.object(sys,'argv',['run_gates','--output',directory,'--reads','8']),patch.object(run_gates,'minimize',minimize),patch.object(run_gates,'qaoa_state',return_value=probability),patch.object(run_gates,'design_statevector',return_value=probability),patch.object(run_gates.np.random,'default_rng',return_value=rng):
            run_gates.main()
            paths=list(Path(directory).rglob('result.json'));self.assertEqual(len(paths),2)
            for path in paths:
                result=json.loads(path.read_text());inp=json.loads((path.parent/'input.json').read_text())
                self.assertEqual(result['num_samples'],8);self.assertEqual(result['metrics']['optimum_hits'],8)
                self.assertIn('sampling_rng_state',inp['settings'])
                self.assertEqual(len(result['probabilities']),2**20)
                self.assertEqual(result['optimization']['actual_evaluations'],1)
                self.assertGreaterEqual(result['timing']['wallclock_seconds'],result['timing']['optimization_seconds']+result['timing']['sampling_seconds'])
                self.assertEqual({p.name for p in path.parent.iterdir()},{'input.json','result.json'})

    def test_serialized_sample_counts_and_scores(self):
        scorer=Scorer()
        from run_io import check_result
        folders=list((ROOT/'results').rglob('result.json'))
        self.assertGreaterEqual(len(folders),6)
        for path in folders:
            result=json.loads(path.read_text())
            if result['status']!='COMPLETED':continue
            check_result(result);summary=result['metrics']
            counts={''.join(map(str,x)):c for x,c in zip(result['solutions'],result['num_occurrences'],strict=True)}
            self.assertEqual(sum(counts.values()),summary['reads'])
            actual=Counter();opt=feasible=0
            for row in result['scores']:
                count=row['occurrences'];actual[row['sample_bits']]+=count;y=np.array(list(map(int,row['design_bits'])))
                i=index_of(y)
                self.assertAlmostEqual(row['repaired_energy'],scorer.problem.repair(y),places=8)
                self.assertEqual(row['feasible'],bool(scorer.problem.feasible(y)))
                opt+=count*(i in scorer.optima);feasible+=count*bool(scorer.feasible[i])
            groups=result['selection_counts']
            self.assertEqual(len(groups),summary['distinct_feasible_selections'])
            self.assertEqual(sum(g['occurrences'] for g in groups),feasible)
            self.assertEqual(len({g['selection_bits'] for g in groups}),len(groups))
            self.assertEqual(summary['total_feasible_selections'],84)
            self.assertEqual(dict(actual),counts);self.assertEqual(opt,summary['optimum_hits']);self.assertEqual(feasible,summary['feasible_hits'])
            # Check actual embedded input coefficients independently of result summaries.
            inp=json.loads((path.parent/'input.json').read_text());model=inp['model']
            for bits,energy in zip(result['solutions'],result['energy'],strict=True):
                expected=model['offset']+sum(v*x for v,x in zip(model['linear'],bits,strict=True))+sum(v*bits[i]*bits[j] for i,j,v in model['quadratic'])
                self.assertAlmostEqual(expected,energy,places=8)
            if 'probabilities' in result:
                prob=np.array(result['probabilities']);self.assertAlmostEqual(float(prob.sum()),1,places=9)
                self.assertAlmostEqual(float(prob[scorer.feasible].sum()),summary['metadata']['exact_feasible_probability'],places=10)

if __name__=='__main__':unittest.main()
