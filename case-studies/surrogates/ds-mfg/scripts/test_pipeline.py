"""Regression tests for conversion conventions and saved sample recoverability."""
import csv
import json
import tempfile
import unittest
from pathlib import Path

import dimod
import numpy as np

from run_sa import Model, ROOT, Scorer, assignments, load_models


class PipelineTests(unittest.TestCase):
    def test_nonunit_scale_diagonal_and_both_triangles(self):
        # A different convention than the real upper-triangular export catches
        # accidental factor-of-two and diagonal/offset handling mistakes.
        model = Model({'scalars.csv': b'n,scale,offset\n2,3,7\n',
                       'L_vector.csv': b'2\n-4\n',
                       'Q_matrix.csv': b'5,6\n8,9\n'})
        x = assignments(2)
        expected = 3 * (7 + 7*x[:, 0] + 5*x[:, 1] + 14*x[:, 0]*x[:, 1])
        np.testing.assert_allclose(model.bqm.energies((x, [0, 1])), expected)
        np.testing.assert_allclose(model.energies(x), expected)

    def test_saved_pilot_samples_reproduce_scores_and_counts(self):
        full, surrogate = load_models()
        scorer = Scorer(full, surrogate)
        folders=list((ROOT/'results/sa').glob('*/result.json'))
        self.assertGreaterEqual(len(folders),2)
        for path in folders:
            result=json.loads(path.read_text())
            if result['status']!='COMPLETED':continue
            summary=result['metrics']
            model=full if result['formulation']=='full_qubo' else surrogate
            expected={}
            native=result['info'].get('native_sampleset')
            if native:
                samples=dimod.SampleSet.from_serializable(native)
                for datum in samples.data(['sample','energy','num_occurrences']):
                    bits=[int(datum.sample[i]) for i in range(model.n)]
                    key=''.join(map(str,bits));expected[key]=expected.get(key,0)+int(datum.num_occurrences)
                    self.assertAlmostEqual(float(datum.energy),model.energies(np.array([bits]))[0],places=8)
            else:
                for bits,energy,count in zip(result['solutions'],result['energy'],result['num_occurrences'],strict=True):
                    key=''.join(map(str,bits));expected[key]=expected.get(key,0)+count
                    self.assertAlmostEqual(energy,model.energies(np.array([bits]))[0],places=8)
            actual={};opt=feasible=0
            for row in result['scores']:
                key=row['sample_bits'];actual[key]=actual.get(key,0)+row['occurrences']
                score=scorer.score(list(map(int,key)))
                self.assertAlmostEqual(row['repaired_energy'],score['repaired_energy'],places=8)
                self.assertEqual(row['feasible'],score['feasible'])
                if result['formulation']=='full_qubo':self.assertGreaterEqual(row['auxiliary_gap'],-1e-8)
                opt+=row['occurrences']*score['optimum_hit'];feasible+=row['occurrences']*score['feasible']
            self.assertEqual(actual,expected)
            self.assertEqual(sum(actual.values()),result['num_samples'])
            self.assertEqual((opt,feasible),(summary['optimum_hits'],summary['feasible_hits']))


if __name__ == '__main__':
    unittest.main()
