#!/usr/bin/env python3
"""Local ILRS SA, one dated input/result pair per model and repeat."""
import argparse
import hashlib
from collections import Counter
import importlib.metadata
import inspect
import time
from pathlib import Path

import dimod
import neal
import numpy as np
from scoring import Scorer
from ilrs import ROOT,index_of
from run_io import (config,model_dict,create_run,read_json,result_template,save_result,
                    add_samples,check_result,timing,utc_now,fail_run,implementation)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model',type=Path,default=ROOT/'model')
    parser.add_argument('--output',type=Path,default=ROOT/'results',help='Results root; run names generated automatically')
    parser.add_argument('--seed',type=int)
    parser.add_argument('--repeats',type=int,default=1)
    args=parser.parse_args()
    if args.repeats<1:parser.error('--repeats must be positive')
    if args.seed is not None and not 0<=args.seed<=2**31-args.repeats:parser.error('Invalid seed/repeat range')
    setup=time.perf_counter();scorer=Scorer(args.model);p=scorer.problem
    models={'full':dimod.BinaryQuadraticModel.from_qubo(p.q+np.diag(p.linear),offset=p.offset),
            'surrogate':dimod.BinaryQuadraticModel.from_qubo(scorer.q+np.diag(scorer.linear),offset=scorer.offset)}
    provenance=dict(packages={k:importlib.metadata.version(k) for k in ['numpy','dimod','dwave-neal','dwave-samplers']},
                    model_summary=scorer.summary,implementation=implementation(__file__),
                    source_hashes={n:hashlib.sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in ['scoring.py','ilrs.py','selections.py']},
                    sampler_source=inspect.getsource(neal.SimulatedAnnealingSampler.sample))
    shared_seconds=time.perf_counter()-setup;shared_id=utc_now()
    for repeat in range(args.repeats):
        for name,model in models.items():
            kwargs={'num_reads':1000}
            if args.seed is not None:kwargs['seed']=args.seed+repeat
            started=time.perf_counter()
            inp=config('ilrs','full_qubo' if name=='full' else 'surrogate','sa',model_dict(model),kwargs,provenance=provenance)
            folder=create_run(args.output,inp)
            result=result_template(read_json(folder/'input.json'),'RUNNING')
            result['timestamps']['started_utc']=utc_now();save_result(folder,result)
            try:
                call=time.perf_counter();samples=neal.SimulatedAnnealingSampler().sample(model,**kwargs);seconds=time.perf_counter()-call
                post=time.perf_counter();counts=Counter()
                for datum in samples.data(['sample','energy','num_occurrences']):
                    x=[int(datum.sample[i]) for i in range(len(model.variables))]
                    expected=p.energy(np.array(x)) if name=='full' else scorer.surrogate[index_of(x)]
                    np.testing.assert_allclose(datum.energy,expected,atol=1e-8,rtol=0)
                    counts[''.join(map(str,x))]+=int(datum.num_occurrences)
                if sum(counts.values())!=1000:raise AssertionError('Read mismatch')
                rows,groups,summary=scorer.score_counts(counts,'SA_'+name,seconds,dict(requested_parameters=kwargs,sampler_info=samples.info),full=name=='full')
                add_samples(result,rows,summary)
                result.update(status='COMPLETED',selection_counts=groups,info={'sampler_info':samples.info},
                              execution_time=seconds,execution_time_basis='local sampler-call wallclock',
                              validation={'sample_energy_checked':True,'source_model_hashes_checked':True})
                result['timing']=timing(sampling_seconds=seconds,postprocessing_seconds=time.perf_counter()-post,
                                       wallclock_seconds=time.perf_counter()-started,
                                       shared_setup_seconds=shared_seconds,shared_setup_id=shared_id,
                                       measurement_notes='Run wallclock excludes shared batch setup and final JSON write; count shared setup once.')
                result['timestamps']['completed_utc']=utc_now()
                save_result(folder,check_result(result))
                print(folder)
            except Exception as error:
                fail_run(folder,started,error);raise


if __name__=='__main__':main()
