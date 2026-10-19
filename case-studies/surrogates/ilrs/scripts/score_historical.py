#!/usr/bin/env python3
"""Import and rescore historical samples; original computation time stays distinct."""
import argparse
from collections import Counter
import hashlib
from pathlib import Path
import time

import dimod
import numpy as np
from scoring import Scorer
from ilrs import ROOT,INPUT
from run_io import (config,model_dict,create_run,read_json,result_template,save_result,
                    add_samples,check_result,timing,utc_now,fail_run)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'results/historical',help='Historical results root')
    args=parser.parse_args();setup=time.perf_counter();scorer=Scorer();p=scorer.problem
    model=dimod.BinaryQuadraticModel.from_qubo(p.q+np.diag(p.linear),offset=p.offset)
    shared_seconds=time.perf_counter()-setup;shared_id=utc_now()
    for source in sorted((INPUT/'discrete_qubo/result_raw').glob('*.json')):
        data=read_json(source)
        if not data.get('solutions'):continue
        solver='qa-unknown' if data['solver'].startswith('dwave') else 'sa'
        started=time.perf_counter()
        inp=config('ilrs','full_qubo',solver,model_dict(model),{},
                   provenance={'source_file':source.name,'sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
                               'original_solver':data['solver'],
                               'source_hashes':{n:hashlib.sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in ['scoring.py','ilrs.py','selections.py','score_historical.py']}},execution_context='historical rescoring; no sampler executed',
                   original_run_date='2026-01-14T16:35:01' if source.stem=='SA_results_20260114_163501' else None,
                   original_timezone=None)
        # The dated directory describes this import execution; the historical execution date remains explicit.
        folder=create_run(args.output,inp);result=result_template(read_json(folder/'input.json'),'RUNNING')
        result['timestamps']['started_utc']=utc_now()
        try:
            variables=data['variables']
            if sorted(map(int,variables))!=list(range(40)):raise ValueError('Historical variable order invalid')
            order=[variables.index(str(i)) for i in range(40)]
            counts=Counter();max_error=0.
            for solution,energy,count in zip(data['solutions'],data['energy'],data['num_occurrences'],strict=True):
                x=np.asarray(solution)[order]
                if not np.isin(x,[0,1]).all():raise ValueError('Nonbinary historical sample')
                error=abs(float(p.energy(x))-float(energy));max_error=max(error,max_error)
                if error>1e-8:raise AssertionError('Historical energy mismatch')
                counts[''.join(str(int(i)) for i in x)]+=int(count)
            if sum(counts.values())!=data['num_samples']:raise ValueError('Historical read mismatch')
            metadata=dict(source_file=source.name,solver=data['solver'],historical=True,original_info=data.get('info',{}),
                          source_sample_energy_max_error=max_error,execution_status='saved samples rescored')
            rows,groups,summary=scorer.score_counts(counts,data['solver']+'_historical',data['execution_time'],metadata,full=True)
            add_samples(result,rows,summary)
            result.update(status='COMPLETED',selection_counts=groups,
                          execution_time=data['execution_time'],execution_time_basis='original historical duration; boundary unverified',
                          info={'original_info':data.get('info',{}),'original_meta':data.get('meta',{})},
                          validation={'source_sample_energy_max_error':max_error})
            elapsed=time.perf_counter()-started
            result['timing']=timing(wallclock_seconds=elapsed,postprocessing_seconds=elapsed,
                                   shared_setup_seconds=shared_seconds,shared_setup_id=shared_id,
                                   measurement_notes='Wallclock measures this import/rescoring, not the historical solver execution.')
            result['timestamps']['completed_utc']=utc_now();save_result(folder,check_result(result));print(folder)
        except Exception as error:fail_run(folder,started,error);raise


if __name__=='__main__':main()
