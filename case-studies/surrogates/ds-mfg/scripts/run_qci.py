#!/usr/bin/env python3
"""Prepare a DS-MFG full-QUBO/surrogate QCI run, submit once, or collect by job ID."""
import argparse
import importlib.metadata
import json
import os
import shutil
from pathlib import Path
import sys
import time

import numpy as np
from run_sa import ROOT, Scorer, assignments, digest, load_models, validate
from run_io import (config, model_dict, create_run, read_json, result_template, save_result,
                    timing, utc_now, begin_submission, add_samples, check_result, implementation, run_lock)


def prepare(output, device='dirac-3', num_samples=None, formulation='surrogate'):
    if device not in ('dirac-1','dirac-3'):
        raise ValueError('Unsupported device')
    if formulation not in ('full_qubo','surrogate'):
        raise ValueError('Unsupported formulation')
    num_samples = (100 if device=='dirac-1' else 5) if num_samples is None else num_samples
    if not 1 <= num_samples <= 100:
        raise ValueError('Pilot sample count must be in [1, 100]')
    started = time.perf_counter()
    full, surrogate = load_models()
    validation = validate(full, surrogate, Scorer(full, surrogate))
    model = full if formulation=='full_qubo' else surrogate
    if formulation=='surrogate':
        x = assignments(model.n)
        scope = 'All binary states'
    else:
        # Zero, singleton and pair states determine every degree-two coefficient.
        # Add all-ones and fixed random probes without enumerating 2**36 states.
        singletons = np.eye(model.n,dtype=np.int8)
        pairs = np.array([singletons[i]+singletons[j] for i in range(model.n) for j in range(i+1,model.n)])
        x = np.vstack([np.zeros((1,model.n),dtype=np.int8),singletons,pairs,
                       np.ones((1,model.n),dtype=np.int8),
                       np.random.default_rng(74018).integers(0,2,size=(100,model.n),dtype=np.int8)])
        scope = 'Complete degree-two coefficient basis, all-ones, and 100 fixed random probes'
    file_name = f'ds_mfg_{formulation}_{model.n}'
    settings = {'device_type':device,'num_samples':num_samples}
    if device=='dirac-1':
        from qci_client.optimization.client import data_to_json
        matrix = np.diag([float(model.bqm.linear[i]) for i in range(model.n)])
        for (i,j),value in model.bqm.quadratic.items():
            matrix[i,j] = matrix[j,i] = float(value)/2
        qubo = {'file_name':file_name+'_dirac1','file_config':{
            'qubo':{'num_variables':model.n,'data':matrix.tolist()}}}
        wire_file = data_to_json(file={'file_name':qubo['file_name'],
                                      'file_config':{'qubo':{'data':matrix}}})
        # Validate the actual sparse SDK representation, not only the dense matrix.
        energy = np.zeros(len(x))
        for term in wire_file['file_config']['qubo']['data']:
            energy += term['val']*x[:,term['i']]*x[:,term['j']]
        submission = {'job_type':'sample-qubo','qubo':qubo,'wire_file':wire_file}
        validation_key = 'qci_qubo_maximum_error'
        solver = 'qci-dirac1'
    else:
        terms = [{'idx':[0,i+1],'val':float(model.bqm.linear[i])} for i in range(model.n)]
        for (i,j),value in model.bqm.quadratic.items():
            i,j = sorted((i,j))
            terms.append({'idx':[i+1,j+1],'val':float(value)})
        energy = np.zeros(len(x))
        for term in terms:
            energy += term['val']*np.prod([x[:,i-1] for i in term['idx'] if i],axis=0)
        polynomial = {'file_name':file_name,'file_config':{'polynomial':{
            'num_variables':model.n,'min_degree':1,'max_degree':2,'data':terms}}}
        submission = {'job_type':'sample-hamiltonian-integer','polynomial':polynomial}
        settings.update(relaxation_schedule=1,num_levels=[2]*model.n)
        validation_key = 'qci_polynomial_maximum_error'
        solver = 'qci-dirac3'
    error = float(np.max(np.abs(energy+model.bqm.offset-model.energies(x))))
    if error > 1e-8:
        raise AssertionError('QCI export energy mismatch')
    validation[validation_key] = error
    validation.update(qci_export_states_checked=len(x),qci_export_validation_scope=scope)
    inp = config('ds-mfg',formulation,solver,model_dict(model.bqm),settings,
                 provenance={'implementation':implementation(__file__),
                             'qci_client_version':importlib.metadata.version('qci-client'),
                             'scorer_sha256':digest(Path(__file__).with_name('run_sa.py').read_bytes()),
                             'full_coefficient_hashes':full.hashes,'surrogate_coefficient_hashes':surrogate.hashes},
                 submission=submission)
    path = create_run(output,inp)
    result = read_json(path/'result.json')
    result['validation'] = validation
    result['timing'] = timing(setup_seconds=time.perf_counter()-started,
                             wallclock_seconds=time.perf_counter()-started,
                             measurement_notes='Accumulated runner-invocation wallclock; excludes gaps between invocations and final JSON writes.')
    save_result(path,result)
    print('Prepared:',path)
    return path


def client():
    from qci_client import QciClient
    return QciClient(url='https://api.qci-prod.com',api_token=os.environ['QCI_TOKEN'],timeout=30)


def add_local_time(result, started):
    old = result['timing'].get('wallclock_seconds')
    result['timing']['wallclock_seconds'] = (old or 0) + time.perf_counter()-started
    if old is None:
        result['timing']['measurement_notes'] = 'Partial local wallclock: earlier invocations were not timed; gaps between invocations excluded.'


def submit(path, api=None):
    started = time.perf_counter()
    inp = read_json(path/'input.json')
    if inp['solver'] not in ('qci-dirac1','qci-dirac3'):
        raise ValueError('Not a QCI run')
    result = begin_submission(path)
    job_requested = False
    try:
        api = api or client()
        allocation = api.get_allocations()['allocations']['dirac']
        if allocation['metered'] and allocation['seconds'] <= 0:
            raise RuntimeError('No available Dirac allocation')
        if inp['solver']=='qci-dirac1':
            from qci_client.optimization.client import data_to_json
            payload = inp['submission']['qubo']
            payload = {'file_name':payload['file_name'],'file_config':{'qubo':{'data':np.asarray(payload['file_config']['qubo']['data'],dtype=float)}}}
            if data_to_json(file=payload) != inp['submission']['wire_file']:
                raise ValueError('SDK serialization differs from validated upload')
            file_key = 'qubo_file_id'
        else:
            payload = inp['submission']['polynomial']
            file_key = 'polynomial_file_id'
        uploaded = api.upload_file(file=payload)
        result['info'][file_key] = uploaded['file_id']
        save_result(path,result)
        body = api.build_job_body(job_type=inp['submission']['job_type'],job_params=inp['settings'],
                                  **{file_key:uploaded['file_id']},job_name=f"ds-mfg-{inp['formulation']}-{inp['settings']['device_type']}-{inp['run_id']}")
        result['info']['job_submission'] = body
        result['timestamps']['submitted_utc'] = utc_now()
        save_result(path,result)
        job_requested = True
        submitted = api.submit_job(job_body=body)
        result.update(job_id=submitted['job_id'],status='SUBMITTED')
        print('Submitted:',submitted['job_id'])
    except Exception as error:
        # A timeout may hide a successful submission. Never automatically submit again.
        # qci-client 4.5 wraps HTTPError and leaves response on the original cause.
        cause = error
        response = None
        while cause is not None and response is None:
            response = getattr(cause,'response',None)
            cause = cause.__cause__
        code = response.status_code if response is not None else None
        result['status'] = 'FAILED' if not job_requested else 'REJECTED' if code is not None and 400<=code<500 else 'SUBMISSION_UNCERTAIN'
        result['info']['failure_type'] = type(error).__name__
        if code is not None:
            result['info']['http_status'] = code
            # Keep only provider error body, never request headers or credentials.
            message = response.text[:2000]
            token = os.environ.get('QCI_TOKEN')
            result['info']['provider_error'] = message.replace(token,'[REDACTED]') if token else message
        # Known no-job failures are not experiments and do not belong in results.
        # Preserve ambiguous requests: a timed-out POST may already have created a job.
        error.qci_status = result['status']
        error.qci_http_status = result['info'].get('http_status')
        error.qci_provider_error = result['info'].get('provider_error')
        if result['status'] in ('FAILED','REJECTED') and not result.get('job_id'):
            shutil.rmtree(path)
        raise
    finally:
        if path.exists():
            add_local_time(result,started)
            save_result(path,result)


def _timestamp(value):
    from datetime import datetime
    return datetime.fromisoformat(value.replace('Z','+00:00'))


def provider_timing(result, response):
    status = response['job_info'].get('job_status',{})
    result['timestamps'].update(status)
    def duration(a,b):
        if not status.get(a) or not status.get(b):return None
        return max(0.0,(_timestamp(status[b])-_timestamp(status[a])).total_seconds())
    result['timing']['queue_seconds'] = duration('queued_at_rfc3339nano','running_at_rfc3339nano')
    result['timing']['remote_turnaround_seconds'] = duration('submitted_at_rfc3339nano','completed_at_rfc3339nano')
    result['timing']['device_seconds'] = response['job_info'].get('job_result',{}).get('device_usage_s')


def retrieve_metrics(api, result):
    """Get detailed timing for an existing job; timing failure never triggers resubmission."""
    if (result['info'].get('provider_metrics') or {}).get('job_metrics'):
        return
    method = getattr(api, 'get_job_metrics', None)
    if method is None:
        result['info']['provider_metrics_status'] = 'not supported by client'
        return
    try:
        response = method(job_id=result['job_id'])
        result['info']['provider_metrics'] = response
        result['info']['provider_metrics_status'] = 'available' if response.get('job_metrics') else 'not returned by provider'
        result['info'].pop('provider_metrics_failure_type', None)
    except Exception as error:
        result['info']['provider_metrics_status'] = 'unavailable; collection can be retried'
        result['info']['provider_metrics_failure_type'] = type(error).__name__


def collect(path, api=None):
    # Serialize collectors so an older QUEUED response cannot overwrite completed samples.
    with run_lock(path):
        return _collect(path, api)


def _collect(path, api=None):
    started = time.perf_counter()
    inp = read_json(path/'input.json')
    result = read_json(path/'result.json')
    if not result.get('job_id'):
        raise RuntimeError('No saved job ID; reconcile submission before collection')
    if result['status'] == 'COMPLETED':
        if not (result['info'].get('provider_metrics') or {}).get('job_metrics'):
            retrieve_metrics(api or client(), result)
            add_local_time(result, started)
        save_result(path, result)
        print(json.dumps({'run':str(path),'status':'COMPLETED','metrics':result['metrics']}))
        return result
    api = api or client()
    try:
        response = api.get_job_results(job_id=result['job_id'])
        result['info']['provider_response'] = response
        provider_timing(result,response)
        result['status'] = response['status']
        if response['status'] == 'COMPLETED':
            retrieve_metrics(api, result)
            result['status'] = 'COLLECTING'
            # Save the provider evidence before local scoring; collection remains retryable.
            save_result(path,result)
            post = time.perf_counter()
            full, surrogate = load_models()
            provenance = inp['provenance']
            if full.hashes != provenance['full_coefficient_hashes'] or surrogate.hashes != provenance.get('surrogate_coefficient_hashes',provenance.get('coefficient_hashes')):
                raise RuntimeError('Scoring coefficients changed since preparation')
            model = full if inp['formulation']=='full_qubo' else surrogate
            if model_dict(model.bqm) != inp['model']:
                raise RuntimeError('Sampling model changed since preparation')
            scorer = Scorer(full,surrogate)
            rows = []
            data = response['results']
            for bits, energy, count in zip(data['solutions'],data['energies'],data['counts'],strict=True):
                if len(bits)!=model.n or any(v not in (0,1) for v in bits) or int(count)!=count or count<=0:
                    raise ValueError('Invalid returned sample')
                bits=list(map(int,bits));local=float(model.energies(np.array([bits]))[0])
                rows.append(dict(sample_bits=''.join(map(str,bits)),occurrences=int(count),
                                 sampled_model_energy=local,device_energy=float(energy),
                                 device_energy_residual=float(energy)+inp['model']['offset']-local,
                                 **scorer.score(bits)))
                if inp['formulation']=='full_qubo':
                    rows[-1]['auxiliary_gap'] = local-rows[-1]['repaired_energy']
            if sum(r['occurrences'] for r in rows)!=inp['settings']['num_samples']:
                raise ValueError('Sample count mismatch')
            feasible=[r for r in rows if r['feasible']]
            summary=dict(optimum_hits=sum(r['occurrences'] for r in rows if r['optimum_hit']),
                         feasible_hits=sum(r['occurrences'] for r in feasible),
                         unique_feasible_designs=len({r['flow_bits'] for r in feasible}),
                         best_feasible_cost=min((r['feasible_cost'] for r in feasible),default=None),
                         maximum_device_energy_residual=max(abs(r['device_energy_residual']) for r in rows))
            add_samples(result,rows,summary)
            result['execution_time']=result['timing']['device_seconds']
            result['execution_time_basis']='provider device usage seconds; not local wallclock'
            result['timing']['postprocessing_seconds']=time.perf_counter()-post
            result['timestamps'].setdefault('first_result_retrieved_utc',utc_now())
            result['status']='COMPLETED'
            check_result(result)
        print(json.dumps({'run':str(path),'status':result['status'],'metrics':result['metrics']}))
    except Exception as error:
        result['info']['collection_failure_type']=type(error).__name__
        if result['status']=='COLLECTING':result['status']='COLLECTION_FAILED'
        raise
    finally:
        add_local_time(result,started)
        save_result(path,result)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'results',help='Results root for preparation')
    parser.add_argument('--device',choices=['dirac-1','dirac-3'],default='dirac-3',help='Device for preparation; existing runs use saved settings')
    parser.add_argument('--formulation',choices=['full_qubo','surrogate'],default='surrogate',help='Model for preparation; existing runs use saved input')
    parser.add_argument('--num-samples',type=int,help='Pilot sample count: default 100 for Dirac-1, 5 for Dirac-3')
    parser.add_argument('--run',type=Path,help='Existing generated run directory for submission/collection')
    modes=parser.add_mutually_exclusive_group()
    modes.add_argument('--submit',action='store_true')
    modes.add_argument('--collect',action='store_true')
    args=parser.parse_args()
    if args.submit or args.collect:
        if args.run is None:parser.error('--run is required for submission/collection')
        (submit if args.submit else collect)(args.run)
    else:
        if args.run is not None:parser.error('--run is only used with --submit or --collect')
        prepare(args.output,args.device,args.num_samples,args.formulation)


if __name__=='__main__':
    main()
