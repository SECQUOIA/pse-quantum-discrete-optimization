#!/usr/bin/env python3
"""Ideal local QAOA/VQE pilots using exact surrogate energy objectives."""
import argparse
import hashlib
import importlib.metadata
import time
from pathlib import Path
from collections import Counter

import numpy as np
from scipy.optimize import minimize
from qiskit.circuit.library import efficient_su2
from qiskit.quantum_info import Statevector
from scoring import Scorer
from ilrs import ROOT
from run_io import (config,create_run,read_json,result_template,save_result,add_samples,
                    check_result,timing,utc_now,implementation,fail_run)


def qaoa_state(energy,angles):
    n=int(np.log2(len(energy)));state=np.full(len(energy),1/np.sqrt(len(energy)),dtype=complex)
    p=len(angles)//2
    # Angles are beta then gamma. exp(-i beta X) on each logical variable.
    for beta,gamma in zip(angles[:p],angles[p:]):
        state*=np.exp(-1j*gamma*energy)
        for i in range(n):
            blocks=state.reshape(-1,2,1<<i);a=blocks[:,0,:].copy();b=blocks[:,1,:].copy()
            blocks[:,0,:]=np.cos(beta)*a-1j*np.sin(beta)*b
            blocks[:,1,:]=np.cos(beta)*b-1j*np.sin(beta)*a
    return state


def design_statevector(circuit,theta):
    """Qiskit qubit i is design variable i; design strings are big-endian."""
    raw=Statevector.from_instruction(circuit.assign_parameters(theta)).data
    n=circuit.num_qubits
    return raw.reshape([2]*n).transpose(tuple(reversed(range(n)))).reshape(-1)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model',type=Path,default=ROOT/'model')
    parser.add_argument('--output',type=Path,default=ROOT/'results',help='Results root; dated runs generated automatically')
    parser.add_argument('--algorithm',choices=['qaoa','vqe','both'],default='both')
    parser.add_argument('--seed',type=int,default=74018)
    parser.add_argument('--reads',type=int,default=4096)
    parser.add_argument('--qaoa-p',type=int,default=2)
    parser.add_argument('--qaoa-evaluations',type=int,default=25)
    parser.add_argument('--vqe-reps',type=int,default=1)
    parser.add_argument('--vqe-evaluations',type=int,default=100)
    args=parser.parse_args()
    if args.reads<1 or args.qaoa_p<1 or args.vqe_reps<1:parser.error('Reads/depths must be positive')
    if args.qaoa_evaluations<2*args.qaoa_p+2:parser.error('QAOA evaluation budget must cover COBYLA initialization')
    if args.vqe_evaluations<40*(args.vqe_reps+1)+2:parser.error('VQE evaluation budget must cover COBYLA initialization')
    setup_started=time.perf_counter()
    scorer=Scorer(args.model)
    energy=np.asarray(scorer.surrogate);mean=float(energy.mean());scale=float(energy.std());normalized=(energy-mean)/scale
    rng=np.random.default_rng(args.seed);n=20
    shared_setup_seconds=time.perf_counter()-setup_started;shared_setup_id=utc_now()
    sampling_model=dict(vartype='BINARY',variables=list(range(n)),scale=1.0,offset=float(scorer.offset),
                        linear=scorer.linear.tolist(),
                        quadratic=[[i,j,float(scorer.q[i,j]+scorer.q[j,i])] for i in range(n) for j in range(i+1,n) if scorer.q[i,j]+scorer.q[j,i]!=0],
                        convention='offset + linear[i]*x[i] + quadratic[i,j]*x[i]*x[j], i<j')
    for algorithm in ['qaoa','vqe']:
        if args.algorithm not in ['both',algorithm]:continue
        run_started=time.perf_counter();history=[]
        if algorithm=='qaoa':
            initial=rng.uniform(-0.5,0.5,2*args.qaoa_p);budget=args.qaoa_evaluations
            circuit=None
            def state_at(theta):return qaoa_state(normalized,theta)
            circuit_metadata={'p':args.qaoa_p,'cost':'exp(-i gamma normalized_surrogate_energy)',
                              'mixer':'exp(-i beta sum X)','initial_state':'uniform',
                              'angle_order':'beta_then_gamma','implementation':'exact diagonal phase and X rotations'}
        else:
            circuit=efficient_su2(n,su2_gates=['ry','rz'],entanglement='reverse_linear',reps=args.vqe_reps)
            initial=rng.uniform(-np.pi,np.pi,circuit.num_parameters);budget=args.vqe_evaluations
            # Design variable 0 is the leftmost bit/index most significant bit.
            # Circuit qubit 0 represents design variable 0; reverse state axes
            # to reconcile Qiskit's little-endian statevector convention.
            def state_at(theta):
                return design_statevector(circuit,theta)
            circuit_metadata={'ansatz':'EfficientSU2','reps':args.vqe_reps,'su2_gates':['ry','rz'],
                              'entanglement':'reverse_linear','parameters':circuit.num_parameters,
                              'circuit_depth':circuit.depth(),'gate_counts':dict(circuit.count_ops())}
        def objective(theta):
            state=state_at(theta);prob=np.abs(state)**2
            value=float(prob@normalized)
            history.append({'evaluation':len(history)+1,'normalized_expected_energy':value,
                            'expected_surrogate_energy':value*scale+mean})
            return value
        settings=dict(seed=args.seed,reads=args.reads,sampling_rng_state=rng.bit_generator.state,initial_parameters=initial.tolist(),
                      optimizer='COBYLA',max_evaluations=budget,rhobeg=0.3,tolerance=1e-3,
                      circuit=circuit_metadata,objective_shift=mean,objective_scale=scale)
        inp=config('ilrs','surrogate',algorithm,sampling_model,settings,
                   provenance={'implementation':implementation(__file__),
                               'source_hashes':{n:hashlib.sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in ['scoring.py','ilrs.py','selections.py']},
                               'packages':{k:importlib.metadata.version(k) for k in ['numpy','scipy','qiskit']},
                               'model_summary':scorer.summary},execution_context='ideal local statevector')
        folder=create_run(args.output,inp)
        saved=result_template(read_json(folder/'input.json'),'RUNNING')
        saved['timestamps']['started_utc']=utc_now();save_result(folder,saved)
        try:
            started=time.perf_counter()
            result=minimize(objective,initial,method='COBYLA',options={'maxiter':budget,'rhobeg':0.3,'tol':1e-3})
            optimization_seconds=time.perf_counter()-started
            started=time.perf_counter();state=state_at(result.x);prob=np.abs(state)**2
            if abs(prob.sum()-1)>1e-8:raise AssertionError('Statevector normalization failed')
            indices=rng.choice(len(prob),size=args.reads,p=prob/prob.sum());sampling_seconds=time.perf_counter()-started
            post=time.perf_counter()
            counts=Counter(format(int(i),'020b') for i in indices)
            metadata=dict(backend='ideal local statevector',objective='exact expected surrogate energy',
                          objective_shift=mean,objective_scale=scale,seed=args.seed,
                          optimization_seconds=optimization_seconds,circuit=circuit_metadata,
                          optimizer_converged=bool(result.success),optimizer_evaluations=int(result.nfev),
                          exact_optimum_probability=float(sum(prob[i] for i in scorer.optima)),
                          exact_feasible_probability=float(prob[scorer.feasible].sum()),
                          expected_surrogate_energy=float(prob@energy))
            rows,groups,summary=scorer.score_counts(counts,algorithm.upper()+'_surrogate',sampling_seconds,metadata)
            add_samples(saved,rows,summary)
            saved.update(status='COMPLETED',selection_counts=groups,probabilities=prob.tolist(),
                         optimization={'optimized':result.x.tolist(),'converged':bool(result.success),
                                       'message':str(result.message),'actual_evaluations':int(result.nfev),'history':history},
                         info=metadata,validation={'probability_normalization_checked':True,'source_model_hashes_checked':True},
                         execution_time=optimization_seconds+sampling_seconds,
                         execution_time_basis='local parameter optimization plus final sampling wallclock')
            saved['timing']=timing(optimization_seconds=optimization_seconds,sampling_seconds=sampling_seconds,
                                   postprocessing_seconds=time.perf_counter()-post,
                                   wallclock_seconds=time.perf_counter()-run_started,
                                   shared_setup_seconds=shared_setup_seconds,shared_setup_id=shared_setup_id,
                                   measurement_notes='Per-algorithm wallclock excludes shared setup and final JSON write; shared setup counted once.')
            saved['timestamps']['completed_utc']=utc_now();save_result(folder,check_result(saved));print(folder,flush=True)
        except Exception as error:fail_run(folder,run_started,error);raise


if __name__=='__main__':main()
