#!/usr/bin/env python3
"""Enumerate repaired ILRS landscape and fit its uniform quadratic projection."""
import argparse
import csv
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from ilrs import ROOT, INPUT, NAMES, PRODUCT_NAMES, Problem, bits, index_of, write_json
from selections import selection_catalog, save_catalog


def walsh_transform(values):
    values=np.asarray(values,dtype=float)
    coefficients=values.copy()
    h=1
    while h<len(values):
        blocks=coefficients.reshape(-1,2*h)
        left=blocks[:,:h].copy();right=blocks[:,h:].copy()
        blocks[:,:h]=left+right;blocks[:,h:]=left-right;h*=2
    return coefficients/len(values)


def quadratic_projection(values,n=20):
    coeff=walsh_transform(values)
    first=np.array([coeff[1<<(n-1-i)] for i in range(n)])
    pairs=np.zeros((n,n))
    for i in range(n):
        for j in range(i+1,n):pairs[i,j]=coeff[(1<<(n-1-i))|(1<<(n-1-j))]
    offset=coeff[0]+first.sum()+pairs.sum()
    linear=-2*first-2*(pairs.sum(axis=0)+pairs.sum(axis=1))
    return offset,linear,4*pairs


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'model')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    started=time.perf_counter();p=Problem()
    energies=np.empty(2**20);surrogate=np.empty_like(energies)
    feasible=np.empty(2**20,dtype=bool);objectives=np.empty_like(energies)
    for start in range(0,len(energies),8192):
        stop=start+8192;y=bits(np.arange(start,stop))
        energies[start:stop]=p.repair(y);feasible[start:stop]=p.feasible(y);objectives[start:stop]=p.objective(y)
    reference=energies[feasible]-objectives[feasible]
    np.testing.assert_allclose(reference,0,atol=1e-8,rtol=0)
    optimum=float(objectives[feasible].min())
    opt_indices=np.flatnonzero(feasible & np.isclose(objectives,optimum,atol=1e-9,rtol=0))
    feasible_indices=np.flatnonzero(feasible)
    # Compare all recorded pool assignments with independently enumerated model.
    pool=json.loads((INPUT/'discrete_ip/result_gurobi/ilrs_gpool_n228_20260126_101442.json').read_text())
    pool_indices=[]
    for solution in pool['solutions']:
        values=solution['values'];y=np.array([int(round(values[name])) for name in NAMES])
        index=index_of(y);pool_indices.append(index)
        if not feasible[index]:raise AssertionError('Pool design violates original constraints')
        np.testing.assert_allclose(objectives[index],solution['ip_obj_value'],atol=1e-8,rtol=0)
        products=[y[5]*y[7],y[6]*y[7],y[5]*y[8],y[6]*y[8]]
        if products!=[int(round(values[k])) for k in PRODUCT_NAMES]:raise AssertionError('Pool product map mismatch')
    if energies.min()<optimum-1e-8:raise AssertionError('Repaired penalty minimum undercuts true optimum')
    selections=selection_catalog(feasible,objectives)
    selection_definition=save_catalog(args.output,selections)
    enumeration_seconds=time.perf_counter()-started
    offset,linear,q=quadratic_projection(energies)
    for start in range(0,len(energies),8192):
        y=bits(np.arange(start,start+8192))
        surrogate[start:start+8192]=offset+y@linear+np.einsum('bi,ij,bj->b',y,q,y)
    error=surrogate-energies
    # Uniform projection residual must be orthogonal to constant, linear, pair features.
    residual_coeff=walsh_transform(error)
    degree2=[0]+[1<<(19-i) for i in range(20)]+[(1<<(19-i))|(1<<(19-j)) for i in range(20) for j in range(i+1,20)]
    orthogonal_error=float(np.max(np.abs(residual_coeff[degree2])))
    if orthogonal_error>1e-8:raise AssertionError('Fit is not uniform least-squares projection')
    # Full brute-force auxiliary checks for two independently chosen design vectors.
    brute_errors=[]
    aux=bits(np.arange(2**20))
    for index in [int(opt_indices[0]),0]:
        best=float('inf');y=bits(index)
        for start in range(0,len(aux),8192):
            full=np.column_stack([np.tile(y,(8192,1)),aux[start:start+8192]])
            best=min(best,float(p.energy(full).min()))
        brute_errors.append(abs(best-energies[index]))
    if max(brute_errors)>1e-8:raise AssertionError('Pair repair disagrees with brute force')
    for name,array in [('repaired_energies',energies),('surrogate_energies',surrogate),('feasible',feasible),('objectives',objectives)]:np.save(args.output/f'{name}.npy',array)
    np.savetxt(args.output/'linear.csv',linear,delimiter=',');np.savetxt(args.output/'quadratic.csv',q,delimiter=',')
    with (args.output/'feasible_designs.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['design_bits','objective','optimum'])
        for index in sorted(feasible_indices,key=lambda i:(objectives[i],i)):
            writer.writerow([''.join(map(str,bits(index))),objectives[index],bool(index in opt_indices)])
    mapping=[dict(qubo_index=i,variable=name,role='design' if i<20 else 'derived_product') for i,name in enumerate(NAMES+PRODUCT_NAMES)]
    mapping.extend(dict(qubo_index=i,variable=f'aux_{i}',role='reformulation_auxiliary') for i in range(24,40))
    variable_map=dict(index_base=0,variables=mapping,design_order=NAMES,repair_components=p.components,
               proof='Independent coefficient-by-coefficient reconstruction from original objective, logic penalties, slack squares, and Rosenberg product reduction')
    min_index=int(np.argmin(surrogate));top=np.argsort(energies)[:50]
    summary=dict(design_variables=20,discrete_variables=24,full_qubo_variables=40,eliminated_variables=20,
        auxiliary_components=p.components,total_design_assignments=len(energies),feasible_count=int(feasible.sum()),feasible_flow_assignment_count=int(feasible.sum()),
        feasible_selection_count=len(selections),selection_variables=9,
        optimum_objective=optimum,optimum_indices=opt_indices.tolist(),
        optimum_design_bits=[''.join(map(str,bits(i))) for i in opt_indices],
        recorded_pool_count=len(pool_indices),recorded_pool_selection_count=len({i >> 11 for i in pool_indices}),missing_from_recorded_pool=len(set(feasible_indices)-set(pool_indices)),
        rmse=float(np.sqrt(np.mean(error**2))),maximum_absolute_error=float(np.max(np.abs(error))),
        r_squared=float(1-np.var(error)/np.var(energies)),fit_orthogonality_error=orthogonal_error,
        optimum_surrogate_ranks=[int(1+np.sum(surrogate<surrogate[i]-1e-8)) for i in opt_indices],
        surrogate_minimum_design_bits=''.join(map(str,bits(min_index))),surrogate_minimum_energy=float(surrogate[min_index]),
        surrogate_minimum_repaired_energy=float(energies[min_index]),surrogate_minimum_feasible=bool(feasible[min_index]),
        top50_maximum_absolute_error=float(np.max(np.abs(error[top]))),offset=float(offset),
        brute_force_repair_maximum_error=max(brute_errors),enumeration_seconds=enumeration_seconds,
        total_construction_validation_seconds=time.perf_counter()-started,
        coefficients_sha256={name:hashlib.sha256((args.output/name).read_bytes()).hexdigest() for name in ['linear.csv','quadratic.csv']},
        model_artifacts_sha256={name:hashlib.sha256((args.output/name).read_bytes()).hexdigest() for name in ['repaired_energies.npy','surrogate_energies.npy','feasible.npy','objectives.npy','feasible_selections.csv']},
        input_provenance_sha256=hashlib.sha256((INPUT/'provenance.json').read_bytes()).hexdigest())
    model_metadata=dict(summary=summary,variable_map=variable_map,selection_definition=selection_definition)
    model_metadata['metadata_sections_sha256']={k:hashlib.sha256(json.dumps(model_metadata[k],sort_keys=True).encode()).hexdigest() for k in ['variable_map','selection_definition']}
    write_json(args.output/'model.json',model_metadata);print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
