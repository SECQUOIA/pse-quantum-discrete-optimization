"""ILRS discrete model, coefficient reconstruction, and exact pair repair."""
from pathlib import Path
import csv
import json
import hashlib
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / 'inputs'
import sys
sys.path.insert(0, str(ROOT.parent/'scripts'))
FLOW_PAIRS = [('srce','r1'),('r1','s1'),('r1','s2'),('r1','s3'),
              ('srce','r2'),('r2','s1'),('r2','s2'),('r2','s3'),
              ('s1','sink'),('s2','sink'),('s3','sink')]
NAMES = ['y(r1)','y(r2)','y(s1)','y(s2)','y(s3)',
         'z_cat(1)','z_cat(2)','z_an(1)','z_an(2)'] + [f'f({a}_{b})' for a,b in FLOW_PAIRS]
PRODUCT_NAMES = ['w(1_1)','w(2_1)','w(1_2)','w(2_2)']


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')


def bits(index, n=20):
    return ((np.asarray(index, dtype=np.uint64)[..., None] >> np.arange(n-1,-1,-1,dtype=np.uint64)) & 1).astype(np.int8)


def index_of(y):
    return int(np.asarray(y, dtype=np.int64) @ (1 << np.arange(len(y)-1,-1,-1)))


class Problem:
    def __init__(self):
        provenance=json.loads((INPUT/'provenance.json').read_text())
        for name,expected in provenance['sha256'].items():
            if hashlib.sha256((INPUT/name).read_bytes()).hexdigest()!=expected:
                raise ValueError(f'Archived input changed: {name}')
        base=INPUT/'discrete_qubo/julia_exports'
        self.n,self.scale,self.offset=np.loadtxt(base/'scalars_1.csv',delimiter=',',skiprows=1)
        self.n=int(self.n)
        self.linear=np.loadtxt(base/'L_vector_1.csv',delimiter=',')
        self.q=np.loadtxt(base/'Q_matrix_1.csv',delimiter=',')
        self.c=self.q+self.q.T
        if self.n!=40 or self.scale!=1 or np.any(np.diag(self.q)):
            raise ValueError('Unexpected exported convention')
        with (INPUT/'data/parameter_cost.csv').open() as f:
            costs={r['unit'][0]+str(int(r['unit'][1:])+1):r for r in csv.DictReader(f)}
        with (INPUT/'data/parameter_alpha.csv').open() as f:
            alpha={int(r['reactor']):float(r['alpha']) for r in csv.DictReader(f)}
        self.unit_cost=np.array([float(costs[u]['fixed']) + (2*float(costs[u]['operating'])*alpha[int(u[1:])] if u.startswith('r') else 0) for u in ['r1','r2','s1','s2','s3']])
        self.product_cost=np.zeros((3,4))
        with (INPUT/'data/parameter_beta.csv').open() as f:
            for r in csv.DictReader(f):
                s,c,a=map(int,[r['separator-k'],r['cation-c'],r['anion-a']])
                self.product_cost[s-1,(a-1)*2+c-1]=2*float(costs[f's{s}']['operating'])*float(r['beta'])
        self.components=[]
        unseen=set(range(20,40))
        while unseen:
            group=[];stack=[min(unseen)];unseen.remove(stack[0])
            while stack:
                i=stack.pop();group.append(i)
                neighbors={j for j in unseen if self.c[i,j]!=0};unseen-=neighbors;stack.extend(sorted(neighbors))
            self.components.append(sorted(group))
        if [len(c) for c in self.components]!=[2]*10:
            raise ValueError('Expected ten repair pairs')
        self.reconstruct()

    def energy(self,x):
        return self.scale*(self.offset + x@self.linear + np.einsum('...i,ij,...j->...',x,self.q,x,optimize=False))

    def reconstruct(self):
        """Independently rebuild the 40-bit polynomial, verifying index semantics."""
        l=np.zeros(40);q=np.zeros((40,40));offset=0.;rho=63.68
        def square(coeff, constant=0):
            nonlocal offset
            offset+=rho*constant**2
            for i,v in coeff.items(): l[i]+=rho*(v*v+2*constant*v)
            for i,v in coeff.items():
                for j,w in coeff.items():
                    if i<j:q[i,j]+=2*rho*v*w
        l[:5]+=self.unit_cost
        for s in range(3):
            for w in range(4):q[2+s,20+w]+=self.product_cost[s,w]
        # Unit selection/entry-exit consistency.
        for i,j in [(0,9),(1,13),(2,17),(3,18),(4,19)]:square({i:1,j:-1})
        square({17:1,18:1,19:1,24:-1,25:-1},-1)
        square({9:-1,10:1,11:1,12:1,26:-1,27:-2})
        square({13:-1,14:1,15:1,16:1,28:-1,29:-2})
        for flow1,flow2,out,s1,s2 in [(10,14,17,30,31),(11,15,18,32,33),(12,16,19,34,35)]:
            square({flow1:1,flow2:1,out:-1,s1:-1,s2:-1})
        square({5:1,6:1},-1);square({7:1,8:1},-1)
        # Polynomial penalties for AND products, with Rosenberg auxiliaries.
        for w,c,a,t in [(20,5,7,39),(21,6,7,36),(22,5,8,37),(23,6,8,38)]:
            l[w]+=rho;l[t]+=4*rho;q[w,t]-=2*rho
            q[c,a]+=rho;q[c,t]-=2*rho;q[a,t]-=2*rho
        # At least one source reactor, f1+f2-f1*f2=1, and flow implications.
        # Binary penalties reduce to (1-f1)*(1-f2) and edge*(1-entry/exit).
        offset+=rho;l[9]-=rho;l[13]-=rho;q[9,13]+=rho
        for edge,entry,out in [(10,9,17),(11,9,18),(12,9,19),(14,13,17),(15,13,18),(16,13,19)]:
            l[edge]+=2*rho;q[min(edge,entry),max(edge,entry)]-=rho;q[edge,out]-=rho
        np.testing.assert_allclose(l,self.linear,atol=1e-10,rtol=0)
        np.testing.assert_allclose(q,self.q,atol=1e-10,rtol=0)
        np.testing.assert_allclose(offset,self.offset,atol=1e-10,rtol=0)

    def objective(self,y):
        products=np.stack([y[...,5]*y[...,7],y[...,6]*y[...,7],y[...,5]*y[...,8],y[...,6]*y[...,8]],axis=-1)
        return y[...,:5]@self.unit_cost + np.einsum('...s,sw,...w->...', y[...,2:5],self.product_cost,products)

    def feasible(self,y):
        f=y[...,9:];ok=(y[...,5]+y[...,6]==1)&(y[...,7]+y[...,8]==1)
        ok &= (f[...,0]==y[...,0])&(f[...,4]==y[...,1])
        for s in range(3):ok &= f[...,8+s]==y[...,2+s]
        ok &= (f[...,0]+f[...,4]>=1)&(np.sum(f[...,8:11],axis=-1)>=1)
        for r in range(2):
            edges=f[...,r*4+1:r*4+4];entry=f[...,r*4]
            ok &= np.sum(edges,axis=-1)>=entry
            ok &= np.all(edges<=entry[...,None],axis=-1)
        for s in range(3):
            edges=np.stack([f[...,1+s],f[...,5+s]],axis=-1);out=f[...,8+s]
            ok &= np.sum(edges,axis=-1)>=out
            ok &= np.all(edges<=out[...,None],axis=-1)
        return ok

    def repair(self,y,return_bits=False):
        y=np.asarray(y);single=y.ndim==1
        if single:y=y[None,:]
        energy=self.offset+y@self.linear[:20]+np.einsum('bi,ij,bj->b',y,self.q[:20,:20],y)
        full=np.zeros((len(y),40),dtype=np.int8);full[:,:20]=y
        for group in self.components:
            a=bits(np.arange(4),2)
            local=a@self.linear[group]+np.einsum('bi,ij,bj->b',a,self.q[np.ix_(group,group)],a)
            values=(y@self.c[np.ix_(range(20),group)])@a.T+local
            choice=np.argmin(values,axis=1);energy+=values[np.arange(len(y)),choice];full[:,group]=a[choice]
        if single:energy=energy[0];full=full[0]
        return (energy,full) if return_bits else energy
