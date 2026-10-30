"""Shared bilateral conductances: one parameter per mirrored edge pair."""
import re
import numpy as np

def canonical(label):
    return re.sub(r'_[rl](?=$|[^a-zA-Z0-9])','_r',label)

def tied_solutions(matrix,rhs,lookup,names,pre,post,solver,upper):
    groups={};weights={};fits={}
    for target in matrix:groups.setdefault(canonical(names[target]),[]).append(target)
    for targets in groups.values():
        keys=sorted({canonical(names[pre[i]]+' -> '+names[post[i]]) for t in targets for i in lookup[t]})
        blocks=[];ys=[]
        for t in targets:
            raw=matrix[t][:len(rhs[t]):5];y=rhs[t][::5];X=np.zeros((len(y),len(keys)))
            for j,i in enumerate(lookup[t]):X[:,keys.index(canonical(names[pre[i]]+' -> '+names[post[i]]))]+=raw[:,j]
            blocks.append(X);ys.append(y)
        X=np.vstack(blocks);y=np.concatenate(ys);regularizer=.03*np.sqrt(len(y))
        solution=solver(np.vstack([X,regularizer*np.eye(len(keys))]),np.r_[y,np.zeros(len(keys))],bounds=(0.,upper),method='bvls',tol=1e-8,max_iter=200)
        values=dict(zip(keys,solution.x))
        for t,block,desired in zip(targets,blocks,ys):
            for i in lookup[t]:
                edge=names[pre[i]]+' -> '+names[post[i]];weights[edge]=float(values[canonical(edge)])
            fits[names[t]]=dict(current_rmse=float(np.sqrt(np.mean((block@solution.x-desired)**2))),shared_parameters=len(keys))
    return weights,fits

def assert_tied(edges):
    seen={}
    for key,value in edges.items():
        c=canonical(key)
        if c in seen and abs(seen[c]-value)>1e-12:raise ValueError('Bilateral conductance mismatch: '+key)
        seen[c]=value
