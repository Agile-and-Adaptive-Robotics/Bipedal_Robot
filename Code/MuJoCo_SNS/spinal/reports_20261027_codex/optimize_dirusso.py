"""Parallel PSO recovery experiment; every named dimension is consumed.

This uses an explicitly labeled engineering survival/speed objective,
not the paper's metabolic model or a claimed full equation-13 replication.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import concurrent.futures
import json
from pathlib import Path
import sys
import time
import numpy as np
import dirusso_recovery as R
from controller_numba import GAIN_KEYS

HERE=Path(__file__).resolve().parent
NAMES=[f'W_{m}_{j}' for m in R.C.MUSCLES for j in range(5)]
LO=[-1.]*45; HI=[1.]*45
for key in GAIN_KEYS:
    for muscle in R.C.MUSCLES:
        NAMES.append(key+'_'+muscle);LO.append(0.);HI.append(2.)
for name,lo,hi in (('omega',3.,8.),('gamma',0.,2.),
                   ('initial_vx',.0,1.8),('initial_height_offset',-.08,.02),
                   ('initial_pitch_offset',-.15,.15),('initial_vz_offset',-.2,.2)):
    NAMES.append(name);LO.append(lo);HI.append(hi)
for key,lo,hi in (('kp',-30.,30.),('kd',-5.,5.),('target',-.15,.15)):
    for muscle in R.C.BALANCE_MUSCLES:
        NAMES.append('balance_'+key+'_'+muscle);LO.append(lo);HI.append(hi)
for side in 'rl':
    for joint in ('hip_flexion','knee_angle','ankle_angle'):
        for prefix,bound in (('initial_q_',.15),('initial_v_',.5)):
            NAMES.append(prefix+joint+'_'+side);LO.append(-bound);HI.append(bound)
LO=np.array(LO);HI=np.array(HI)


def unpack(x):
    values=LO+x*(HI-LO)
    p=dict(zip(NAMES,values.tolist()))
    w=np.array([[p.pop(f'W_{m}_{j}') for j in range(5)] for m in R.C.MUSCLES])
    p['W']=np.concatenate([w,w]).tolist()
    return p


def evaluate(job):
    x,support,duration=job
    p=unpack(np.array(x))
    try:
        result=R.simulate(p,support=support,duration=duration)
    except (FloatingPointError,ValueError,RuntimeError) as e:
        result=dict(cost=1e6,error=str(e),survive=False,elapsed=0.)
    return result


def main():
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--resume')
    parser.add_argument('--prefix',default='dirusso_pso')
    parser.add_argument('--generations',type=int,default=12)
    parser.add_argument('--support',type=float,default=.9)
    parser.add_argument('--debug-surrogate',action='store_true',
                        help='Reproduce the rejected engineering experiment; not a gait-replication search')
    args=parser.parse_args()
    if not args.debug_surrogate:
        raise SystemExit('This survival/speed surrogate was rejected after selecting excessive body rotation. '
                         'Establish the original SCONE solution before further optimization. '
                         'See REPORT.md and ZCODE_HANDOFF.md.')
    prior=json.loads(Path(args.resume).read_text())['params'] if args.resume else {}
    rng=np.random.default_rng(20261027)
    W,_=R.imitation_weights()
    seed=np.zeros(len(NAMES))
    defaults=R.C.default_gains()
    for i,name in enumerate(NAMES):
        if name.startswith('W_'):
            _,muscle,j=name.split('_')
            value=np.array(prior.get('W',W))[R.C.MIDX[muscle],int(j)]
        elif name.startswith('balance_'):
            _,key,muscle=name.split('_')
            value=prior.get(name,prior.get('balance_'+key,dict(kp=8.,kd=1.5,target=0.)[key]))
        elif name in ('omega','gamma','balance_kp','balance_kd','balance_target',
                      'initial_vx','initial_height_offset'):
            value=dict(omega=2*np.pi/1.23,gamma=1.,balance_kp=8.,balance_kd=1.5,
                       balance_target=0.,initial_vx=1.175,initial_height_offset=-.03)[name]
        elif name.startswith('initial_'):
            value=0.
        else:
            key,muscle=name.rsplit('_',1)
            value=float(defaults[key][R.C.MIDX[muscle]])*.1
        if name in prior:
            value=prior[name]
        seed[i]=np.clip((value-LO[i])/(HI[i]-LO[i]),0.,1.)
    positions=np.clip(seed+rng.normal(0,.025 if prior else .06,(80,len(seed))),0.,1.)
    positions[0]=seed
    velocities=np.zeros_like(positions)
    personal=positions.copy();scores=np.full(80,np.inf)
    best=seed.copy();bestcost=np.inf
    support=args.support
    started=time.perf_counter()
    records=[]
    with concurrent.futures.ProcessPoolExecutor(max_workers=8) as pool:
        for gen in range(args.generations):
            results=list(pool.map(evaluate,[(x.tolist(),support,10.) for x in positions]))
            for i,result in enumerate(results):
                result.update(generation=gen,particle=i,position=positions[i].tolist())
                records.append(result)
                cost=float(result['cost'])
                if cost<scores[i]:
                    personal[i]=positions[i];scores[i]=cost
                if cost<bestcost:
                    bestcost=cost;best=positions[i].copy()
                    (HERE/(args.prefix+'_best.json')).write_text(json.dumps(dict(
                        params=unpack(best),metrics=result,support=support,
                        objective='survival/speed engineering surrogate'),indent=2))
            with (HERE/(args.prefix+'_trials.jsonl')).open('a') as stream:
                for result in results:
                    stream.write(json.dumps(result)+'\n')
            print('generation',gen,'best cost',round(bestcost,4),
                  'longest',max(x.get('elapsed',0.) for x in results),
                  'survivors',sum(x['survive'] for x in results),
                  'wall',round(time.perf_counter()-started,1),flush=True)
            velocities=.72*velocities+1.49*rng.random(positions.shape)*(personal-positions)+1.49*rng.random(positions.shape)*(best-positions)
            positions=np.clip(positions+np.clip(velocities,-.15,.15),0.,1.)
            np.savez_compressed(HERE/(args.prefix+'_checkpoint.npz'),
                                positions=positions,velocities=velocities,
                                personal=personal,scores=scores,best=best,
                                generation=gen,bestcost=bestcost)
    params=unpack(best)
    result=R.simulate(params,support=support,duration=10.,save=args.prefix+'_replay')
    print('REPLAY',json.dumps(result),flush=True)


if __name__=='__main__':
    main()
