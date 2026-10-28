"""Fresh-process, literature-motivated probes on the existing Gait2392 model.

No physics, production winner, or database edits. Each exploratory trial saves
the effective gain dictionary and compiled edge evidence, including knobs not
handled by historical curriculum loaders.
"""
import os,sys,json,time,hashlib,argparse
from pathlib import Path
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
HERE=Path(__file__).resolve().parent;os.chdir(HERE)
OUT=HERE/'reports_20261028_family_audit/human_trials';OUT.mkdir(parents=True,exist_ok=True)
p=argparse.ArgumentParser();p.add_argument('family',choices=['s3k','w2lvar','syn6']);p.add_argument('case');p.add_argument('--reference',action='store_true');p.add_argument('--preallocated',action='store_true');a=p.parse_args()
os.environ['AARL_NET']='' if a.family=='s3k' else a.family
base=json.loads((HERE/'best_walk_params_v10.json').read_text())['multipliers']
if a.family=='s3k':
    import _curriculum as C
    base.update(json.loads((HERE/'reports_20260923/s3k_trial34_full_params.json').read_text())['params']);base.update(renshaw=.5,full_rules=0.)
    C.set_stage(4,base)
else:
    if a.family=='syn6':import _curriculum_syn6 as C
    else:import _curriculum_w2lvar as C
    seed=HERE/('curriculum_syn6_stage5_reseed_20261003.json' if a.family=='syn6' else 'curriculum_w2lvar_stage5.json')
    base.update(json.loads(seed.read_text())['params']);C.set_stage(5,base)
import params as P
import runner as R
import build_network as BN
sys.path.insert(0,str(HERE/'w2l_cpg'));from exact_fast_eval import ExactFastForward
updates={
 's3k':{'baseline':{},'rg_pf_gate':{},'ib_half':{'ib_e_central':P.G['ib_e_central']*.5},'afferent_half':{'ia_f_central':P.G['ia_f_central']*.5,'ii_f_central':P.G['ii_f_central']*.5}},
 'syn6':{'baseline':{},'central05':{'ib_e_central':.05},'central15':{'ib_e_central':.15},'pf_gate':{'iain_pf_gate':.1,'pf_mixed_recip':.1},'central_gate':{'ib_e_central':.05,'iain_pf_gate':.1,'pf_mixed_recip':.1}},
 'w2lvar':{'baseline':{},'central_half':{},'comm_half':{},'central_quarter':{}}
}[a.family][a.case]
P.G.update(updates);base.update(updates)
capture={};original=BN.build
def build(acts,dt,interleg):
    net=original(acts,dt=dt,interleg=interleg)
    if a.family=='s3k' and a.case=='rg_pf_gate':
        # Shinohara Table A.1: opposite RG inhibitory IN -> PF phase.
        # Add to the actual compiled matrices only, preserving existing edges.
        for side in ('r','l'):
            for pre,phases,g in [('InF',('E1','E2'),.21),('InE',('F1','F2'),.03)]:
                for phase in phases:
                    post=net.idx['PF_'+phase+'_'+side];source=net.idx[pre+'_'+side]
                    assert net.compiled.g_max_non[post,source]==0
                    net.compiled.g_max_non[post,source]=g;net.compiled.del_e[post,source]=BN.E_REV_INH
                    net.compiled.e_lo[post,source]=0.;net.compiled.e_hi[post,source]=P.E_HI
    capture.update(true_neurons=int(net.compiled.V.size),inputs=len(net.inputs),nonzero_edges=int(__import__('numpy').count_nonzero(net.compiled.g_max_non)),effective_gains=dict(P.G),effective_tau=dict(P.TAU))
    if a.preallocated:
        from preallocated_exact_eval import PreallocatedExactForward
        net.compiled=PreallocatedExactForward(net.compiled)
    elif not a.reference:net.compiled=ExactFastForward(net.compiled)
    return net
if a.family=='w2lvar':
    import build_network_w2lvar as V
    if a.case in ('central_half','central_quarter'):
        scale=.5 if a.case=='central_half' else .25
        for k in ('w2l_ib_central','w2l_ia_central','w2l_ii_central'):V.GAINS[k]*=scale
    if a.case=='comm_half':
        # w2l_comm is a boolean topology gate, not a conductance multiplier.
        for k in ('w2l_v2a','w2l_v0v','w2l_ini','w2l_ini_inh','w2l_v0d','w2l_v0d_inh','w2l_v3e','w2l_v3e_contra','w2l_v3e_inE'):V.GAINS[k]*=.5
    capture['variant_gains']=dict(V.GAINS)
BN.build=build
label=a.family+'_'+a.case+('_reference' if a.reference else '_preallocated' if a.preallocated else '')
os.environ['AARL_NPZ']=str(OUT/(label+'.npz'))
start=time.perf_counter();metrics=R.main(['--eval','--drive',repr(base['drive'])])
result=dict(family=a.family,case=a.case,machine='easteregg2',backend='reference' if a.reference else 'preallocated_exact' if a.preallocated else 'exact_fast',
    wall_s=time.perf_counter()-start,metrics=metrics,configuration=base,probe_updates=updates,compiled=capture,
    source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),HERE/'runner.py',HERE/'build_network.py',HERE/'build_network_w2lvar.py',HERE/'build_network_syn6.py',HERE/'params.py',HERE/'w2l_cpg/exact_fast_eval.py',HERE/'w2l_cpg/preallocated_exact_eval.py']},
    scope='Exploratory supported Gait2392 evaluation, not an accepted hardware or CAD-robot gait result')
(OUT/(label+'.json')).write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
