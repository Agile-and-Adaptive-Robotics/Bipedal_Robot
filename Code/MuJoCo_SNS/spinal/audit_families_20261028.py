"""Compiled inventories and exact-fast parity for 92- and 54-actuator families."""
import os,sys,json,csv,copy,time,hashlib,argparse
from pathlib import Path
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent);os.environ.pop('AARL_NET',None)
import numpy as np
import build_network as bn
import params
HERE=Path(__file__).resolve().parent
OUT=HERE/'reports_20261028_family_audit';OUT.mkdir(exist_ok=True)
sys.path.insert(0,str(HERE/'w2l_cpg'))
from exact_fast_eval import ExactFastForward

def main(family,body):
    if body=='robot54':
        import xml.etree.ElementTree as ET
        root=ET.parse(HERE.parent/'robot_torso/ah001_abduction_diagnostic_v2.xml').getroot()
        acts=[a.get('name') for a in root.find('actuator')]
        params.G.update(c1_gain=.15,v3_gain=0.,heel_rge=.15,toe_rge=.05)
    else:
        data=np.load(HERE/'synergy_basis.npz',allow_pickle=True)
        acts=sorted([str(n) for s in ('r','l') for n in data['muscle_names_'+s]]+
                    ['quad_fem_r','quad_fem_l','gem_r','gem_l','peri_r','peri_l'])
        assert len(acts)==92 and len(set(acts))==92
    if family=='s3k':net=bn.build(acts,dt=.002,interleg=True);module=bn
    elif family=='w2lvar':
        import build_network_w2lvar as module
        net=module.build(acts,dt=.002,interleg=True)
    else:
        import build_network_syn6 as module
        net=module.build(acts,dt=.002,interleg=True)
    names=[p['name'] for p in net.net.populations];edges=[]
    for c in net.net.connections:
        p=c['params'];edges.append(dict(source=names[c['source']],target=names[c['destination']],
            g=float(p['max_conductance']),reversal=float(p['reversal_potential'])))
    key=family+'_'+body
    with (OUT/(key+'_edges.csv')).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(edges[0]));writer.writeheader();writer.writerows(edges)
    pf_mn=[e for e in edges if e['source'].startswith('PF_') and not e['source'].startswith('PF_IN') and e['target'].startswith('MN_') and e['g']>0]
    covered={e['target'][3:] for e in pf_mn}
    pairs={}
    for e in edges:pairs.setdefault((e['source'],e['target']),[]).append(e)
    duplicate=[v for v in pairs.values() if len(v)>1]
    details={}
    if family=='w2lvar':
        weights=bn._jpf_weights();reintroduced=[];present=[]
        for act in acts:
            base=act.rsplit('_',1)[0];mi=module.classify(act)
            for cell in module.W2L_HCS:
                vals=[weights.get(row,{}).get(base) for row in module.W2L_CELL_ROWS[cell]]
                if any(v is not None for v in vals):
                    present.append([act,cell,vals])
                    if sum(v for v in vals if v is not None)<=0 and net._w2l_pf_weight(mi,base,cell,weights)>0:reintroduced.append([act,cell,vals])
        details=dict(explicit_fitted_zero_reintroduced=reintroduced,table_entries_present=len(present),
                     variant_gains=dict(module.GAINS))
    if family=='syn6':
        tables=module._load_synergy_tables()
        details=dict(families=net.syn_families,mixed=net.syn_mixed,
                     basis_sha256=hashlib.sha256(module.SYNERGY_NPZ.read_bytes()).hexdigest(),
                     basis_covered_actuators=sorted(a for a in acts if a in tables[a[-1]]['names']),
                     dress=dict(module.DRESS))
    # Use matched cloned state and seeded inputs; never regroup reductions.
    ref=copy.deepcopy(net.compiled);fast=ExactFastForward(copy.deepcopy(net.compiled))
    rng=np.random.default_rng(20261028);inputs=[]
    for i in range(200):
        x=net.make_inputs()
        for j,port in enumerate(net.inputs):x[j]=3. if port=='DRIVE' else .15 if port=='POSTURE' else rng.uniform(0,.3)
        inputs.append(x)
    maximum=0.;gate_maximum=0.
    for x in inputs:
        ref.forward(x);fast.forward(x)
        maximum=max(maximum,float(np.max(abs(ref.V-fast.V))))
        if ref.gated:gate_maximum=max(gate_maximum,float(np.max(abs(ref.b_gate-fast.b_gate))))
    assert maximum==gate_maximum==0.
    # Warm JIT outside timing, then use the same input sequence and initial state.
    a=copy.deepcopy(net.compiled);b=ExactFastForward(copy.deepcopy(net.compiled));b.forward(inputs[0])
    times={}
    for label,compiled in [('reference',a),('exact_fast',b)]:
        start=time.perf_counter()
        for x in inputs:compiled.forward(x)
        times[label]=(time.perf_counter()-start)/len(inputs)
    result=dict(family=family,body=body,actuators=acts,true_neurons=int(net.compiled.V.size),
        named_indices=len(net.idx),inputs=len(net.inputs),builder_edges=len(edges),
        compiled_nonzero_edges=int(np.count_nonzero(net.compiled.g_max_non)),
        pf_cells=[n for n in names if n.startswith('PF_') and not n.startswith('PF_IN')],
        pf_to_mn_edges=len(pf_mn),actuators_without_pf=sorted(set(acts)-covered),
        renshaw_cells=[n for n in names if n.startswith('RC_')],duplicate_endpoint_edges=duplicate,
        direct_rg_inhibition=[e for e in edges if e['source'].startswith('RG_') and e['target'].startswith('RG_') and e['reversal']<0],
        direct_pf_inhibition=[e for e in edges if e['source'] in [p['source'] for p in pf_mn] and e['target'] in [p['source'] for p in pf_mn] and e['reversal']<0],
        per_muscle_central_edges=[e for e in edges if e['source'].startswith(('Ia_','II_','Ib_')) and e['target'].startswith(('RG_','PF_','InE_','InF_'))],
        gains=dict(params.G),time_constants=dict(params.TAU),details=details,
        parity=dict(steps=len(inputs),voltage_max_difference=maximum,h_gate_max_difference=gate_maximum,bit_exact=True),
        timing_seconds_per_neural_step=times,neural_speedup=times['reference']/times['exact_fast'],
        source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(module.__file__),HERE/'params.py',HERE/'w2l_cpg/exact_fast_eval.py']})
    (OUT/(key+'.json')).write_text(json.dumps(result,indent=2))
    print(json.dumps({k:result[k] for k in ('family','body','true_neurons','inputs','builder_edges','compiled_nonzero_edges','actuators_without_pf','neural_speedup')}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('family',choices=['s3k','w2lvar','syn6']);p.add_argument('body',choices=['human92','robot54']);a=p.parse_args();main(a.family,a.body)
