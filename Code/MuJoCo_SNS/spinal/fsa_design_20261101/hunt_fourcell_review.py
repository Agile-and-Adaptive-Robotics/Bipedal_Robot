"""Gait neural pathway composition reviewer: current Hunt four-cell audit.

Reuses existing composition evidence; only reviews its invariant PF/FF/bi
buffers and extracted four-cell operators. No body or training calculations.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;sys.path[:0]=[str(HERE),str(HERE.parent),str(HERE.parent/'w2l_cpg')]
from hunt_pathway_audit import C,MODEL,options_from_physical_condition,mujoco
import ibexc_coincidence as IB


def fourcell_names(side):
    return [f'PF_BI_JOINT_{f}_{side}' for f in ('HAM','RECT')]+[f'PF_IN_BI_PAIR_{f}_{side}' for f in ('HAM','RECT')]


def operator(n,side,use_original_template):
    from sns_toolbox.networks import Network
    from sns_toolbox.neurons import NonSpikingNeuron
    from sns_toolbox.connections import NonSpikingSynapse
    c=n.compiled;core=c.c;cells=fourcell_names(side);signals=[f'PF_{joint}-{phase}_{side}' for joint in ('HIP','KNEE') for phase in ('E','F')];all_names=signals+cells;lookup={name:i for i,name in enumerate(n.names)};local={name:i for i,name in enumerate(all_names)};edges=[]
    for i,(post,pre) in enumerate(zip(c.post,c.pre)):
        a,b=n.names[pre],n.names[post]
        if b not in cells:continue
        if a not in all_names:
            if c.gmax[i]>1e-12:raise ValueError('Four-cell block has unsupported active incoming source')
            continue
        g=float(n.original_g[i] if use_original_template else c.gmax[i]);edges.append({'source':a,'target':b,'g':g,'E':float(c.reversal[i]),'lo':float(core.e_lo[post,pre]),'hi':float(core.e_hi[post,pre])})
    cases=[]
    for pattern in ('off','RF','HAM'):
        tiny=Network(name='extracted four-cell RF HAM operator')
        for name in all_names:
            j=lookup[name]
            if core.g_ion[:,j].any():raise ValueError('Expected passive PF/relay cells,not an NaP oscillator')
            tiny.add_neuron(NonSpikingNeuron(membrane_capacitance=float(core.c_m[j]),membrane_conductance=float(core.g_m[j]),resting_potential=float(core.V_rest[j]),bias=float(core.i_b[j])),name=name)
        for e in edges:
            if e['g']>0:tiny.add_connection(NonSpikingSynapse(max_conductance=e['g'],reversal_potential=e['E'],e_lo=e['lo'],e_hi=e['hi']),e['source'],e['target'])
        t=tiny.compile(dt=core.dt,backend='numpy');drive={'RF':{f'PF_HIP-F_{side}',f'PF_KNEE-E_{side}'},'HAM':{f'PF_HIP-E_{side}',f'PF_KNEE-F_{side}'}}.get(pattern,set())
        for step in range(400):
            for name in signals:t.V[local[name]]=5. if name in drive else 0.
            t.forward([])
        ham=float(t.V[local[cells[0]]]);rf=float(t.V[local[cells[1]]]);cases.append({'pattern':pattern,'HAM_mV':ham,'RF_mV':rf,'relay_mV':{name:float(t.V[local[name]]) for name in cells[2:]},'finite':bool(np.isfinite(t.V).all())})
    if not all(r['finite'] for r in cases):raise AssertionError('Extracted four-cell response is nonfinite')
    if not (cases[1]['RF_mV']>cases[1]['HAM_mV']+.1 and cases[2]['HAM_mV']>cases[2]['RF_mV']+.1):raise AssertionError('Four-cell operator does not select requested antagonist')
    return {'side':side,'parameters':'original allocated source template; separate expression fixture' if use_original_template else 'current fitted actual command1 conductances','edges':edges,'cases':cases,'not_intrinsic_rhythm_or_body_claim':True}


def review(composition_path,out):
    proof=json.loads(composition_path.read_text());condition=Path(proof['condition']);cfg=json.loads(condition.read_text());model=mujoco.MjModel.from_xml_path(str(MODEL));n=C.build([model.actuator(i).name for i in range(model.nu)],dt=.002,options=options_from_physical_condition(cfg));excluded_PT=n.__dict__.pop('balance',None) is not None;n.apply_weights(1.)
    identity=IB.hunt_identity(n);pairing=IB.audit_hunt_gain_pairing(n)
    if len(n.compiled.V)!=1033:raise ValueError('Primary design must be current1033 Hunt/PT,not historical1027')
    for buffer in ('pf_gain','ff_gain','bi_gain'):
        expected=proof['source_gain_pairing_before_profiles']['buffers'][buffer]['full_buffer_sha256'];composed=proof['composed_gain_pairing']['buffers'][buffer]['full_buffer_sha256']
        if pairing['buffers'][buffer]['full_buffer_sha256']!=expected or expected!=composed:raise ValueError('Review buffers differ from actual composition proof')
    c=n.compiled;names=n.names;cells=sum((fourcell_names(s) for s in ('r','l')),[])
    if any(name not in names for name in cells):raise ValueError('Missing four-cell source cell')
    info=[];links=[];motor=[]
    for name in cells:
        i=n.idx[name];info.append({'name':name,'index':int(i),'capacitance_backend_microF':float(c.c.c_m[i]),'leak_uS':float(c.c.g_m[i]),'bias_nA':float(c.c.i_b[i]),'gated':bool(c.c.g_ion[:,i].any())})
    for i,(post,pre) in enumerate(zip(c.post,c.pre)):
        a,b=names[pre],names[post]
        if b in cells:links.append({'source':a,'target':b,'actual_g_uS':float(c.gmax[i]),'template_g_uS':float(n.original_g[i]),'reversal_mV':float(c.reversal[i]),'authoritative_bi_fit_parameter':bool(i in n.bi_edges),'above_tolerance':bool(c.gmax[i]>1e-12)})
        if a in cells and b.startswith('MN_'):motor.append({'slot':int(i),'source':a,'target':b,'authoritative_ff_gain_uS':float(n.ff_gain[i]),'command1_g_uS':float(c.gmax[i]),'template_g_uS':float(n.original_g[i]),'above_tolerance':bool(c.gmax[i]>1e-12),'role':'Fitted command-gated PF-to-MN parameter,not injected torque'})
    command=[]
    for gate in (0.,.5,1.):
        n.apply_weights(gate);command.append({'command':gate,'allocated_motor_slots':len(motor),'above_tolerance_motor_slots':sum(c.gmax[r['slot']]>1e-12 for r in motor),'maximum_motor_conductance_uS':max((float(c.gmax[r['slot']]) for r in motor),default=0.)})
    current=[operator(n,s,False) for s in ('r','l')];template=[operator(n,s,True) for s in ('r','l')]
    for name in ('pf_gain','ff_gain','bi_gain'):
        if hashlib.sha256(getattr(n,name).tobytes()).hexdigest()!=pairing['buffers'][name]['full_buffer_sha256']:raise AssertionError('Four-cell review modified fitted parameters')
    active=[r for r in motor if r['above_tolerance']];row={'primary_source':'Current1033-neuron Hunt/PT,already four-cell RF/HAM per leg','source_condition':str(condition),'composition_proof':str(composition_path),'fourcell_cells':cells,'current_MN_output_status':'Command-enabled coefficients active above tolerance' if active else 'Allocated fit variables only; no meaningful current RF/HAM-to-MN expression','meaningful_motor_slots':len(active),'maximum_ff_gain_uS':max((r['authoritative_ff_gain_uS'] for r in motor),default=0.),'command_and_cache_composition_verified':True,'physical_biarticular_contribution_proved':False,'historical1027':'Untied comparison only; not promoted as primary'}
    result={'role':'Gait neural pathway composition reviewer','identity':identity,'family_review_rows':{'Gait2392_biw2lvar_Hunt_PT':row,'Gait2392_fourcell_RF_HAM':row},'comparison_tolerance_uS':1e-12,'separate_SAR_GRAC_TFL_groups':{k:v for k,v in n.extra_pairs.items() if k.startswith(('BI_SAR','BI_GRAC','BI_TFL'))},'composition_source_sha256':proof['source_sha256'],'composition_proof_sha256':hashlib.sha256(composition_path.read_bytes()).hexdigest(),'review_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'cell_parameters':info,'incoming_and_reciprocal_edges':links,'motor_outputs':motor,'command_policy':command,'actual_fitted_operator':current,'separate_template_operator':template,'source_gain_pairing':pairing,'unbound_PT_excluded':excluded_PT,'mechanics_unchanged':True,'physics_steps':0,'training':False,'user_approval_required':True,'global_release':False}
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2,default=lambda v:int(v) if isinstance(v,np.integer) else bool(v) if isinstance(v,np.bool_) else str(v)),encoding='utf-8');print(json.dumps({'neurons':1033,'cells':len(cells),'motor_slots':len(motor),'meaningful_motor_slots':len(active),'max_ff_uS':row['maximum_ff_gain_uS'],'output':str(out)}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--composition',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();review(a.composition,a.out)
