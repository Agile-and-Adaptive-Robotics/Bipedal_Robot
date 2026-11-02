"""Static actual-Hunt sparse-step current accounting; no body integration."""
from pathlib import Path
import argparse,copy,hashlib,json,sys
from types import SimpleNamespace
import xml.etree.ElementTree as ET
import numpy as np

H=Path(__file__).resolve().parent;S=H.parent
sys.path[:0]=[str(H),str(S),str(S/'w2l_cpg'),str(S/'functional_gait2392')]
import controller as C
from verification_capture import sparse_current_snapshot,AlignedCapture,requested_cells


def fingerprint(network):
    h=hashlib.sha256()
    for scope,obj in(('wrapper',network.compiled),('native',network.compiled.c)):
        for name,value in sorted(vars(obj).items()):
            if isinstance(value,np.ndarray):
                h.update((scope+name+str(value.dtype)+str(value.shape)).encode());h.update(value.tobytes())
    return h.hexdigest()


def clone_wrapper(wrapper):
    result=object.__new__(type(wrapper));result.__dict__.update(copy.deepcopy(wrapper.__dict__));return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--configuration',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    cfg=json.loads(a.configuration.read_text());model=S.parents[2]/'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody_cvt3.xml'
    acts=[node.attrib['name']for node in ET.parse(model).getroot().find('actuator')]
    # Same actual Hunt controller/gain source, but no body bind or task-balance update.
    options=dict(weights_directory=cfg['weights_directory'],use_timing_afferents=cfg.get('timing_afferents',True),
                 axis_feedback=cfg.get('axis_feedback',0.),balance=cfg.get('balance'),timing_scale=cfg.get('timing_scale',1.),
                 reciprocal_gain=cfg.get('reciprocal_gain',0.),commissural_gain=cfg.get('commissural_gain',0.),
                 biarticular_pair_gain=cfg.get('biarticular_pair_gain',0.),bilateral_midline_drive=cfg.get('bilateral_midline_drive',False),
                 cross_contact_gain=cfg.get('cross_contact_gain',0.),learn_reciprocal=cfg.get('learn_reciprocal',False))
    n=C.build(acts,dt=.002,options=options);w=n.compiled;c=w.c
    unbound_balance=n.__dict__.pop('balance',None)
    try:n.apply_weights(1.);n.sync()
    finally:
        if unbound_balance is not None:n.balance=unbound_balance
    names=requested_cells();ids=np.array([n.idx[name]for name in names]);missing=[name for name in names if name not in n.idx]
    assert not missing
    before_first=sparse_current_snapshot(n,ids,'sparse');assert before_first['status']==5
    # Declared static fixture seeds one real chemical source; no claim of natural rhythm.
    dst=n.idx['MN_soleus_r'];edge=next(i for i,(post,pre)in enumerate(zip(w.post,w.pre))if post==dst and w.gmax[i]>0 and pre!=dst)
    seeded_source=int(w.pre[edge]);seed=float(w.elo[edge]+.6*w.erange[edge]);c.V[seeded_source]=seed
    u=n.make_inputs();u[n.ports['POST_soleus_r']]=.7;u[n.ports['POST_tib_ant_l']]=.3
    w.forward(u) # First declared fixture step; observation adds no step.
    before=fingerprint(n);observed=sparse_current_snapshot(n,ids,'sparse');after=fingerprint(n)
    assert before==after and observed['available'],observed
    direct_before=before;direct_after=after
    fields=observed['fields'];mn=[i for i,name in enumerate(names)if name.startswith('MN_')]
    closure=float(np.max(abs(fields['neural_update_residual_mV'][mn])))
    assert closure<=1e-12
    assert fields['neural_applied_current_nA'][names.index('MN_soleus_r')]==.7
    assert fields['neural_applied_current_nA'][names.index('MN_tib_ant_l')]==.3
    assert fields['neural_chemical_current_nA'][names.index('MN_soleus_r')]!=0.
    class Poison:
        def __getattribute__(self,name):raise AssertionError('Disabled/native observation touched a buffer')
    assert sparse_current_snapshot(Poison(),ids,None)['status']==0
    assert sparse_current_snapshot(Poison(),ids,'native')['status']==2
    stale=clone_wrapper(w);stale.forward=stale.c.forward;stale.sums[:]=np.nan
    native=sparse_current_snapshot(SimpleNamespace(compiled=stale),ids,'sparse');assert native['status']==3 and native['fields']is None
    general=clone_wrapper(w);general.input_onehot=False
    assert sparse_current_snapshot(SimpleNamespace(compiled=general),ids,'sparse')['status']==4
    bad=clone_wrapper(w);bad.c.V_last=bad.last_voltage;bad.i_app[dst]+=1.
    rejected=sparse_current_snapshot(SimpleNamespace(compiled=bad),ids,'sparse');assert rejected['status']==7
    # Exercise the real optional diagnostic wrapper, whose native closure must
    # point to this sparse instance. No live mechanics => only delegation runs.
    from swing_diagnostics import SwingDiagnostics
    diagnostic=SwingDiagnostics(n);w.forward(u)
    wrapped_before=fingerprint(n);wrapped=sparse_current_snapshot(n,ids,'sparse');wrapped_after=fingerprint(n)
    assert wrapped['available']and wrapped['forward_scope']=='audited_SwingDiagnostics_delegating_to_actual_sparse'
    assert wrapped_before==wrapped_after
    wrapped_native=clone_wrapper(w);wrapped_native.forward=wrapped_native.c.forward;wrapped_native.sums[:]=np.nan
    native_net=SimpleNamespace(compiled=wrapped_native,targets_data=n.targets_data,mn_edges=n.mn_edges,names=n.names,idx=n.idx,live=None,step_count=0)
    native_diagnostic=SwingDiagnostics(native_net)
    assert sparse_current_snapshot(native_net,ids,'sparse')['status']==3
    # Observation never advances the wrapped recorder or its native delegate.
    before_default=fingerprint(n)
    # Construction/default-off storage check against an existing exact compiled model;
    # no MjData, forward, kinematics, or physics operation is required.
    import mujoco
    campaign=a.configuration.parents[1]
    m=mujoco.MjModel.from_binary_path(str(campaign/'figures/source_baseline_s2/model.mjb'))
    default_capture=AlignedCapture(m,n,2)
    assert default_capture.current_buffers is None and default_capture._saved_current_arrays(0)=={}
    assert before_default==fingerprint(n)
    fields=wrapped['fields'];closure=float(np.max(abs(fields['neural_update_residual_mV'][mn])))
    a.out.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(a.out/'ACTUAL_SPARSE_CURRENT_FIXTURE.npz',neural_names=np.array(names),V_new=c.V[ids],**fields)
    result=dict(neurons=len(c.V),selected_cells=len(ids),input_onehot=bool(w.input_onehot),
                explicit_neural_fixture_steps=2,observation_additional_neural_steps=0,physics_steps=0,
                complete_state_and_buffer_fingerprint_before=wrapped_before,complete_state_and_buffer_fingerprint_after=wrapped_after,
                direct_observation_fingerprint_before=direct_before,direct_observation_fingerprint_after=direct_after,
                observation_nonmutating=direct_before==direct_after and wrapped_before==wrapped_after,maximum_MN_update_residual_mV=closure,
                maximum_all_selected_update_residual_mV=wrapped['maximum_closure_error_mV'],
                actual_direct_POST_currents_nA={'MN_soleus_r':.7,'MN_tib_ant_l':.3},
                seeded_static_chemical_source=dict(source=n.names[seeded_source],voltage_mV=seed,target='MN_soleus_r'),
                actual_soleus_chemical_current_nA=float(fields['neural_chemical_current_nA'][names.index('MN_soleus_r')]),
                default_off_does_not_read_buffers=True,default_off_adds_no_numeric_fields=True,
                native_config_missing=True,native_forward_with_stale_buffers_rejected=True,
                non_onehot_buffer_missing=True,inconsistent_current_closure_rejected=True,
                actual_SwingDiagnostics_sparse_wrapper_supported=True,diagnostic_wrapper_observation_nonmutating=wrapped_before==wrapped_after,
                diagnostic_wrapper_maximum_update_residual_mV=wrapped['maximum_closure_error_mV'],
                diagnostic_wrapping_native_forward_rejected=True,
                unsupported_fields='Missing, never synthesized zero; valid passive ionic zero is an actual zero.',
                activation_scope='This electrical fixture contains no mechanical activation data; live capture uses actual data.act separately from data.ctrl.',
                source_sha256={str(H/name):hashlib.sha256((H/name).read_bytes()).hexdigest()for name in('verification_capture.py','verify_sparse_current_capture.py')},
                scope='Actual compiled Hunt circuit; two explicit electrical fixture steps (direct and diagnostic-wrapped), no observation-induced step. No body performance, natural rhythm, or training claim.')
    (a.out/'STATIC_CURRENT_CLOSURE.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
