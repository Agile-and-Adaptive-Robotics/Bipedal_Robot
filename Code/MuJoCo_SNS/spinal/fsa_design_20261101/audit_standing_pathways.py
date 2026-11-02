"""Static inherited BAL and posture audit; no integration or training."""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import argparse,json,sys
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np,mujoco
sys.path.insert(0,str(Path(__file__).resolve().parent))
import integrated_hunt_design as I


def run(root,out):
    cfgpath=root/'configurations/candidate_ii02_s3.json';cfg=json.loads(cfgpath.read_text());saved=json.loads(cfgpath.with_suffix('.result.json').read_text());capture=Path(cfg['verification_capture']);m=mujoco.MjModel.from_binary_path(str(capture/'model.mjb'));d=mujoco.MjData(m);raw=np.load(cfgpath.with_suffix('.npz'));cap=np.load(capture/'aligned_trace.npz');mech=np.load(cfgpath.with_suffix('.mechanics.npz'),allow_pickle=True)
    d.qpos[:]=raw['qfull'][0];mujoco.mj_fwdPosition(m,d)
    # Archive was saved after the declared harness lowering. Reconstruct the
    # ORIGINAL bind reference, not a new rig or a different posture.
    displacement=np.array([0.,0.,cfg.get('harness_lowering_m',0.)+cfg.get('harness_vertical_offset_m',0.)])
    for axis in ('tx','ty','tz'):
        j=m.joint('pelvis_'+axis);m.qpos_spring[j.qposadr[0]]-=float(d.xaxis[j.id]@displacement)
    acts=[q.attrib['name']for q in ET.parse(I.MODEL).getroot().find('actuator')]
    n=I.build(acts,options=I.condition_options(cfg),profile=cfg['integrated_fsa_profile'],protocol_stage=3);n.bind_mujoco(m,d);n.apply_weights(0.);I.verify(n,full=True,require_bound_PT=True)
    fit=np.asarray(n.balance.calibration['local_acceleration_response']);expected=np.asarray(saved['balance_calibration']['local_acceleration_response']);cal_error=float(np.max(abs(fit-expected)))
    if cal_error>1e-8:raise AssertionError(('Original bound response not reproduced',cal_error,fit.tolist()))
    R=n.balance.response.copy();c=n.compiled;names=mech['names'].tolist();order=mech['acts'].tolist();V=mech['voltage'];mn=np.array([n.idx['MN_'+a]for a in order]);edges=n.mn_edges;dest=np.array([order.index(n.names[c.post[e]][3:])for e in edges])
    # Capture is post-update/control; mechanics gives exact pre-update voltage.
    cap_index=np.array([int(round((t-.6)/.002))for t in mech['t']]);valid=(cap_index>=0)&(cap_index<len(cap['t_state']));ix=np.flatnonzero(valid);ci=cap_index[valid]
    if not np.allclose(cap['t_neural_input'][ci],mech['t'][ix],rtol=0,atol=1e-9):raise AssertionError('Neural timestamps do not pair')
    ctrl=cap['control'][ci];postV=5.*ctrl;known=(ctrl>1e-9)&(ctrl<1.-1e-9)
    for j,a in enumerate(order):
        label='MN_'+a
        if label in cap['neural_names'].tolist():postV[:,j]=cap['neural_voltage'][ci,cap['neural_names'].tolist().index(label)];known[:,j]=True
    app=(c.c.c_m[mn]*(postV-V[ix][:,mn])/.002-mech['leak_and_bias'][ix]-mech['currents'][ix].sum(axis=1))
    direct=np.full(len(order),np.nan);spread=np.full(len(order),np.nan)
    for j in range(len(order)):
        vals=app[known[:,j],j]
        if len(vals):direct[j]=np.median(vals);spread[j]=np.ptp(vals)
    if np.nanmax(spread)>1e-9:raise AssertionError(('Direct input is not constant/reconstruction failed',float(np.nanmax(spread))))
    source_record={};isolated={};probe=mujoco.MjData(m);probe.qpos[:]=d.qpos;probe.act[:]=.05;mujoco.mj_forward(m,probe);jac=np.zeros((3,m.nv));mujoco.mj_jacSubtreeCom(m,probe,jac,0);native_base=jac@probe.qacc
    for source in ('BAL_PF','BAL_DF','POSTURE'):
        idx=n.idx[source];ids=np.flatnonzero(c.pre[edges]==idx);gain=c.gmax[edges[ids]];vm=V[:,idx];alpha=np.clip((vm[:,None]-c.elo[edges[ids]])/c.erange[edges[ids]],0,1)
        current=gain*alpha*(c.reversal[edges[ids]]-V[:,c.post[edges[ids]]]);delta=np.zeros(len(order))
        isolated_alpha=np.clip((1./c.c.g_m[idx]-c.elo[edges[ids]])/c.erange[edges[ids]],0,1);geff=gain*isolated_alpha
        delta[dest[ids]]=(geff*(c.reversal[edges[ids]]-2.5)/(c.c.g_m[c.post[edges[ids]]]+geff))/5.
        isolated[source]={'input_nA':1.,'MN_baseline_mV':2.5,'actual_eligible_MN_edges':len(ids),'gain_min_max_uS':[float(gain.min()),float(gain.max())],'predicted_MN_activation_increment_L1':float(delta.sum()),'ground_reference_acceleration_delta_native_xyz_m_s2':(R@delta).tolist(),'scope':'Frozen initial ground-bound activation response; isolated steady source/MN operator including shunting,not a body rollout.'}
        matched=np.zeros(len(order));matched[dest[ids]]=(geff*(c.reversal[edges[ids]]-.25)/(c.c.g_m[c.post[edges[ids]]]+geff))/5.
        probe.act[:]=.05+matched;mujoco.mj_forward(m,probe)
        isolated[source]['matched_MN_baseline_mV']=.25
        isolated[source]['matched_mechanical_activation_baseline']=.05
        isolated[source]['matched_reference_linear_acceleration_delta_m_s2']=(R@matched).tolist()
        isolated[source]['exact_native_static_forward_acceleration_delta_m_s2']=(jac@probe.qacc-native_base).tolist()
        isolated[source]['static_forward_scope']='Actual native nonlinear contact/constraint forward dynamics at the original fixed pose,zero velocity,activation0.05 plus the isolated electrical increment. No integration or neural closed-loop rollout.'
        keep=(mech['t']>=5.)&(mech['t']<8.)
        source_record[source]={'late_source_voltage_mean_mV':float(vm[keep].mean()),'late_source_voltage_min_max_mV':[float(vm[keep].min()),float(vm[keep].max())],'late_total_MN_current_mean_nA':float(current[keep].sum(axis=1).mean()),'MN_edges':[{'actuator':order[dest[i]],'g_uS':float(gain[j]),'late_current_nA':float(current[keep,j].mean())}for j,i in enumerate(ids)]}
    # All MN incoming conductances, including recorded pose-dependent recruitment.
    gs=np.broadcast_to(c.gmax[edges],(len(V),len(edges))).copy();locations={int(e):i for i,e in enumerate(edges)}
    for e,act,row,sign in zip(n.geom_edges,n.geom_acts,n.geom_rows,n.geom_signs):
        gs[:,locations[int(e)]]=n.ff_gain[e]*np.clip(sign*mech['rho'][:,row,act]/n.geom_scale[row],0,1)
    alpha=np.clip((V[:,c.pre[edges]]-c.elo[edges])/c.erange[edges],0,1);geff=gs*alpha
    syn=np.zeros((len(V),len(order)));shunt=np.zeros_like(syn);reversal_sum=np.zeros_like(syn)
    for j,e in enumerate(edges):
        a=dest[j];syn[:,a]+=geff[:,j]*(c.reversal[e]-V[:,c.post[e]]);shunt[:,a]+=geff[:,j];reversal_sum[:,a]+=geff[:,j]*c.reversal[e]
    syn_error=float(np.max(abs(syn-mech['currents'].sum(axis=1))))
    if syn_error>1e-9:raise AssertionError(('Full recorded synaptic-current replay mismatch',syn_error))
    selected={};late=(mech['t']>=5)&(mech['t']<8)
    posture_idx=n.idx['POSTURE']
    for act in ('iliacus_r','bifemlh_r','vas_int_r','soleus_r','tib_ant_r'):
        j=order.index(act);a_static=direct[j]/6.;pids=np.flatnonzero((c.post==n.idx['MN_'+act])&(c.pre==posture_idx));pcurrent=sum(c.gmax[e]*np.clip((V[:,posture_idx]-c.elo[e])/c.erange[e],0,1)*(c.reversal[e]-V[:,n.idx['MN_'+act]])for e in pids)
        selected[act]={'inferred_actual_POST_direct_current_nA':float(direct[j]),'constant_input_reconstruction_spread_nA':float(spread[j]),'implied_static_activation_target':float(a_static),'direct_only_steady_activation_without_any_synapses':float(np.clip(direct[j]/(c.c.g_m[n.idx['MN_'+act]]*5.),0,1)),'global_POSTURE_chemical_current_late_mean_nA':float(np.mean(pcurrent[late])),'actual_activation_late_mean':float(np.mean(mech['act'][late,j])),'actual_control_late_mean':float(np.mean(mech['ctrl'][late,j]))}
    result={'schema':'inherited_balance_and_posture_static_audit_v1','configuration':str(cfgpath),'source_hashes':{'configuration':I.sha(cfgpath),'audit':I.sha(__file__),'model':I.sha(capture/'model.mjb')},'original_ground_PT_calibration_reproduction_max_error_m_s2':cal_error,'initial_ground_activation_response_matrix':R.tolist(),'actuator_order':order,'full_recorded_synaptic_current_replay_max_error_nA':syn_error,'actual_inherited_source_records':source_record,'isolated_current_to_MN_to_ground_acceleration':isolated,'posture_inputs':{'actual_direct_POST_current_decoded_count':int(np.isfinite(direct).sum()),'largest_time_variation_in_decoded_direct_POST_nA':float(np.nanmax(spread)),'selected_MNs':selected,'all_direct_POST_currents_nA':direct.tolist(),'meaning':'POST_<act> is direct applied current6*staticactivation at standing sf1; global POSTURE chemical input is additionally active. MN leak and all synaptic shunts affect the delivered activation. Even direct-only current6*a gives1.2*a at1uS leak and5mV readout,not a. This audit reconstructs actual POST current from pre/post neural states and actualctrl; it does not confuse legacy act=ctrl with physicalactivation.'},'physics_steps':0,'neural_steps':0,'training':False,'scope':'Native forward-dynamics evaluations recreate the saved fixed-pose activation response,then algebra evaluates actual compiled neural operators. No integration or production parameter change. Initial bound response is local; extrapolation to moving forefoot-only states is not a closed-loop proof.'}
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps({'output':str(out),'R_reproduction_error':cal_error,'synaptic_replay_error':syn_error,'isolated_BAL_accel':{k:v['ground_reference_acceleration_delta_native_xyz_m_s2']for k,v in isolated.items()}}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();run(a.root,a.out)
