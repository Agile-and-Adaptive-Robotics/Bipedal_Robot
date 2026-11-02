"""Read-only additive RG/Ib attribution and conditional recruitment replay."""
import json
import numpy as np
import cad_pair_pf_audit as A
OUT=A.OUT;C=A.C
def main():
 with np.load(OUT/'TEACHER_NEURAL.npz') as z:d={k:z[k] for k in z.files}
 names=d['neuron_names'].tolist();ports=d['input_names'].tolist();v=d['V_before'];gmax=d['static_g_max_non'];lo=d['static_e_lo'];hi=d['static_e_hi'];rev=d['static_del_e'];span=np.where(hi==lo,1,hi-lo);upstream=[]
 for side in ('r','l'):
  i=names.index('IBEXC_knee_ext_'+side);g=np.clip(gmax[i]*(v-lo[i])/span[i],0,gmax[i]);I=g*(rev[i]-v[:,i,None]);rows=[]
  for j in np.flatnonzero(gmax[i]):rows.append(dict(source=names[j],gmax_uS=float(gmax[i,j]),mean_current_nA=float(I[:,j].mean()),min_current_nA=float(I[:,j].min()),max_current_nA=float(I[:,j].max())))
  rg=np.array([n.startswith('RG_') for n in names]);ib=np.array([n.startswith('Ib_') for n in names]);rgI=I[:,rg].sum(axis=1);ibI=I[:,ib].sum(axis=1);upstream.append(dict(node=names[i],incoming=rows,mean_RG_current_nA=float(rgI.mean()),mean_Ib_current_nA=float(ibI.mean()),Ib_fraction_positive_drive=float(np.maximum(ibI,0).sum()/np.maximum(I,0).sum()),zero_Ib_instantaneous_RG_only_equilibrium_mV=(g[:,rg]*rev[i,rg]).sum(axis=1).mean()/(d['static_g_m'][i]+g[:,rg].sum(axis=1).mean())))
 ids=np.array([names.index('MN_vas_med_'+s) for s in ('r','l')]);target=d['V_after'][:,ids];g=np.clip(gmax[ids][None]*(v[:,None]-lo[ids][None])/span[ids][None],0,gmax[ids][None]);mask=np.array([n.startswith(('PF_','IBEXC_')) for n in names]);drive=(g[:,:,mask]*(rev[ids][:,mask][None]-target[:,:,None])).sum(axis=2);app=d['inputs']@d['static_input_connectivity'][ids].T;app-=np.column_stack([d['inputs'][:,ports.index('POST_vas_med_'+s)] for s in ('r','l')]);other=(g[:,:,~mask]*(rev[ids][:,~mask][None]-target[:,:,None])).sum(axis=2)+d['static_i_b'][ids]-d['static_g_m'][ids]*(target-d['static_V_rest'][ids])+app;gain=float(np.clip(-np.sum(drive*other)/np.sum(drive**2),0,1));scaled=g.copy();scaled[:,:,mask]*=gain
 # Conditional isolated-neuron replay, not a free-network or live model test.
 vm=d['static_V_rest'][ids].copy();trace=[]
 def advance(index):
  nonlocal vm
  I=-d['static_g_m'][ids]*(vm-d['static_V_rest'][ids])+d['static_i_b'][ids]+(scaled[index]*(rev[ids]-vm[:,None])).sum(axis=1)+app[index];vm+=d['static_time_factor_membrane'][ids]*I
 for _ in range(100):advance(0)
 for k in range(len(v)):advance(k);trace.append(vm.copy())
 trace=np.array(trace);cfg=json.loads((A.OLD/'ASSUMED_PAIR_CONFIGURATION.json').read_text());bank=C.RoutedBPABank(C.PAIR,cfg,.002);m=C.mujoco.MjModel.from_xml_path(str(C.ROBOT/'ah001_neural.xml'));aids=[m.actuator(n).id for n in C.PAIR];offset=np.array([cfg['actuators'][n]['constant_remote_path_offset_m'] for n in C.PAIR]);physical=[];error=None
 with np.load(OUT/'teacher/mixed/ACTUATOR_BRIDGE.npz') as z:L=z['length'][:,aids];actual=z['force'][:,aids]
 for k in range(len(trace)):
  try:bank.step(np.clip(trace[k]/5,0,1),L[k]+offset)
  except ValueError as e:error=str(e);break
  physical.append(np.r_[bank.command,bank.pressure,bank.force,bank.flags])
 physical=np.asarray(physical);report=dict(IBEXC_upstream=upstream,source_rule='build_network.py939–943: RG_E→IBEXC excitatoryg1, Ib→IBEXC excitatory, IBEXC→MN excitatoryg0.6. Additiveinputs; no multiplicativeloadgate.',conditional_PF_and_IBEXC_common_gain=gain,POSTURE_preserved=True,conditional_current_RMS_residual_nA=float(np.sqrt(np.mean((drive*gain+other)**2))),isolated_MN_voltage_range_mV=np.column_stack((trace.min(axis=0),trace.max(axis=0))).tolist(),conditional_pressure_failure=error,conditional_pressure_peak_kpa=float(physical[:,2:4].max()) if len(physical) else None,conditional_force_RMS_error_N=np.sqrt(np.mean((physical[:,4:6]-actual[:len(physical)])**2,axis=0)).tolist() if len(physical) else None,conditional_inverse_flags=np.unique(physical[:,6:8]).tolist() if len(physical) else [],scope='Read-only conditional fit extendsbeyond PF-only scope by scalingexisting stanceoutputweights aswell. POSTURE/leak/tonic/sensory unchanged. Teacherpresynapticvoltages+geometry frozen; no sourceweights or free-network/live controller changed. Requires explicit parentapprovalbeforeuse.')
 np.savez_compressed(OUT/'CONDITIONAL_CENTRAL_REPLAY.npz',t=d['time'],MN_mV=trace,pressure_force=physical);(OUT/'STANCE_PATH_AUDIT.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__=='__main__':main()
