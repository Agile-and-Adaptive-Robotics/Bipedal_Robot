"""Measured CAD teacher currents, pressure-map replay, and shared PF-only fit."""
import json,hashlib,shutil,sys
from pathlib import Path
import numpy as np
from scipy.optimize import lsq_linear
import cad_mixed_pressure_followup as K
J=K.J;C=K.C;L=K.L;OLD=J.OUT;OUT=OLD.parent/'pair_neural_recruitment';HERE=Path(__file__).resolve().parent

def save(name,d):OUT.mkdir(exist_ok=True);(OUT/name).write_text(json.dumps(d,indent=2),encoding='utf-8')
def teacher():
 OUT.mkdir(exist_ok=True);(OUT/'teacher/reference').mkdir(parents=True,exist_ok=True)
 for name in ('ASSUMED_PAIR_CONFIGURATION.json','STATIC_AUDIT.json'):shutil.copyfile(OLD/name,OUT/name)
 ref=json.loads((OLD/'reference/BRIDGE_SUMMARY.json').read_text());ref.update(reference_fixture_only_not_new_run=True,original_reference_summary=str(OLD/'reference/BRIDGE_SUMMARY.json'));(OUT/'teacher/reference/BRIDGE_SUMMARY.json').write_text(json.dumps(ref,indent=2))
 K.BASE=OUT;original=C.R.bn.build;records=[];holder=[]
 def observe_build(*args,**kwargs):
  n=original(*args,**kwargs);holder.append(n);original_step=n.step
  def step(x):
   before=n.compiled.V.copy();out=original_step(x);a=L.ACTIVE;records.append((float(a.d.time),before,n.compiled.V.copy(),x.copy()));return out
  n.step=step;return n
 C.R.bn.build=observe_build
 registration=dict(machine='easteregg2',phase='Read-only teacher audit before any experimental fit/live release',controller='Originalhybrid all54 MNclampedteacher, twoactualpressureactuators; NOT experimentalneuralpair',teacher_duration_s=2.,fit_scope='Shared bilateral multiplier per existingPFphase→vas_medMNedge; no allocator target imposed at experimentalpairMN; other52remainexternalbaseline',fit_limits=[0.,20.],target_source='ActualSNS MN outputs from meaningfulground-supported pressurepair record, not measuredEMG;5%normalizedcommand is not an actual5%forcefloor',pressure_map_fixed=True,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),Path(J.__file__),Path(C.__file__),Path(L.__file__))});save('AUDIT_PROTOCOL.json',registration)
 K.child('teacher');n=holder[0];ordered=sorted(n.idx,key=n.idx.get);t=np.array([r[0] for r in records]);last=np.r_[np.flatnonzero(np.diff(t)>0),len(t)-1];np.savez_compressed(OUT/'TEACHER_NEURAL.npz',time=t[last],V_before=np.array([r[1] for r in records])[last],V_after=np.array([r[2] for r in records])[last],inputs=np.array([r[3] for r in records])[last],neuron_names=np.array(ordered),input_names=np.array(n.inputs),**{'static_'+k:getattr(n.compiled.c,k).copy() for k in ('g_max_non','e_lo','e_hi','del_e','g_m','V_rest','i_b','input_connectivity','time_factor_membrane')})

def audit():
 if not (OUT/'TEACHER_NEURAL.npz').exists():teacher()
 with np.load(OUT/'TEACHER_NEURAL.npz') as z:d={k:z[k] for k in z.files}
 names=d['neuron_names'].tolist();ports=d['input_names'].tolist();ids=np.array([names.index('MN_vas_med_'+s) for s in ('r','l')]);mnv=d['V_before'][:,ids];target=d['V_after'][:,ids];gmax=d['static_g_max_non'][ids];lo=d['static_e_lo'][ids];hi=d['static_e_hi'][ids];rev=d['static_del_e'][ids];span=np.where(hi==lo,1.,hi-lo);g=np.clip(gmax[None,:,:]*(d['V_before'][:,None,:]-lo[None,:,:])/span[None,:,:],0,gmax[None,:,:]);curr=g*(rev[None,:,:]-mnv[:,:,None]);groups={}
 for prefix in ('PF_','Ia_','II_','Ib_','RC_','IaIN_','IIIN_'):
  selected=np.array([n.startswith(prefix) for n in names]);groups[prefix]=curr[:,:,selected].sum(axis=2)
 groups['other_synapses']=curr.sum(axis=2)-sum(groups.values());groups['tonic']=np.broadcast_to(d['static_i_b'][ids],mnv.shape).copy();groups['leak']=-d['static_g_m'][ids]*(mnv-d['static_V_rest'][ids]);allapp=d['inputs']@d['static_input_connectivity'][ids].T;postports=np.array([ports.index('POST_vas_med_'+s) for s in ('r','l')]);post=np.column_stack([d['inputs'][:,p] for p in postports]);groups['direct_POST_clamp']=post;groups['other_applied_input']=allapp-post
 stats={k:dict(mean_R_L=v.mean(axis=0).tolist(),min_R_L=v.min(axis=0).tolist(),max_R_L=v.max(axis=0).tolist()) for k,v in groups.items()};phases=('E1','E2','F1','F2');A=[]
 # Linear conductance equilibrium at the desired measured membrane waveform.
 for phase in phases:
  mask=np.array([n.startswith('PF_'+phase+'_') for n in names]);A.append((g[:,:,mask]*(rev[None,:,mask]-target[:,:,None])).sum(axis=2))
 A=np.stack(A,axis=2);nonpf=~np.array([n.startswith('PF_') for n in names]);other=(g[:,:,nonpf]*(rev[None,:,nonpf]-target[:,:,None])).sum(axis=2)+groups['tonic']-d['static_g_m'][ids]*(target-d['static_V_rest'][ids])+groups['other_applied_input'];b=-other
 active=np.any(np.abs(A.reshape(-1,4))>1e-12,axis=0);fit=lsq_linear(A.reshape(-1,4)[:,active],b.reshape(-1),bounds=(0.,20.),tol=1e-12);gain=np.ones(4);gain[active]=fit.x;residual=(A@gain+other);steady0=groups['tonic']+groups['other_applied_input']+(g*rev[None,:,:]).sum(axis=2);den=d['static_g_m'][ids]+g.sum(axis=2);unclamped_equilibrium=steady0/den
 # Verify the fixed physical pressure interface using the real ground record.
 with np.load(OUT/'teacher/mixed/ACTUATOR_BRIDGE.npz') as z:p={k:z[k] for k in z.files}
 cfg=json.loads((OLD/'ASSUMED_PAIR_CONFIGURATION.json').read_text());bank=C.RoutedBPABank(C.PAIR,cfg,.002);m=C.mujoco.MjModel.from_xml_path(str(C.ROBOT/'ah001_neural.xml'));aids=[m.actuator(n).id for n in C.PAIR];offset=np.array([cfg['actuators'][n]['constant_remote_path_offset_m'] for n in C.PAIR]);physical=[]
 for u,length in zip(p['neural_output'][:,aids],p['length'][:,aids]):bank.step(u,length+offset);physical.append(np.r_[bank.pressure,bank.force])
 physical=np.array(physical);pressure_exact=bool(np.array_equal(physical[:,:2],p['pressure']));force_exact=bool(np.array_equal(physical[:,2:],p['force'][:,aids]));assert pressure_exact and force_exact,'Actual neural→pressure→force replay mismatch'
 report=dict(machine='easteregg2',target_MN_voltage_range_mV=np.column_stack((target.min(axis=0),target.max(axis=0))).tolist(),target_normalized_command_range=np.column_stack((target.min(axis=0)/5,target.max(axis=0)/5)).tolist(),original_unclamped_instantaneous_equilibrium_range_mV=np.column_stack((unclamped_equilibrium.min(axis=0),unclamped_equilibrium.max(axis=0))).tolist(),current_components_nA=stats,PF_phases=list(phases),sensitive_phase_columns=active.tolist(),shared_phase_multiplier_fit=gain.tolist(),fit_RMS_current_residual_nA=float(np.sqrt(np.mean(residual**2))),fit_max_abs_current_residual_nA=float(np.max(np.abs(residual))),required_other_current_without_PF_range_nA=np.column_stack((b.min(axis=0),b.max(axis=0))).tolist(),pressure_replay_arrays_exact=pressure_exact,actual_tension_replay_arrays_exact=force_exact,fit_scope='Offlinefixedpresynapticrecord conductance-equilibrium fit only; no claim fullfreeMNnetwork matches waveform. Existing presynapticPFonly, multiplierssharedbothlegs, bounds0..20. Originaltonic/leak/inhibition remain and directPOSTclamp excluded.',physical_trial_authorized_by_this_audit=False)
 np.savez_compressed(OUT/'CURRENT_AUDIT.npz',time=d['time'],target_MN_mV=target,free_equilibrium_mV=unclamped_equilibrium,PF_current_design=A,other_current_at_target=other,fit_current_residual=residual,**groups);save('CURRENT_AUDIT.json',report);print(json.dumps(report,indent=2));return report
if __name__=='__main__':audit()
