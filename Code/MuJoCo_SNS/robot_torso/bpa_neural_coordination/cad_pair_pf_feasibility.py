"""Read-only source attribution and conditional descending-gain feasibility."""
import json
import numpy as np
from scipy.optimize import lsq_linear
import cad_pair_pf_audit as A
OUT=A.OUT
def main():
 with np.load(OUT/'TEACHER_NEURAL.npz') as z:d={k:z[k] for k in z.files}
 names=d['neuron_names'].tolist();ids=[names.index('MN_vas_med_'+s) for s in ('r','l')];vm=d['V_after'][:,ids];old=d['V_before'];gmax=d['static_g_max_non'][ids];lo=d['static_e_lo'][ids];hi=d['static_e_hi'][ids];rev=d['static_del_e'][ids];g=np.clip(gmax[None]*(old[:,None,:]-lo[None])/np.where(hi==lo,1,hi-lo)[None],0,gmax[None]);curr=g*(rev[None]-vm[:,:,None]);rows=[]
 for side in range(2):
  incoming=[]
  for index in np.flatnonzero(gmax[side]):incoming.append(dict(source=names[index],gmax_uS=float(gmax[side,index]),voltage_range_mV=[float(old[:,index].min()),float(old[:,index].max())],mean_current_nA=float(curr[:,side,index].mean()),min_current_nA=float(curr[:,side,index].min()),max_current_nA=float(curr[:,side,index].max())))
  rows.append(dict(target=names[ids[side]],incoming=sorted(incoming,key=lambda r:-abs(r['mean_current_nA']))))
 pf=np.array([n.startswith('PF_') for n in names]);posture=np.array([n=='POSTURE' for n in names]);descending=pf|posture;ports=d['input_names'].tolist();app=d['inputs']@d['static_input_connectivity'][ids].T;direct=np.column_stack([d['inputs'][:,ports.index('POST_vas_med_'+s)] for s in ('r','l')]);non=curr[:,:,~descending].sum(axis=2)+d['static_i_b'][ids]-d['static_g_m'][ids]*(vm-d['static_V_rest'][ids])+app-direct;drive=curr[:,:,descending].sum(axis=2);gain=float(np.clip(-np.sum(drive*non)/np.sum(drive**2),0,1));residual=drive*gain+non
 # Keep all phase/tonic ratios with one shared normalization rather than fit
 # indistinguishable saturated phase columns independently.
 c=np.column_stack((curr[:,:,pf].sum(axis=2).reshape(-1),curr[:,:,posture].sum(axis=2).reshape(-1)));s=np.linalg.svd(c,compute_uv=False);two=lsq_linear(c,-non.reshape(-1),bounds=(0,1));r=dict(incoming_by_side=rows,PF_and_POSTURE_design_singular_values=s.tolist(),two_weight_nonnegative_fit=two.x.tolist(),shared_all_descending_gain_diagnostic=gain,shared_gain_RMS_current_residual_nA=float(np.sqrt(np.mean(residual**2))),shared_gain_max_abs_current_residual_nA=float(np.max(np.abs(residual))),scope='Conditional read-only diagnosis, not an authorized candidate or physicaltrial. No circuit/sourceweights changed. POSTURE is existingneuralpresynapse, not intrinsicbias/directMNclamp. Shared gain would preservePF/POSTURE ratiosandbilateral symmetry.',proposed_changed_edges=[dict(source=names[j],targets=[names[ids[i]] for i in range(2) if gmax[i,j]>0]) for j in np.flatnonzero(descending & np.any(gmax>0,axis=0))]);(OUT/'TONIC_AND_PF_FEASIBILITY.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
if __name__=='__main__':main()
