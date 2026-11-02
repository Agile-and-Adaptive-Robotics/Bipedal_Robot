"""Export the actual CAD core, verified against saved fresh numerical arrays."""
import json
import numpy as np
import cad_mixed_coactivation05 as J
R=J.C.R;OUT=J.OUT
def main():
 m=J.C.mujoco.MjModel.from_xml_path(str(J.C.ROBOT/'ah001_neural.xml'));names=[m.actuator(i).name for i in range(m.nu)];R.params.G.update(c1_gain=.15,v3_gain=0.,heel_rge=.15,toe_rge=.05);R.params.TAU['rg_nap_h']=.35;n=R.bn.build(names,dt=m.opt.timestep,interleg=True);c=n.compiled
 with np.load(OUT/'fresh_a/mixed/FINAL_COMPLETE_STATE.npz') as z:
  checks={k:bool(np.array_equal(getattr(c,k),z['backend_'+k])) for k in ('g_max_non','del_e','e_lo','e_hi','input_connectivity','output_voltage_connectivity')}
 assert all(checks.values()),'Rebuilt connectome differs from actual saved run'
 ordered=sorted(n.idx,key=n.idx.get);nodes=[dict(id=i,name=name,group='MN' if name.startswith('MN_') else name.split('_')[0]) for i,name in enumerate(ordered)];edges=[]
 for post,pre in zip(*np.nonzero(c.g_max_non)):
  edges.append(dict(source=int(pre),target=int(post),gmax_uS=float(c.g_max_non[post,pre]),reversal_mV=float(c.del_e[post,pre]),e_lo_mV=float(c.e_lo[post,pre]),e_hi_mV=float(c.e_hi[post,pre])))
 output=dict(machine='easteregg2',scope='Actual262cell54actuator CAD SNS core, numericalmatrixverified against fixed5% pressurepair run; engineering allocator/MNcurrentclamp external to this graph must be shown.',nodes=nodes,edges=edges,inputs=[dict(port=i,name=name,targets=[dict(neuron=int(j),weight=float(c.input_connectivity[j,i])) for j in np.flatnonzero(c.input_connectivity[:,i])]) for i,name in enumerate(n.inputs)],verification=checks,external_blocks=[dict(name='Posture torque allocation',role='Produces descending MN excitation targets'),dict(name='MN current clamp',role='Cancels instantaneous background MN synaptic current and imposes membrane target; retain in schematics'),dict(name='Two pressure actuators',routes=['vas_med_r','vas_med_l'],conversion='SNSMN/5mV→normalizedexcitation→MIF4530N reference→nominalstaticinverse(L)→4msdelay→40/60ms pressure→BPAseriesforce'),dict(name='Other52 reference tensions',role='Normalizedexcitation→30msactivation→masterMIF-scaled reference force'),dict(name='Full moving-route forces',role='−J_length.T@positive_tensions, applied once; nativegain/bias0'),dict(name='Pelvis apparatus',role='Independent measured translationalforces and rotationaltorques; ~91.7%BW upwardassistance')])
 (OUT/'ACTUAL_CORE_CONNECTOME.json').write_text(json.dumps(output,indent=2));print(json.dumps(dict(nodes=len(nodes),edges=len(edges),inputs=len(n.inputs),verified=True),indent=2))
if __name__=='__main__':main()
