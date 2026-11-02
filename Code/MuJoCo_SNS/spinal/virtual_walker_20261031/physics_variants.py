"""Explicit opt-in runtime mechanical changes; source XML and geometry stay intact."""
import math
import numpy as np
import mujoco

def apply_segment_mass_scale(model,data,scale=1.0):
 """Scale bilateral limb density at fixed shape, including its principal inertia.

 Call after model compilation and before first forward/measurement. The default
 returns without mj_setConst or any array assignment, preserving the baseline.
 """
 scale=float(scale)
 if not math.isfinite(scale) or scale<=0:raise ValueError('segment_mass_scale must be finite and positive')
 names=[segment+'_'+side for side in ('L','R') for segment in ('femur','tibia','foot','toe')]
 ids=[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,name) for name in names]
 if any(i<0 for i in ids):raise ValueError('Mass variation requires all eight named bilateral limb bodies')
 ledger=[dict(body=name,body_id=int(i),mass_before_kg=float(model.body_mass[i]),inertia_before_kg_m2=model.body_inertia[i].copy().tolist(),mass_after_kg=float(model.body_mass[i]*scale),inertia_after_kg_m2=(model.body_inertia[i]*scale).tolist()) for name,i in zip(names,ids)]
 before=float(np.sum(model.body_mass[ids]))
 if scale!=1.0:
  model.body_mass[ids]*=scale;model.body_inertia[ids]*=scale
  mujoco.mj_setConst(model,data)
 return dict(segment_mass_scale=scale,applied=scale!=1.,body_ledger=ledger,selected_mass_before_kg=before,selected_mass_after_kg=before*scale,
  geometry_changed=False,source_modified=False,scope='Fixed-geometry bilateral limb-density perturbation only; pelvis, apparatus, contact-only bodies, neural parameters, and tendon paths unchanged')

def pool_inherited_hip_law(controller,bank,force_params,model,data):
 """Explicit arithmetic pooling by verified MN function, not actuator labels.

 Both source hips contribute equally. This is a symmetry counterfactual, not
 an identified hardware law. All other actuator fields are already paired.
 """
 outputs={controller.names[pre]:(pre,index) for pre,index,coeff in controller.outputs};ledger=[];indices={}
 for action,expected_sign in [('flx',1),('ext',-1)]:
  pair=[outputs[side+' Hip MN '+action][1] for side in 'LR'];before=force_params[pair].copy();other=[i for i in range(force_params.shape[1]) if i not in (7,8)];assert np.array_equal(before[0,other],before[1,other]),'Unexpected additional hip-law asymmetry needs its own audited policy'
  pooled=np.mean(before[:,[7,8]],axis=0);force_params[np.ix_(pair,[7,8])]=pooled[None,:];rows=[]
  for side,index,old in zip('LR',pair,before):
   joint=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,'hip_'+side);tid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_TENDON,'route_'+bank[index]['name']);arm=-float(data.ten_J.reshape(model.ntendon,model.nv)[tid,model.jnt_dofadr[joint]]);assert arm*expected_sign>0,'Hip MN function does not match actual route action'
   rows.append(dict(side=side,motor_neuron=side+' Hip MN '+action,physical_route=bank[index]['name'],physical_route_id=bank[index]['id'],home_route_length_m=float(data.ten_length[tid]),home_tension_moment_arm_m=arm,sigmoid_B_before_N=float(old[7]),sigmoid_C_before_per_V=float(old[8]),sigmoid_B_after_N=float(pooled[0]),sigmoid_C_after_per_V=float(pooled[1])));indices[(side,action)]=index
  ledger.append(dict(function='hip_flexor' if action=='flx' else 'hip_extensor',rule='Equal arithmetic mean of the two source parameter fields',routes=rows))
 input_indices=[]
 for side in 'LR':
  matches=[i for i,row in enumerate(controller.inputs) if controller.names[row[0]]==side+' Hip ext Ib' and row[2]=='Tension'];assert len(matches)==1;i=matches[0];assert controller.inputs[i][1]==indices[(side,'ext')];input_indices.append(i)
 original=[controller.inputs[i][3].copy() for i in input_indices];pooled=np.mean(original,axis=0);sensor=[]
 for side,i,before in zip('LR',input_indices,original):
  target,index,kind,coeff,body=controller.inputs[i];controller.inputs[i]=(target,index,kind,pooled.copy(),body);sensor.append(dict(side=side,target=controller.names[target],physical_route=bank[index]['name'],source_data_type=kind,polynomial_before_A_per_input=before.tolist(),polynomial_after_A_per_input=pooled.tolist()))
 return dict(enabled=True,actuator_ledger=ledger,legacy_hip_extensor_Ib_ledger=sensor,source_modified=False,neural_weights_changed=False,geometry_or_mass_changed=False,scope='Only functionally paired hip sigmoid B/C and legacy hip-extensor force-to-Ib polynomial pooled. Other actuator fields, geometry, masses, passive springs, neural gains, and new geometry-normalized II converters retained. Not a claim that the entire source mechanics is perfectly symmetric.')
