"""Independent actual OpenSim Thelen static force/fiber equilibrium check."""
import os,sys,json
from pathlib import Path
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np
import opensim as osim
from series_elastic import SeriesElasticBank
OUT=Path(__file__).resolve().parent/'runs/audit_20261030';OUT.mkdir(parents=True,exist_ok=True)
rows=[]
for activation in (.05,.3,.8):
 for angle in (0.,.3):
  model=osim.Model();model.setGravity(osim.Vec3(0))
  muscle=osim.Thelen2003Muscle('test',1000.,.1,.2,angle)
  muscle.set_FmaxTendonStrain(.033)
  muscle.set_KshapeActive(.5);muscle.set_KshapePassive(4.);muscle.set_FmaxMuscleStrain(.6)
  muscle.addNewPathPoint('origin',model.getGround(),osim.Vec3(0))
  L=.2+.1*np.cos(angle)+.002
  muscle.addNewPathPoint('insertion',model.getGround(),osim.Vec3(float(L),0,0))
  model.addForce(muscle);model.finalizeConnections();state=model.initSystem()
  muscle.setActivation(state,activation);model.equilibrateMuscles(state);model.realizeDynamics(state)
  actual_force=float(muscle.getTendonForce(state));actual_length=float(muscle.getFiberLength(state))
  params=np.array([[1000.,.1,.2,.1*np.sin(angle),.033,.05]])
  bank=SeriesElasticBank(params);bank.initialize(np.array([L]),np.array([activation]))
  prototype_length=float(np.hypot(bank.state[0],params[0,3]))
  muscle.setFiberLength(state,prototype_length);model.realizeDynamics(state)
  same_state_tendon=float(muscle.getTendonForce(state))
  static_fiber=(1000.*activation*muscle.getActiveForceLengthMultiplier(state)+muscle.getPassiveFiberForce(state))*muscle.getCosPennationAngle(state)
  row=dict(activation=activation,pennation_opt_rad=angle,opensim_force_N=actual_force,prototype_force_N=float(bank.force[0]),force_difference_N=float(bank.force[0]-actual_force),fiber_length_difference_m=prototype_length-actual_length,
      same_fiber_state_tendon_error_N=float(bank.force[0]-same_state_tendon),opensim_static_balance_at_prototype_state_N=float(static_fiber-same_state_tendon))
  rows.append(row)
report=dict(opensim_version=osim.GetVersion(),scope='Static isometric equilibrium only; no dynamic force-velocity or integration equivalence.',cases=rows,max_force_error_N=max(abs(r['force_difference_N']) for r in rows),max_fiber_length_error_m=max(abs(r['fiber_length_difference_m']) for r in rows))
report['maximum_same_state_tendon_error_N']=max(abs(r['same_fiber_state_tendon_error_N']) for r in rows)
report['maximum_opensim_static_balance_error_N']=max(abs(r['opensim_static_balance_at_prototype_state_N']) for r in rows)
# OpenSim's default initialization is approximate; compare the actual curves
# at identical fiber states more tightly than that initialization difference.
assert report['max_force_error_N']<1.,report
assert report['maximum_same_state_tendon_error_N']<1e-6,report
assert report['maximum_opensim_static_balance_error_N']<1e-6,report
(OUT/'OPENSIM_STATIC_CHECK.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
