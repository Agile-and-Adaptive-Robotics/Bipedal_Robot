"""Opt-in all-external CAD force bridge; no reference-force/BPA double counting."""
import numpy as np
import mujoco
from cad_interface import RoutedBPABank

class CADMovingForce:
 def __init__(self,model,dynamic,names):
  self.model=model;self.dynamic=dynamic;self.names=list(names)
  # All native force channels must be off, not just the selected subset.
  if np.any(model.actuator_gainprm) or np.any(model.actuator_biasprm):raise ValueError('Disable all native reference gains/biases explicitly before applying external CAD tensions')
  self.aids=np.array([model.actuator(n).id for n in names]);self.tids=model.actuator_trnid[self.aids,0]
  if np.any(model.actuator_trntype[self.aids]!=mujoco.mjtTrn.mjTRN_TENDON) or not np.allclose(model.actuator_gear[self.aids,0],1.):raise ValueError('This bridge supports unit-gear native tendon transmissions only')
  self.last_generalized_force=np.zeros(model.nv)
 def update_geometry(self,data):
  self.dynamic.update(data);return data.ten_length[self.tids].copy(),self.dynamic.velocities[self.tids].copy()
 def apply_tensions(self,data,positive_tension):
  """Add once per physics tick after caller reset; never call native correction."""
  tension=np.asarray(positive_tension)
  if tension.shape!=(len(self.names),) or not np.isfinite(tension).all() or np.any(tension<0):raise ValueError('Require finite nonnegative tensile forces')
  self.last_generalized_force[:]=RoutedBPABank.generalized_force(self.dynamic.length_jacobian[self.tids],tension);data.qfrc_applied[:]+=self.last_generalized_force
  return self.last_generalized_force.copy()
 def advance_pressure(self,data,bank,normalized_neural_output):
  if bank.names!=self.names:raise ValueError('BPA bank/geometry identity mismatch')
  lengths,velocity=self.update_geometry(data);tension=bank.step(normalized_neural_output,lengths);force=self.apply_tensions(data,tension)
  return dict(lengths=lengths,velocities=velocity,pressure_kpa=bank.pressure.copy(),positive_tensions=tension.copy(),generalized_force=force)
