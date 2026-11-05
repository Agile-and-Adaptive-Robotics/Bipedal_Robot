"""Opt-in stance-gated pelvic posture command to antagonistic hip PF cells.

Explicit engineering PD task estimator; no direct MN drive or root force.
The existing PF pair provides the premotor channel. This is a diagnostic,
not the full added postural subnetwork proposed for Gait2392.
"""
import numpy as np
import mujoco
class PelvisPFPosture:
 def __init__(self,controller,configuration):
  self.c=controller;self.cfg=dict(configuration);self.base=controller.tonic.copy();self.ids=np.array([[controller.names.index(side+' Hip PF '+phase)for phase in ('flx','ext')]for side in ('L','R')]);self.gate=np.zeros(2);self.last=np.zeros(8)
  self.velocity=np.zeros(6)
  assert set(self.cfg)<=set(('target_deg','kp','kd','limit_nA','start_s'))
  assert 0<=self.cfg.get('kp',1.)<=8 and 0<=self.cfg.get('kd',.1)<=1 and 0<self.cfg.get('limit_nA',.3)<=1
 def observe(self,m,d,measured,t,dt):
  root=m.body('Root').id;R=d.xmat[root].reshape(3,3);angle=np.arctan2(R[0,1],R[1,1]);mujoco.mj_objectVelocity(m,d,mujoco.mjtObj.mjOBJ_BODY,root,self.velocity,0);updot=np.cross(self.velocity[:3],R[:,1]);rate=float((R[1,1]*updot[0]-R[0,1]*updot[1])/max(R[0,1]**2+R[1,1]**2,1e-9));target=np.radians(self.cfg.get('target_deg',4.));weight=float(m.body_subtreemass[root])*9.81
  loads=np.zeros(2)if measured is None else np.maximum(measured['foot_force_world'][:,1],0.)
  self.gate+=(1-np.exp(-dt/.02))*((loads>.05*weight)-self.gate)
  command=np.clip(self.cfg.get('kp',1.)*(target-angle)-self.cfg.get('kd',.1)*rate,-self.cfg.get('limit_nA',.3),self.cfg.get('limit_nA',.3))if t>=self.cfg.get('start_s',1.)else 0.
  commands=np.column_stack((np.maximum(command,0.)*self.gate,np.maximum(-command,0.)*self.gate))
  self.c.tonic[self.ids]=self.base[self.ids]+commands;self.last[:]=[np.degrees(angle),np.degrees(rate),*self.gate,*commands.ravel()]
