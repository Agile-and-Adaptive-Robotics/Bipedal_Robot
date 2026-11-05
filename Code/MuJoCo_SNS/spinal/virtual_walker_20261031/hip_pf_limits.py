"""Soft thigh-excursion sensory feedback through existing inhibitory PF relays.

Adds neural current, never clips joint states or actuator forces. Engineering
tilt/proprioceptive signal, not a claimed direct muscle-spindle measurement.
"""
import numpy as np
class HipPFLimits:
 def __init__(self,c,cfg):
  self.c=c;self.cfg=dict(cfg);assert set(cfg)<=set(('lower_deg','upper_deg','gain_nA_deg','limit_nA'));self.base=c.tonic.copy();self.targets=[];self.bindings=[];self.last=np.zeros(6)
  assert -20<=cfg['lower_deg']<cfg['upper_deg']<=25 and 0<cfg['gain_nA_deg']<=1 and 0<cfg.get('limit_nA',12.)<=12
  for side in ('L','R'):
   for phase in ('flx','ext'):
    target=c.names.index(side+' Hip PF '+phase);edges=[e for e in c.edges if int(e[1])==target and c.names[int(e[0])].startswith(side+' Hip PF ')and c.names[int(e[0])].endswith(' IN')and e[5]<c.params[target,0]]
    assert len(edges)==1,'Resolve the actual inhibitory relay from wiring, not its label';source=int(edges[0][0]);self.targets.append(source);self.bindings.append(dict(sensor_side=side,suppressed_PF=c.names[target],excited_IN=c.names[source],gmax_uS=float(edges[0][2]),reversal_mV=float(edges[0][5])))
 def observe(self,m,d,dt):
  angles=[];currents=[]
  for side in ('L','R'):
   h=m.joint('hip_'+side).id;k=m.joint('knee_'+side).id;v=d.xanchor[k]-d.xanchor[h];angle=float(np.degrees(np.arctan2(v[0],-v[1])));angles.append(angle)
   currents.extend((max(0.,angle-self.cfg['upper_deg']),max(0.,self.cfg['lower_deg']-angle)))
  current=np.minimum(self.cfg.get('limit_nA',12.),self.cfg['gain_nA_deg']*np.array(currents));self.c.tonic[self.targets]=self.base[self.targets]+current;self.last[:]=[*angles,*current]
