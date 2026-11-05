"""Opt-in, load-confirmed heel reset to contralateral RG-F only.

Preserves ipsilateral extensor and PF contact projections. Event dwell is an
engineering handover condition, not an asserted neural conduction delay.
"""
import numpy as np
class HeelHandover:
 def __init__(self,c,cfg):
  self.c=c;self.cfg=dict(cfg);self.mode=cfg.get('mode','event');assert self.mode in('event','ablate','tonic_scale');assert set(cfg)<=set(('dwell_s','gain','mode'));assert 0<=cfg.get('dwell_s',0.)<=.2 and 0<cfg.get('gain',1.)<=4
  self.ids=[];self.base=[];self.bindings=[];self.filtered=np.zeros(2);self.dwell=np.zeros(2);self.unloaded=np.zeros(2);self.armed=np.zeros(2,dtype=bool);self.pulse=np.zeros(2);self.last=np.zeros(8)
  for side,other in(('L','R'),('R','L')):
   candidates=[i for i,e in enumerate(c.edges)if c.names[int(e[0])]==side+' heel contact'and c.names[int(e[1])]==other+' RG flx'];assert len(candidates)==1;idx=candidates[0];assert c.edges[idx,5]>c.params[int(c.edges[idx,1]),0];self.ids.append(idx);self.base.append(c.edges[idx,2]);self.bindings.append(dict(source=side+' heel contact',target=other+' RG flx',gmax_uS=float(c.edges[idx,2]),reversal_mV=float(c.edges[idx,5])))
  self.base=np.array(self.base)
  if self.mode=='ablate':c.edges=np.delete(c.edges,self.ids,axis=0)
  elif self.mode=='tonic_scale':c.edges[self.ids,2]=self.base*cfg['gain']
  else:c.edges[self.ids,2]=0.
 def observe(self,m,d,measured,dt):
  if self.mode!='event':return
  bw=float(m.body_subtreemass[m.body('Root').id])*9.81;load=np.zeros(2)if measured is None else np.maximum(measured['foot_force_world'][:,1],0.)/bw
  self.filtered+=(1-np.exp(-dt/.02))*(load-self.filtered);self.pulse*=np.exp(-dt/.05)
  for i in range(2):
   if self.filtered[i]<.05:
    self.unloaded[i]+=dt;self.dwell[i]=0.
    if self.unloaded[i]>=.08:self.armed[i]=True
   elif self.filtered[i]>.1:
    self.unloaded[i]=0.;self.dwell[i]+=dt
    if self.armed[i]and self.dwell[i]>=self.cfg['dwell_s']:
     self.pulse[i]=self.cfg['gain'];self.armed[i]=False
  self.pulse[self.pulse<1e-6]=0.;self.c.edges[self.ids,2]=self.base*self.pulse;self.last[:]=np.r_[self.filtered,self.dwell,self.pulse,self.armed.astype(float)]
