"""Reference-pose recruitment for Hunt's geometry-dependent PF readout.

The live controller updates these gates every five neural steps. This adapter
uses reference moment arms at those same update times; it does not advance body
physics. The original ungated replay remains available for archived diagnostics.
"""
import numpy as np
from replay_adapter import Replay

class ContractReplay(Replay):
    def __init__(self,manifest):
        super().__init__(manifest)
        self.pose_enabled=manifest.get('pose_recruitment_replay',False)
        self.feedforward_only=manifest.get('feedforward_only',False)
        if self.pose_enabled and manifest['family']!='hunt':raise ValueError('Only Hunt has this live pose-readout contract')
        z=np.load(manifest['targets'],allow_pickle=False);n=self.net
        if self.pose_enabled:
            query=manifest['target_start_s']+self.t[::5]%self.T
            self.pose_scale=np.empty((len(query),len(n.geom_edges)))
            for i,(row,act,sign)in enumerate(zip(n.geom_rows,n.geom_acts,n.geom_signs)):
                self.pose_scale[:,i]=np.clip(sign*np.interp(query,z['t'],z['rho'][:,row,act])/n.geom_scale[row],0.,1.)
        from controller_profiles_20261104.profile import central
        self.local_slots=np.array([i for i,(pre,post)in enumerate(zip(self.c.pre,self.c.post))if
            (self.names[pre].startswith(('Ia_','II_','Ib_'))and not central(self.names[post]))or
            (self.names[pre].startswith(('IBIN_','IBEXC_','IIX_','IIIN_'))and self.names[post].startswith('MN_'))],dtype=int)
        self.pose_step=0;self.pose_base=np.empty(len(n.geom_edges))if self.pose_enabled else None;forward=self.c.forward
        def pose_forward(u):
            if self.pose_step==0:
                if self.pose_enabled:self.pose_base[:]=self.c.gmax[n.geom_edges]
                if self.feedforward_only:self.c.gmax[self.local_slots]=0.
            if self.pose_enabled and self.pose_step%5==0:self.c.gmax[n.geom_edges]=self.pose_base*self.pose_scale[self.pose_step//5]
            result=forward(u);self.pose_step+=1;return result
        self.c.forward=pose_forward
    def evaluate(self,*args,**kwargs):
        self.pose_step=0
        return super().evaluate(*args,**kwargs)

# Compatibility for the explicitly recorded first pose-replay experiment.
PoseReplay=ContractReplay
