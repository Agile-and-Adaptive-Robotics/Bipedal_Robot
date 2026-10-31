"""Once-counted hindfoot/forefoot normal forces from actual contact locations.

The configurable split lies between calcaneus and toe body origins. These are
regional virtual sensors, not calibrated anatomical heel/toe sensor pads.
"""
import numpy as np
import mujoco

class RegionalContacts:
    def __init__(self,m,split=.5):
        if not 0.<split<1.:raise ValueError('Foot-region split must be between 0 and 1')
        self.m=m;self.split=split;self.force=np.zeros(6)
        self.calcn=[m.body('calcn_'+s).id for s in ('r','l')]
        self.toes=[m.body('toes_'+s).id for s in ('r','l')]
        self.owners={m.body(b+'_'+s).id:i for i,s in enumerate(('r','l')) for b in ('talus','calcn','toes')}
        self.loads=np.zeros((2,2))
    def sample(self,d):
        self.loads[:]=0.
        for k in range(d.ncon):
            c=d.contact[k];a,b=c.geom1,c.geom2
            if self.m.geom_type[a]==mujoco.mjtGeom.mjGEOM_PLANE:g=b
            elif self.m.geom_type[b]==mujoco.mjtGeom.mjGEOM_PLANE:g=a
            else:continue
            body=int(self.m.geom_bodyid[g])
            while body and body not in self.owners:body=int(self.m.body_parentid[body])
            if not body:continue
            side=self.owners[body];origin=d.xpos[self.calcn[side]]
            direction=d.xpos[self.toes[side]]-origin
            fraction=float((c.pos-origin)@direction/max(direction@direction,1e-12))
            region=int(body==self.toes[side] or fraction>=self.split)
            mujoco.mj_contactForce(self.m,d,k,self.force)
            self.loads[side,region]+=max(float(self.force[0]),0.)
        return self.loads
