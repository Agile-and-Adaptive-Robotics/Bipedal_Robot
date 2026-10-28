"""Gait2392 task coordinates and equality-aware, preallocated muscle Jacobians."""
import numpy as np
from numba import njit
import mujoco

JOINTS=tuple(n+'_'+s for s in ('r','l') for n in
    ('hip_flexion','hip_adduction','hip_rotation','knee_angle','ankle_angle','subtalar_angle','mtp_angle'))+(
    'lumbar_extension','lumbar_bending','lumbar_rotation')
BI_FAMILIES={'HAM':('semimem','semiten','bifemlh'),'RECT':('rect_fem',),
             'SAR':('sar',),'GRAC':('grac',),'TFL':('tfl',),'GAS':('med_gas','lat_gas')}

@njit(cache=True)
def fill_projection(E,driver_rows,dependent,independent,coefficients,offsets,qpos,qadr):
    E[:]=0.
    for k in range(len(driver_rows)):E[driver_rows[k],k]=1.
    for i in range(len(dependent)):
        x=qpos[qadr[i]]-offsets[i]
        derivative=coefficients[i,1]+x*(2*coefficients[i,2]+x*(3*coefficients[i,3]+4*x*coefficients[i,4]))
        for j in range(E.shape[1]):E[dependent[i],j]=derivative*E[independent[i],j]

class CoupledArms:
    def __init__(self,model):
        self.model=model
        self.jids=np.array([mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,n) for n in JOINTS])
        assert np.all(self.jids>=0)
        self.rows=model.jnt_dofadr[self.jids].copy();self.addresses=model.jnt_qposadr[self.jids].copy()
        links=[]
        for eq in range(model.neq):
            if model.eq_type[eq]!=mujoco.mjtEq.mjEQ_JOINT:continue
            dep,ind=int(model.eq_obj1id[eq]),int(model.eq_obj2id[eq])
            if ind<0:continue
            links.append((int(model.jnt_dofadr[dep]),int(model.jnt_dofadr[ind]),model.eq_data[eq,:5].copy(),
                          float(model.qpos0[model.jnt_qposadr[ind]]),int(model.jnt_qposadr[ind])))
        # Dependents can themselves drive other followers. Resolve that order.
        ordered=[];known=set(self.rows.tolist());remaining=list(links)
        while remaining:
            ready=[r for r in remaining if r[1] in known]
            if not ready:
                # Equality branches unrelated to the selected coordinates have zero derivative.
                external=[r for r in remaining if r[1] not in {x[0] for x in remaining}]
                if not external:raise ValueError('Cyclic equality dependency')
                known.update(r[1] for r in external);continue
            for r in ready:ordered.append(r);known.add(r[0]);remaining.remove(r)
        self.dependent=np.array([r[0] for r in ordered],dtype=np.int64)
        self.independent=np.array([r[1] for r in ordered],dtype=np.int64)
        self.coefficients=np.array([r[2] for r in ordered]);self.offsets=np.array([r[3] for r in ordered])
        self.qadr=np.array([r[4] for r in ordered],dtype=np.int64)
        self.E=np.zeros((model.nv,len(JOINTS)));self.arms=np.zeros((len(JOINTS),model.nu))

    def update(self,data):
        fill_projection(self.E,self.rows,self.dependent,self.independent,self.coefficients,self.offsets,data.qpos,self.qadr)
        # Muscle force is negative along MuJoCo's dl/dq transmission.
        np.matmul(self.E.T,-data.actuator_moment.T,out=self.arms)
        return self.arms

    def project(self,data):
        for dep,ind,c,offset,qa in zip(self.dependent,self.independent,self.coefficients,self.offsets,self.qadr):
            jid=int(np.flatnonzero(self.model.jnt_dofadr==dep)[0]);ad=self.model.jnt_qposadr[jid]
            x=data.qpos[qa]-offset
            data.qpos[ad]=self.model.qpos0[ad]+np.polynomial.polynomial.polyval(x,c)

    def finite_difference(self,data,h=1e-5):
        saved=data.qpos.copy();result=np.zeros_like(self.arms)
        for k,adr in enumerate(self.addresses):
            for sign in (-1.,1.):
                data.qpos[:]=saved;data.qpos[adr]+=sign*h;self.project(data);mujoco.mj_forward(self.model,data)
                result[k]-=sign*data.actuator_length/(2*h)
        data.qpos[:]=saved;mujoco.mj_forward(self.model,data)
        return result
