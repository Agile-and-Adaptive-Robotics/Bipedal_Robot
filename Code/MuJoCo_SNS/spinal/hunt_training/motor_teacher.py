"""Diagnostic motor teachers only. Their clock/feedback is never a neural gait result."""
import numpy as np,mujoco
from scipy.optimize import lsq_linear

class MotorTeacher:
    def __init__(self,n,options):
        self.n=n;self.options=options;self.z=n.targets_data;self.mn=np.array([n.idx['MN_'+a] for a in self.z['act_names']]);self.probe=None;self.a=None;self.log=[]
        self.target_t=self.z['t'];self.target_a=self.z['activation'];self.target_q=self.z['qpos'];self.target_v=self.z['qvel']
        self.samples=[]
    def update(self,time):
        if time<2:return
        n=self.n;m,d=n.live;period=self.options.get('period_s',1.233);query=.617+((time-2)%period)*1.233/period
        i=int(np.clip(np.searchsorted(self.target_t,query)-1,0,len(self.target_t)-2));fraction=(query-self.target_t[i])/(self.target_t[i+1]-self.target_t[i])
        desired=self.target_a[i]+fraction*(self.target_a[i+1]-self.target_a[i])
        if self.options['mode']=='feedback' and (self.a is None or (n.step_count-1)%self.options.get('stride',5)==0):
            if self.probe is None:self.probe=mujoco.MjData(m)
            p=self.probe;p.qpos[:]=d.qpos;p.qvel[:]=d.qvel;p.act[:]=0.;mujoco.mj_forward(m,p);f0=-p.actuator_force.copy()
            p.act[:]=1.;mujoco.mj_forward(m,p);capacity=np.maximum(-p.actuator_force-f0,0.)
            rho=n.arm_model.arms;addresses=n.arm_model.addresses;rows=n.arm_model.rows
            qr=(self.target_q[i]+fraction*(self.target_q[i+1]-self.target_q[i]))[addresses];vr=(self.target_v[i]+fraction*(self.target_v[i+1]-self.target_v[i]))[rows]*1.233/period
            qr[[1,8]]=.08;qr[[2,5,6,9,12,13]]=0.;vr[[1,2,5,6,8,9,12,13]]=0.
            kp=self.options.get('kp',180.);kd=self.options.get('kd',12.)
            tau=n.arm_model.E.T@(d.qfrc_bias-d.qfrc_passive-d.qfrc_constraint)+kp*(qr-d.qpos[addresses])+kd*(vr-d.qvel[rows])
            A=rho*capacity[None,:];rhs=tau-rho@f0;norm=np.linalg.norm(A,axis=0);lam=.02*np.median(norm[norm>1e-5])
            if self.options.get('sagittal_only',False):
                selected=[0,3,4,7,10,11];A=A[selected];rhs=rhs[selected]
            solution=lsq_linear(np.vstack([A,lam*np.eye(m.nu)]),np.r_[rhs,np.zeros(m.nu)],bounds=(0.,.98),method='bvls',tol=1e-6,max_iter=100)
            self.a=solution.x;self.log.append((time,float(np.linalg.norm(A@self.a-rhs)),float(np.sqrt(np.mean((qr-d.qpos[addresses])**2)))))
            self.samples.append(dict(t=time,qpos=d.qpos.copy(),qvel=d.qvel.copy(),activation=d.act.copy(),command=self.a.copy(),
                force=-d.actuator_force.copy(),passive=f0,capacity=capacity,rho=rho.copy(),
                features=np.stack([n.timing_features(d,side) for side in ('r','l')]),target_q=qr.copy()))
        if self.options['mode']=='feedback':desired=self.a
        n.compiled.c.V[self.mn]=5.*desired
    def save(self,path):
        if not self.samples:return
        np.savez_compressed(path,**{key:np.array([sample[key] for sample in self.samples]) for key in self.samples[0]},act_names=self.z['act_names'],period_s=self.options.get('period_s',1.233),scope='Clock/feedback teacher, not a neural controller. Training targets only; never a curriculum pass.')
