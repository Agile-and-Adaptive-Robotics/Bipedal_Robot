"""Optional force/current recorder; separate MjData leaves the live model unchanged."""
import numpy as np,mujoco

class SwingDiagnostics:
    def __init__(self,net,capacity=350,stride=25):
        self.net=net;self.stride=stride;self.used=0;self.probe=None;self.capacity=capacity
        c=net.compiled;self.mn=np.array([net.idx['MN_'+a] for a in net.targets_data['act_names']]);self.to_act={int(p):i for i,p in enumerate(self.mn)}
        self.edges=net.mn_edges;self.dest=np.array([self.to_act[int(c.post[e])] for e in self.edges])
        def group(name):
            if name.startswith('PF_'):return 0
            if name.startswith('IaIN_'):return 2
            if name.startswith(('Ia_','IBIN_','IBEXC_')):return 1
            if name.startswith(('POST','BAL','PT_')):return 3
            return 4
        self.group=np.array([group(net.names[c.pre[e]]) for e in self.edges]);self.currents=np.zeros((capacity,5,len(self.mn)))
        self.voltage=np.zeros((capacity,len(c.V)));self.t=np.zeros(capacity);self.leak=np.zeros((capacity,len(self.mn)))
        self.qpos=[];self.qvel=[];self.act=[];self.ctrl=[];self.force=[];self.passive=[];self.rho=[];self.terms=[]
        native=c.forward
        def forward(u):
            if net.live is not None and net.step_count%stride==0 and self.used<capacity:
                k=self.used;v=c.V;ids=self.edges
                values=c.gmax[ids]*np.clip((v[c.pre[ids]]-c.elo[ids])/c.erange[ids],0.,1.)*(c.reversal[ids]-v[c.post[ids]])
                for j in range(5):self.currents[k,j]=np.bincount(self.dest[self.group==j],weights=values[self.group==j],minlength=len(self.mn))
                self.voltage[k]=v;self.leak[k]=-c.c.g_m[self.mn]*(v[self.mn]-c.c.V_rest[self.mn])+c.c.i_b[self.mn]
            return native(u)
        c.forward=forward
    def sample(self,time):
        n=self.net
        # step_count has incremented after the recorded neural update.
        if (n.step_count-1)%self.stride or self.used>=self.capacity:return
        m,d=n.live
        if self.probe is None:self.probe=mujoco.MjData(m)
        p=self.probe;p.qpos[:]=d.qpos;p.qvel[:]=d.qvel;p.act[:]=0.;mujoco.mj_forward(m,p)
        rho=n.arm_model.arms.copy();f=np.maximum(-d.actuator_force,0.);f0=np.maximum(-p.actuator_force,0.)
        self.t[self.used]=time;self.qpos.append(d.qpos.copy());self.qvel.append(d.qvel.copy());self.act.append(d.act.copy());self.ctrl.append(d.ctrl.copy())
        self.force.append(f);self.passive.append(f0);self.rho.append(rho)
        self.terms.append(np.stack([rho@(f-f0),rho@f0,n.arm_model.E.T@d.qfrc_passive,n.arm_model.E.T@d.qfrc_constraint,n.arm_model.E.T@d.qfrc_bias,d.qacc[n.arm_model.rows]]))
        self.used+=1
    def save(self,path):
        n=self.net;k=self.used
        np.savez_compressed(path,t=self.t[:k],currents=self.currents[:k],voltage=self.voltage[:k],leak_and_bias=self.leak[:k],
            names=n.names[:len(n.compiled.V)],acts=n.targets_data['act_names'],joints=n.arm_model.jids,qpos=self.qpos,qvel=self.qvel,act=self.act,ctrl=self.ctrl,
            force=self.force,passive_force=self.passive,rho=self.rho,terms=self.terms,
            current_groups=['PF','Ia_Ib','IaIN','Posture_balance','other'],torque_terms=['active_muscle','passive_muscle','nonmuscle_passive','constraint','bias','joint_acceleration'],
            scope='Pre-update synaptic currents; mechanical state at that same physics instant. Probe MjData has activation zero. Passive includes native muscle bias; joint springs/damping are separate. No live state or control clamp.')
