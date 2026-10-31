"""Experimental Peterka/Ting-inspired task feedback, not a paper replication.

Delayed CoM displacement/velocity and separately optional filtered acceleration
drive four directional SNS populations. Fixed, pose-calibrated recruitment
projects onto native MNs; inhibitory effects pass through relay neurons.
No runtime inverse-dynamics allocator or target-voltage cancellation is used.
"""
import numpy as np,mujoco
from scipy.optimize import lsq_linear
import build_network as BN

LABELS=('XP','XN','YP','YN') # Native world x anterior, y left.
def add_populations(net):
    # W2L's logging aliases are not populations; keep them after new true cells.
    for alias in net._aliases:net.idx.pop(alias,None)
    for label in LABELS:
        source='PT_BAL_'+label;relay='PT_IN_'+label
        net._add(source,.02,net.net);net.net.add_input(source);net.inputs.append(source)
        net._add(relay,.02,net.net);net.net.add_connection(BN._syn(1.,True),source,relay)
        for act in net.muscles:
            net.net.add_connection(BN._syn(1e-12,True),source,'MN_'+act)
            net.net.add_connection(BN._syn(1e-12,False),relay,'MN_'+act)
    net._n_true=len(net.idx);net.idx.update(net._aliases)

class TaskBalance:
    def __init__(self,net,options):
        self.n=net;self.options=options;self.dt=.002
        self.delay_steps=max(0,round(options.get('delay_s',.1)/self.dt))
        self.history=np.zeros((max(1,self.delay_steps),6));self.cursor=0;self.initialized=False
        self.last_position=None;self.velocity=np.zeros(2);self.acceleration=np.zeros(2)
        self.gain=np.zeros(len(net.compiled.gmax));self.last_task=np.zeros(6)
        self.last_command=np.zeros(2);self.ports=np.array([net.ports['PT_BAL_'+l] for l in LABELS])
    def bind(self,m,d):
        self.m=m;self.d=d;self.jac=np.zeros((3,m.nv))
        # Same talus reference as the inherited standing screen, captured once.
        feet=np.array([d.xpos[m.body('talus_'+s).id] for s in ('r','l')])
        self.reference=feet.mean(axis=0)[:2];self.reference[1]=d.subtree_com[0,1]
        probe=mujoco.MjData(m);probe.qpos[:]=d.qpos;probe.qvel[:]=0.;probe.act[:]=.05
        mujoco.mj_forward(m,probe);mujoco.mj_jacSubtreeCom(m,probe,self.jac,0)
        baseline=self.jac@probe.qacc;response=np.zeros((3,m.nu));epsilon=.02
        for i in range(m.nu):
            probe.act[:]=.05;probe.act[i]+=epsilon;mujoco.mj_forward(m,probe)
            response[:,i]=(self.jac@probe.qacc-baseline)/epsilon
        # Fixed nonnegative recruitment for +/- horizontal acceleration, with
        # vertical cross-effects discouraged and bounded per-muscle activity.
        matrix=np.vstack([response,np.eye(m.nu)*.1]);basis=[];fit=[]
        for axis,sign in ((0,1.),(0,-1.),(1,1.),(1,-1.)):
            wanted=np.zeros(3);wanted[axis]=sign
            solution=lsq_linear(matrix,np.r_[wanted,np.zeros(m.nu)],bounds=(0.,.3),method='bvls')
            basis.append(solution.x);fit.append((response@solution.x).tolist())
        if self.options.get('bilateral_recruitment',False):
            # Pool sagittal responses and mirror the jointly fitted lateral pair.
            acts=[m.actuator(i).name for i in range(m.nu)]
            mirror=np.array([acts.index(a[:-2]+('_l' if a.endswith('_r') else '_r'))
                             if a.endswith(('_r','_l')) else i for i,a in enumerate(acts)])
            families=list(dict.fromkeys(a[:-2] if a.endswith(('_r','_l')) else a for a in acts))
            tied=np.zeros((m.nu,len(families)))
            for i,a in enumerate(acts):tied[i,families.index(a[:-2] if a.endswith(('_r','_l')) else a)]=1.
            for k,sign in ((0,1.),(1,-1.)):
                wanted=np.array([sign,0.,0.])
                sol=lsq_linear(np.vstack([response@tied,.1*tied]),np.r_[wanted,np.zeros(m.nu)],bounds=(0.,.3),method='bvls')
                basis[k]=tied@sol.x
            sol=lsq_linear(np.vstack([response,response[:,mirror],.1*np.eye(m.nu)]),
                          np.r_[0.,1.,0.,0.,-1.,0.,np.zeros(m.nu)],bounds=(0.,.3),method='bvls')
            basis[2]=sol.x;basis[3]=sol.x[mirror]
            fit=[(response@b).tolist() for b in basis]
            if not all(np.array_equal(basis[k],basis[k][mirror]) for k in (0,1)) or not np.array_equal(basis[2],basis[3][mirror]):
                raise ValueError('Balance recruitment violates bilateral tying')
        self.basis=np.array(basis);self.response=response
        c=self.n.compiled;names=self.n.names;acts=self.n.targets_data['act_names'].tolist()
        for i,(pre,post) in enumerate(zip(c.pre,c.post)):
            source=names[pre];target=names[post]
            if not target.startswith('MN_'):continue
            act=acts.index(target[3:])
            if source.startswith('PT_BAL_'):
                k=LABELS.index(source[7:]);self.gain[i]=self.basis[k,act]*2.
            elif source.startswith('PT_IN_'):
                k=LABELS.index(source[6:]);opposite=k^1
                self.gain[i]=self.basis[opposite,act]*self.options.get('inhibition',.25)
        self.edges=np.flatnonzero(self.gain)
        self.calibration=dict(reference_native_xy_m=self.reference.tolist(),local_acceleration_response=fit,
            maximum_basis_activation=float(self.basis.max()),response_rank=int(np.linalg.matrix_rank(response)),
            bilateral_recruitment=bool(self.options.get('bilateral_recruitment',False)),
            basis_scope='Fixed-pose finite-activation forward-dynamics response; not a physiological synergy identification')
    def update(self,u):
        position=self.d.subtree_com[0,:2].copy()
        mujoco.mj_jacSubtreeCom(self.m,self.d,self.jac,0)
        actual_velocity=(self.jac@self.d.qvel)[:2]
        alpha=self.dt/(self.options.get('filter_s',.04)+self.dt)
        previous=self.velocity.copy();self.velocity+=alpha*(actual_velocity-self.velocity)
        self.acceleration+=alpha*((self.velocity-previous)/self.dt-self.acceleration)
        # Somatic estimate relative to the captured support reference and an
        # orientation-based inverted-pendulum proxy are explicitly weighted.
        up=self.d.xmat[self.m.body('torso').id].reshape(3,3)[:,1]
        height=max(self.d.subtree_com[0,2],.3)
        somatic=position-self.reference
        orientation=height*up[:2]/max(up[2],.3)
        weight=float(self.options.get('somatic_weight',1.))
        state=np.r_[weight*somatic+(1.-weight)*orientation,self.velocity,self.acceleration]
        if not self.initialized:self.history[:]=state;self.initialized=True
        delayed=self.history[self.cursor].copy();self.history[self.cursor]=state
        self.cursor=(self.cursor+1)%len(self.history)
        if self.delay_steps==0:delayed=state
        self.last_task[:]=delayed
        command=-(self.options.get('kp',0.)*delayed[:2]+self.options.get('kd',0.)*delayed[2:4]+self.options.get('ka',0.)*delayed[4:])
        self.last_command[:]=np.clip(command,-self.options.get('max_current',3.),self.options.get('max_current',3.))
        u[self.ports]=[max(self.last_command[0],0.),max(-self.last_command[0],0.),max(self.last_command[1],0.),max(-self.last_command[1],0.)]
        # Replace the inherited balance sensors to isolate the added pathway;
        # tonic standing posture and native muscle reflexes remain present.
        if self.options.get('replace_inherited',True):
            for name in ('BAL_PF','BAL_DF','BAL_TRK_EXT','BAL_TRK_FLX','BAL_LAT_R','BAL_LAT_L'):u[self.n.ports[name]]=0.
    def apply(self):
        self.n.compiled.gmax[self.edges]=self.gain[self.edges]
