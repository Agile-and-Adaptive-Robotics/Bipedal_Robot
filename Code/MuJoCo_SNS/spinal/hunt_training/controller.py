"""Opt-in Hunt-trained, full-3D Gait2392 SNS; no runtime MN target clamp.

Intrinsic and afferent parameters are separately fitted. Phase-filter neurons
supply a causal basis for joint PF fitting. Joint PF/MN and reflex weights
are fitted in later blocks. Joint-position feedback is an explicit engineering
posture pathway; native muscle afferents provide the fitted timing/reflex input.
"""
import os,sys,json
from pathlib import Path
from dataclasses import dataclass,field
H=Path(__file__).resolve().parent;S=H.parent
sys.path[:0]=[str(S/'functional_gait2392'),str(S),str(S/'w2l_cpg')]
import numpy as np
import params as P,build_network as BN,build_network_w2lvar as W
import network as F,deng2023_joint_pf as D
from sns_toolbox.neurons import NonSpikingNeuron
from sparse_pairwise_eval import SparsePairwiseForward
from mechanics import JOINTS

@dataclass
class HuntNetwork(D.JointIntegrationNetwork):
    hunt_pairs:dict=field(default_factory=dict)
    bank_names:dict=field(default_factory=dict)
    bank_connections:list=field(default_factory=list)

    def _build_pf_w2l(self,n,side):
        super()._build_pf_w2l(n,side)
        self._pair(n,'MTP_Z',side)
        n.add_input('RG_E_'+side);self.inputs.append('HUNT_INIT_E_'+side)
        n.add_input('RG_F_'+side);self.inputs.append('HUNT_INIT_F_'+side)
        for i,target in enumerate(('F','E','F','E','E','E')):
            name=f'HUNT_AF_{i}_{side}';self._add(name,.02,n);n.add_input(name);self.inputs.append(name)
            n.add_connection(BN._syn(self.options['sensory'][side][i],True),name,f'RG_{target}_{side}')
        banks=[]
        for phase in ('E','F'):
            for k,tau in enumerate((.02,.05,.1,.2,.4,.7)):
                name=f'HUNT_BANK_{phase}_{k}_{side}';self._add(name,tau,n)
                n.add_connection(BN._syn(1.,True),'RG_'+phase+'_'+side,name);banks.append(name)
        self.bank_names[side]=banks
        self.hunt_pairs.update({
            'hip_flexion_'+side:('PF_HIP-E_'+side,'PF_HIP-F_'+side),
            'hip_adduction_'+side:tuple(self.extra_pairs['HIP_X_'+side]),
            'hip_rotation_'+side:tuple(self.extra_pairs['HIP_Y_'+side]),
            'knee_angle_'+side:('PF_KNEE-E_'+side,'PF_KNEE-F_'+side),
            'ankle_angle_'+side:tuple(reversed(self.extra_pairs['ANKLE_Z_'+side])),
            'subtalar_angle_'+side:tuple(self.extra_pairs['FOOT_X_'+side]),
            'mtp_angle_'+side:tuple(self.extra_pairs['MTP_Z_'+side])})
        combinations={
            'SAR':(('PF_HIP-F_'+side,'PF_KNEE-F_'+side),('PF_HIP-E_'+side,'PF_KNEE-E_'+side)),
            'GRAC':((self.extra_pairs['HIP_X_'+side][1],'PF_KNEE-F_'+side),(self.extra_pairs['HIP_X_'+side][0],'PF_KNEE-E_'+side)),
            'TFL':(('PF_HIP-F_'+side,self.extra_pairs['HIP_X_'+side][0],'PF_KNEE-E_'+side),('PF_HIP-E_'+side,self.extra_pairs['HIP_X_'+side][1],'PF_KNEE-F_'+side))}
        for family,channels in combinations.items():
            for cell,sources in zip(self.extra_pairs['BI_'+family+'_'+side],channels):
                for source in sources:n.add_connection(BN._syn(1./len(sources),True),source,cell)
        if side==self.sides[0]:
            self.hunt_pairs.update(lumbar_extension=tuple(self.extra_pairs['LUM_Z_mid']),
                lumbar_bending=tuple(self.extra_pairs['LUM_X_mid']),lumbar_rotation=tuple(self.extra_pairs['LUM_Y_mid']))
            n.add_neuron(NonSpikingNeuron(membrane_capacitance=.05,membrane_conductance=1.,resting_potential=0.,bias=5.),name='PF_HUNT_TONIC')
            self.idx['PF_HUNT_TONIC']=len(self.idx)
        # Unique excitatory and inhibitory relay populations allow both signed
        # temporal basis effects without duplicate synapses on the same endpoints.
        for k,bank in enumerate(banks):
            relay=f'HUNT_BANK_IN_{k}_{side}';self._add(relay,.02,n)
            n.add_connection(BN._syn(1.,True),bank,relay)
        current=[(key,pair) for key,pair in self.hunt_pairs.items() if key.endswith('_'+side) or (key.startswith('lumbar') and side==self.sides[0])]
        for key,pair in current:
            for cell in pair:
                for k,bank in enumerate(banks):
                    for source,exc in ((bank,True),(f'HUNT_BANK_IN_{k}_{side}',False)):
                        n.add_connection(BN._syn(1e-12,exc),source,cell)
                        self.bank_connections.append((source,cell,key))

    def _wire_muscle_w2l(self,n,act,mi):
        super()._wire_muscle_w2l(n,act,mi)
        row=(0 if mi.side=='r' else 7)+6;key='MTP_Z_'+mi.side
        for sign,cell in zip((1.,-1.),self.extra_pairs[key]):
            self._edge(n,cell,'MN_'+act,1e-12,record=(act,row,sign,key))
        self._edge(n,'PF_HUNT_TONIC','MN_'+act,1e-12)

    def install(self):
        self.compiled=SparsePairwiseForward(self.compiled);c=self.compiled
        self.names=list(self.idx);self.ports={name:i for i,name in enumerate(self.inputs)}
        # Tune inherited joint-PF inhibition only in explicit Hunt candidates.
        scale=float(self.options.get('joint_pf_inhibition_scale',1.))
        if not 0.<scale<=1.:raise ValueError('joint_pf_inhibition_scale must be in (0, 1]')
        if scale!=1.:
            for i,(pre,post) in enumerate(zip(c.pre,c.post)):
                if self.names[pre].startswith(('PF_IN_HIP-','PF_IN_KNEE-')) and self.names[post].startswith(('PF_HIP-','PF_KNEE-')):
                    c.gmax[i]*=scale
        self.targets_data=np.load(H/'runs/targets/targets.npz',allow_pickle=True)
        self.target_meta=json.loads((H/'runs/targets/metadata.json').read_text())
        self.sensor_mid=np.asarray(self.target_meta['sensor_mid']);self.sensor_half=np.asarray(self.target_meta['sensor_half'])
        self.target_fmax=self.targets_data['Fmax'];self.geom_scale=np.maximum(np.max(abs(self.targets_data['rho']),axis=(0,2)),.003)
        self.mn_edges=np.array([i for i,p in enumerate(c.post) if self.names[p].startswith('MN_')],dtype=int)
        self.original_g=c.gmax.copy();self.pf_edges=[];self.reflex_edges=[];self.stand_edges=[]
        for i in self.mn_edges:
            pre=self.names[c.pre[i]]
            if pre.startswith('PF_') and not pre.startswith('PF_IN'):self.pf_edges.append(i)
            elif pre.startswith(('Ia_','IBIN_','IBEXC_')) or (self.options.get('learn_reciprocal',False) and pre.startswith('IaIN_')):self.reflex_edges.append(i)
            elif pre.startswith(('POST','BAL')):self.stand_edges.append(i)
        self.pf_edges=np.asarray(self.pf_edges,dtype=int);self.reflex_edges=np.asarray(self.reflex_edges,dtype=int)
        self.stand_edges=np.asarray(self.stand_edges,dtype=int)
        self.geometry={}
        for source,target,act,row,sign,key in self.extra_edges:
            edge=np.flatnonzero((c.pre==self.idx[source])&(c.post==self.idx[target]))
            assert len(edge)==1
            self.geometry[int(edge[0])]=(self.targets_data['act_names'].tolist().index(act),row,sign)
        self.geom_edges=np.array(list(self.geometry),dtype=int)
        self.geom_acts=np.array([self.geometry[i][0] for i in self.geom_edges]);self.geom_rows=np.array([self.geometry[i][1] for i in self.geom_edges])
        self.geom_signs=np.array([self.geometry[i][2] for i in self.geom_edges])
        self.pfbank_edges=[]
        for source,target,key in self.bank_connections:
            ids=np.flatnonzero((c.pre==self.idx[source])&(c.post==self.idx[target]));assert len(ids)==1
            self.pfbank_edges.append(int(ids[0]))
        self.pfbank_edges=np.asarray(self.pfbank_edges,dtype=int)
        self.pf_internal_edges=np.array([i for i,(pre,post) in enumerate(zip(c.pre,c.post))
            if self.names[pre].startswith(('PF_IN_HIP-','PF_IN_KNEE-')) and self.names[post].startswith(('PF_HIP-','PF_KNEE-'))],dtype=int)
        self.direct_pf_edges=np.array([i for i,(pre,post) in enumerate(zip(c.pre,c.post))
            if self.names[pre].startswith('RG_') and self.names[post].startswith(('PF_HIP-','PF_KNEE-'))],dtype=int)
        self.ankle_rg_edges=np.array([i for i,(pre,post) in enumerate(zip(c.pre,c.post)) if self.names[pre].startswith('RG_') and self.names[post].startswith('PF_FUNC_ANKLE_Z_')],dtype=int)
        direct_gain=float(self.options.get('direct_joint_pf_gain',0.))
        if not 0.<=direct_gain<=20.:raise ValueError('direct_joint_pf_gain must be in [0, 20]')
        if direct_gain>0.:assert len(self.direct_pf_edges)==4*len(self.sides)
        # Remove inherited timing inputs: the six fitted sensors are the only
        # sensory projections into each calibrated RG; local INs are retained.
        for i in range(len(c.post)):
            pre,post=self.names[c.pre[i]],self.names[c.post[i]]
            crossed_contact=self.options.get('cross_contact_gain',0.)>0. and post.startswith('RG_F_') and pre=='HEEL_'+('l' if post.endswith('_r') else 'r')
            ipsi_heel=self.options.get('ipsilateral_heel_reset_gain',0.)>0. and post.startswith('RG_E_') and pre=='HEEL_'+post[-1]
            if post.startswith('RG_') and not pre.startswith(('InE_','InF_','HUNT_AF_','HUNT_COMM_')) and not crossed_contact and not ipsi_heel:c.gmax[i]=0.
            if post.startswith(('InE_','InF_')) and not pre.startswith('RG_'):c.gmax[i]=0.
            if post.startswith('PF_') and pre.startswith('RG_'):c.gmax[i]=0.
            if post.startswith('PF_') and pre.startswith(('HEEL_','TOEDF_')):c.gmax[i]=0.
        c.gmax[self.mn_edges]=0.
        self.pf_gain=np.zeros(len(c.gmax));self.ff_gain=np.zeros(len(c.gmax));self.reflex_gain=np.zeros(len(c.gmax))
        self.pf_gain[self.pf_internal_edges]=c.gmax[self.pf_internal_edges]
        self.geometry_factor=np.ones(len(c.gmax));self.step_count=0;self.training_gate=None
        self.reciprocal_edges=np.array([i for i in self.mn_edges if self.names[c.pre[i]].startswith('IaIN_')],dtype=int)
        self.bi_edges=np.array([i for i,(post,pre) in enumerate(zip(c.post,c.pre)) if self.names[post].startswith(('PF_BI_JOINT_HAM_','PF_BI_JOINT_RECT_')) and self.names[pre].startswith(('PF_HIP','PF_KNEE','PF_IN_BI_PAIR_'))],dtype=int)
        self.bi_input_edges=np.array([i for i in self.bi_edges if self.names[c.pre[i]].startswith(('PF_HIP','PF_KNEE'))],dtype=int)
        self.bi_gain=c.gmax.copy()
        self.load_weights()
        self.sync()

    def sync(self):
        c=self.compiled;c.c.g_max_non[c.post,c.pre]=c.gmax

    def load_weights(self):
        for filename,target in [('pf_weights.json',self.pf_gain),('bi_weights.json',self.bi_gain),('mn_weights.json',self.ff_gain),('reflex_weights.json',self.reflex_gain)]:
            path=Path(self.options.get('weights_directory',H/'runs/fitted'))/filename
            if not path.exists():continue
            values=json.loads(path.read_text())['edges']
            for i,(post,pre) in enumerate(zip(self.compiled.post,self.compiled.pre)):
                key=self.names[pre]+' -> '+self.names[post]
                if key in values:target[i]=float(values[key])

    def apply_weights(self,gate):
        c=self.compiled
        c.gmax[self.pfbank_edges]=self.pf_gain[self.pfbank_edges]*gate
        c.gmax[self.pf_internal_edges]=self.pf_gain[self.pf_internal_edges]
        if self.options.get('direct_joint_pf_gain',0.)>0.:
            c.gmax[self.direct_pf_edges]=self.options['direct_joint_pf_gain']*gate
        if self.options.get('direct_ankle_pf_gain',0.)>0.:
            c.gmax[self.ankle_rg_edges]=self.options['direct_ankle_pf_gain']*gate
        c.gmax[self.bi_edges]=self.bi_gain[self.bi_edges]
        c.gmax[self.mn_edges]=0.
        # Geometry-dependent recruitment stays active for positional feedback;
        # tonic/sagittal pattern expression is descending-command gated.
        c.gmax[self.pf_edges]=self.ff_gain[self.pf_edges]*gate
        c.gmax[self.geom_edges]=self.ff_gain[self.geom_edges]*self.geometry_factor[self.geom_edges]
        c.gmax[self.reflex_edges]=self.reflex_gain[self.reflex_edges]
        c.gmax[self.reciprocal_edges]=self.reflex_gain[self.reciprocal_edges] if self.options.get('learn_reciprocal',False) else self.original_g[self.reciprocal_edges]*self.options.get('reciprocal_gain',0.)
        c.gmax[self.stand_edges]=self.original_g[self.stand_edges]*(1.-gate)
        if hasattr(self,'balance'):self.balance.apply()
        self.sync()

    def bind_mujoco(self,model,data):
        F.FunctionalNetwork.bind_mujoco(self,model,data)
        self.goal[:]=np.mean(self.targets_data['qpos'][:,self.arm_model.addresses],axis=0)
        if 'hip_abduction_goal' in self.options:self.goal[[1,8]]=self.options['hip_abduction_goal']
        self.Fmax=np.maximum(model.actuator_gainprm[:,2],5.)
        if hasattr(self,'balance'):self.balance.bind(model,data)

    def timing_features(self,data,side):
        features=[];groups=self.target_meta['feature_members'][side]
        for i,indices in enumerate(groups):
            if i<2:raw=(data.actuator_length-self.sensor_mid)/self.sensor_half
            elif i<4:raw=data.actuator_velocity/.6
            else:raw=np.maximum(-data.actuator_force,0.)/self.Fmax
            features.append(float(np.maximum(raw[indices],0.).mean()))
        return np.clip(features,0.,1.5)

    def step(self,u):
        gate=float(np.clip(u[self.ports['DRIVE']]/2.5,0.,1.)) if self.training_gate is None else self.training_gate
        if self.live is not None:
            model,data=self.live
            if self.step_count%5==0:
                rho=self.arm_model.update(data)
                self.geometry_factor[self.geom_edges]=np.clip(self.geom_signs*rho[self.geom_rows,self.geom_acts]/self.geom_scale[self.geom_rows],0.,1.)
            for side in self.sides:
                for i,value in enumerate(self.timing_features(data,side)):
                    u[self.ports[f'HUNT_AF_{i}_{side}']]=5.*value*self.options.get('timing_scale',1.) if self.options['use_timing_afferents'] else 0.
            q=data.qpos[self.arm_model.addresses];v=data.qvel[self.arm_model.rows]
            demand=np.clip((self.goal-q)/.2-v/.8,-2.,2.)*self.options['axis_feedback']
            # Existing axis ports use FUNC_<label>_<side>_<channel>.
            for key,cells in self.extra_pairs.items():
                if key.startswith(('BI_','ANKLE_Z')):continue
                side=key[-1];offset=0 if side=='r' else 7
                row=15 if key.startswith('LUM_X') else 16 if key.startswith('LUM_Y') else 14 if key.startswith('LUM_Z') else offset+(1 if key.startswith('HIP_X') else 2 if key.startswith('HIP_Y') else 6 if key.startswith('MTP') else 5)
                for i,sign in enumerate((1.,-1.)):u[self.ports['FUNC_'+key+'_'+str(i)]]=max(0.,sign*demand[row])
        if hasattr(self,'balance'):self.balance.update(u)
        self.apply_weights(gate);self.compiled.forward(u);self.step_count+=1
        return self.compiled.V

def build(acts,dt=.002,interleg=False,options=None):
    options=dict(options or {})
    circuit_file=Path(options.get('weights_directory',H/'runs/fitted'))/'circuit_options.json'
    if circuit_file.exists():
        saved=json.loads(circuit_file.read_text())
        if set(saved)-{'joint_pf_inhibition_scale','direct_joint_pf_gain','require_bilateral_symmetry','direct_ankle_pf_gain'}:raise ValueError('Unknown saved Hunt circuit option')
        options={**saved,**options}
    a=json.loads((H/'runs/step_A/winner.json').read_text())
    sensory={side:json.loads((H/('runs/step_B' if side=='r' else 'runs/step_B_left')/'winner.json').read_text()) for side in ('r','l')}
    timing_file=Path((options or {}).get('weights_directory',H/'runs/fitted'))/'timing_weights.json'
    if timing_file.exists():
        timing=json.loads(timing_file.read_text())
        if options.get('require_bilateral_symmetry',False):
            from bilateral import assert_tied
            if timing['gains']['r']!=timing['gains']['l']:raise ValueError('Bilateral timing gains must be identical')
            for file in ('pf_weights.json','bi_weights.json','mn_weights.json','reflex_weights.json'):
                p=timing_file.parent/file
                if p.exists():assert_tied(json.loads(p.read_text())['edges'])
        if not timing.get('accepted'):raise RuntimeError('Candidate sensory fit has not passed expected-afferent gates')
        sensory={side:dict(accepted=True,afferent_gains=timing['gains'][side]) for side in ('r','l')}
    if not a['accepted'] or not all(v['accepted'] for v in sensory.values()):raise RuntimeError('Intrinsic/sensory fits must pass before motor fitting.')
    P.NAP.update(a['nap']);P.TAU['rg_nap_h']=a['tau_h'];P.G['rg_mutual_inh']=a['mutual_inhibition'];P.G['rg_weak_exc']=0.
    for key in ('w2l_ia_central','w2l_ii_central','w2l_ib_central'):W.GAINS[key]=0.
    for key,value in W.VARIANT_G.items():P.G[key]=float(value)
    opts={**F.DEFAULTS,'axes':True,'biarticular':True,'joint_mix_gain':1.,'axis_feedback':0.,'bi_afferent':0.,
        'sensory':{s:v['afferent_gains'] for s,v in sensory.items()},'use_timing_afferents':True,**(options or {})}
    routing=json.loads((F.OUT/'routing_audit.json').read_text())
    muscles={act:W.classify(act) for act in acts}
    net=HuntNetwork(muscles=muscles,sides=tuple(sorted({m.side for m in muscles.values()})),interleg=False,options=opts,routing=routing)
    if opts.get('commissural_gain',0.)>0.:
        for alias in net._aliases:net.idx.pop(alias,None)
        for side in ('r','l'):net._add('HUNT_COMM_'+side,.05,net.net)
        for side,other in (('r','l'),('l','r')):
            net.net.add_connection(BN._syn(1.,True),'RG_F_'+side,'HUNT_COMM_'+side)
            net.net.add_connection(BN._syn(opts['commissural_gain'],False),'HUNT_COMM_'+side,'RG_F_'+other)
        net._n_true=len(net.idx);net.idx.update(net._aliases)
    if opts.get('commissural_e_gain',0.)>0.:
        for alias in net._aliases:net.idx.pop(alias,None)
        for side in ('r','l'):net._add('HUNT_COMM_E_'+side,.05,net.net)
        for side,other in (('r','l'),('l','r')):
            net.net.add_connection(BN._syn(1.,True),'RG_E_'+side,'HUNT_COMM_E_'+side)
            net.net.add_connection(BN._syn(opts['commissural_e_gain'],False),'HUNT_COMM_E_'+side,'RG_E_'+other)
        net._n_true=len(net.idx);net.idx.update(net._aliases)
    if opts.get('balance'):
        from balance import add_populations
        add_populations(net)
    if opts.get('biarticular_pair_gain',0.)>0.:
        for alias in net._aliases:net.idx.pop(alias,None)
        for side in ('r','l'):
            for family,other in (('HAM','RECT'),('RECT','HAM')):
                relay=f'PF_IN_BI_PAIR_{family}_{side}';net._add(relay,P.TAU['pf'],net.net)
                net.net.add_connection(BN._syn(1.,True),'PF_BI_JOINT_'+family+'_'+side,relay)
                net.net.add_connection(BN._syn(opts['biarticular_pair_gain'],False),relay,'PF_BI_JOINT_'+other+'_'+side)
        net._n_true=len(net.idx);net.idx.update(net._aliases)
    if opts.get('cross_contact_gain',0.)>0.:
        for side,other in (('r','l'),('l','r')):net.net.add_connection(BN._syn(opts['cross_contact_gain'],True),'HEEL_'+side,'RG_F_'+other)
    if opts.get('ipsilateral_heel_reset_gain',0.)>0.:
        for side in ('r','l'):net.net.add_connection(BN._syn(opts['ipsilateral_heel_reset_gain'],True),'HEEL_'+side,'RG_E_'+side)
    net.compile(dt);net.install();net.fast_io=True
    if opts.get('balance'):
        from balance import TaskBalance
        net.balance=TaskBalance(net,opts['balance'])
    return net
