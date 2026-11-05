"""Small graded antagonistic postural module for testing before deployment.

POS/NEG mean joint torque signs, not assumed anatomical flexor labels. Command
encoding matches the isolated static PN->MN transfer in the direct comparator.
Ib inhibition/support excitation are explicit alternatives, not contact rules.
"""
from dataclasses import dataclass
from pathlib import Path
import sys
import numpy as np

S=Path(__file__).resolve().parent.parent
sys.path[:0]=[str(S),str(S/'w2l_cpg')]
from sns_toolbox.networks import Network
from sns_toolbox.neurons import NonSpikingNeuron
from sns_toolbox.connections import NonSpikingSynapse
from compact_pairwise_eval import CompactPairwiseForward

R=5.; REV=8.; GM=1.
G_FEED=R*GM/(REV-R)
NAMES=('PP_POS','PP_NEG','PP_IN_POS','PP_IN_NEG','MN_POS','MN_NEG',
       'Ia_POS','Ia_NEG','II_POS','II_NEG','Ib_POS','Ib_NEG','IbRelay_POS','IbRelay_NEG')
INPUTS=('POST_POS','POST_NEG','Ia_POS','Ia_NEG','II_POS','II_NEG','Ib_POS','Ib_NEG')


def build(architecture='premotor',ib_policy='inhibitory',backend='compact',dt=.001,parameters=None):
    if architecture not in ('premotor','direct_mn') or ib_policy not in ('off','inhibitory','support_positive'):
        raise ValueError('Explicit architecture and Ib policy required')
    parameters=dict(parameters or {})
    allowed={'reciprocal_scale','Ia_scale','II_scale','Ib_scale','tau_scale'}
    if set(parameters)-allowed or any(not np.isfinite(v)or not 0<=v<=4 for v in parameters.values()):raise ValueError('Unknown or unsafe shared neural scale')
    n=Network(name='isolated antagonistic posture '+architecture)
    edges=[]
    for name in NAMES:
        tau=.01 if name.startswith('MN_') else .02
        if name.startswith('PP_'):tau*=parameters.get('tau_scale',1.)
        if tau<.005:raise ValueError('Premotor time constant is below the tested integration range')
        n.add_neuron(NonSpikingNeuron(membrane_capacitance=tau,membrane_conductance=GM,resting_potential=0.,bias=0.),name=name)
        n.add_output(name)
    for name in INPUTS:
        target=('PP_' if architecture=='premotor' else 'MN_')+name[5:] if name.startswith('POST_') else name
        n.add_input(target)
    def edge(source,target,g,exc=True):
        n.add_connection(NonSpikingSynapse(max_conductance=g,reversal_potential=REV if exc else -REV,e_lo=0.,e_hi=R),source,target)
        edges.append(dict(source=source,target=target,gmax_uS=g,reversal_mV=REV if exc else -REV))
    for side,other in (('POS','NEG'),('NEG','POS')):
        edge('PP_'+side,'PP_IN_'+side,G_FEED)
        edge('PP_IN_'+side,'PP_'+other,.45*parameters.get('reciprocal_scale',1.),False)
        edge('PP_'+side,'MN_'+side,G_FEED)
        edge('PP_IN_'+side,'MN_'+other,.45*parameters.get('reciprocal_scale',1.),False)
        target=('PP_' if architecture=='premotor' else 'MN_')+side
        edge('Ia_'+side,target,.03*parameters.get('Ia_scale',1.))
        edge('II_'+side,target,.04*parameters.get('II_scale',1.))
        # Preserve a small local monosynaptic Ia pathway in both cases.
        if architecture=='premotor':edge('Ia_'+side,'MN_'+side,.01)
        else:edges[-2]['local_Ia_combined']=True
        if architecture=='direct_mn':
            # Avoid duplicate endpoint overwrite in SNS-Toolbox: combine it.
            for c in n.connections:
                if n.populations[c['source']]['name']=='Ia_'+side and n.populations[c['destination']]['name']==target:
                    c['params']['max_conductance']+=.01
            for e in edges:
                if e['source']=='Ia_'+side and e['target']==target:e['gmax_uS']+=.01
        if ib_policy!='off':
            edge('Ib_'+side,'IbRelay_'+side,G_FEED)
            edge('IbRelay_'+side,target,.04*parameters.get('Ib_scale',1.),ib_policy=='support_positive')
    core=n.compile(dt=dt,backend='numpy')
    compiled=CompactPairwiseForward(core) if backend=='compact' else core
    return Circuit(compiled,architecture,ib_policy,edges)


@dataclass
class Circuit:
    compiled:object
    architecture:str
    ib_policy:str
    edges:list
    def __post_init__(self):
        self.idx={name:i for i,name in enumerate(NAMES)}
        self.ports={name:i for i,name in enumerate(INPUTS)}
        self.u=np.zeros(len(INPUTS))
    def command(self,activation):
        target=np.clip(activation,0.,1.)*R
        # Invert V_MN=Erev*g(V_PP)/(Gm+g(V_PP)) in the premotor case.
        voltage=target*GM*R/(G_FEED*(REV-target)) if self.architecture=='premotor' else target
        return GM*voltage


XML='''<mujoco model="postural_ankle_bench">
 <compiler angle="radian"/><option timestep="0.001" gravity="0 -9.81 0" integrator="implicitfast"/>
 <visual><global offwidth="640" offheight="480"/></visual>
 <worldbody>
  <light pos="0 1 1" dir="0 -1 -1"/>
  <camera name="side" pos="0 0.30 1.6" xyaxes="1 0 0 0 1 0" fovy="35"/>
  <geom name="ground" type="plane" size="2 2 .1" quat=".70710678 -.70710678 0 0" rgba=".9 .91 .93 1" contype="0" conaffinity="0"/>
  <geom type="box" pos="0 -.025 0" size=".19 .025 .09" rgba=".35 .4 .46 1" contype="0" conaffinity="0"/>
  <site name="back_anchor" pos="-.16 .035 0" size=".008"/>
  <site name="front_anchor" pos=".16 .035 0" size=".008"/>
  <body name="pendulum">
   <joint name="ankle" type="hinge" axis="0 0 1" limited="true" range="-.5 .5" damping=".2"/>
   <geom name="rod" type="capsule" fromto="0 0 0 0 .6 0" size=".035" mass="10" rgba=".15 .6 .65 1" contype="0" conaffinity="0"/>
   <site name="back_insert" pos="-.025 .26 0" size=".008"/>
   <site name="front_insert" pos=".025 .26 0" size=".008"/>
  </body>
 </worldbody>
 <tendon>
  <spatial name="back" width=".004" rgba=".8 .25 .15 1"><site site="back_anchor"/><site site="back_insert"/></spatial>
  <spatial name="front" width=".004" rgba=".2 .35 .8 1"><site site="front_anchor"/><site site="front_insert"/></spatial>
 </tendon>
 <actuator>
  <general name="back" tendon="back" dyntype="filter" dynprm=".03" gainprm="-300" ctrllimited="true" ctrlrange="0 1"/>
  <general name="front" tendon="front" dyntype="filter" dynprm=".03" gainprm="-300" ctrllimited="true" ctrlrange="0 1"/>
 </actuator>
</mujoco>'''
