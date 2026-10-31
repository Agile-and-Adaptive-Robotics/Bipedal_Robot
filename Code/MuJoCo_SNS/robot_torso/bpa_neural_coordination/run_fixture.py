"""Bounded SNS/BPA actuator-interface experiment, not a full walker controller."""
from __future__ import annotations
import os,sys,time,json,hashlib,math
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np
import mujoco
from sns_toolbox.networks import Network
from sns_toolbox.neurons import NonSpikingNeuron
from sns_toolbox.connections import NonSpikingSynapse
from actuator_interface import packed_parameters,force_at_pressure,inverse_pressure,pressure_step,neural_step

HERE=Path(__file__).resolve().parent
OUT=HERE/'runs'/'fixture_20261031_v1'
DT=.002;CENTER=.35;AMPLITUDE=.12;FREQUENCY=.5
NAMES=['PF_F','PF_E','II_positive','II_negative','Ia_positive','Ia_negative','IN_F','IN_E','MN_F','MN_E']
XML='''<mujoco model="BPA_neural_interface_fixture"><compiler angle="radian"/><option timestep="0.002" integrator="implicitfast" gravity="0 -9.81 0"/>
<worldbody><body name="shank"><joint name="knee_fixture" axis="0 0 1" range="0.05 1.1" limited="true" damping="0.08" armature="0.05" solreflimit="0.004 1" solimplimit="0.999 0.9999 0.001"/>
<geom type="capsule" fromto="0 0 0 0 -0.3 0" size="0.025" mass="1"/></body></worldbody>
<tendon><fixed name="flexor"><joint joint="knee_fixture" coef="-0.03"/></fixed><fixed name="extensor"><joint joint="knee_fixture" coef="0.03"/></fixed></tendon></mujoco>'''

def network(cfg):
    net=Network(name='BPA interface fixture: prescribed PF, neural proprioception')
    cell=NonSpikingNeuron(membrane_capacitance=10.,membrane_conductance=1.,resting_potential=0.)
    for name in NAMES:net.add_neuron(cell,name=name)
    for name in NAMES[:6]:net.add_input(name,name=name+'_current')
    for name in NAMES:net.add_output(name,name=name+'_voltage')
    edges=[]
    def connect(src,dst,g,e=20.):
        net.add_connection(NonSpikingSynapse(max_conductance=g,reversal_potential=e,e_lo=0.,e_hi=5.),src,dst)
        edges.append(dict(source=src,destination=dst,conductance_uS=g,reversal_mV=e))
    connect('PF_F','MN_F',.1);connect('PF_E','MN_E',.1)
    for s,m in [('II_positive','MN_E'),('II_negative','MN_F')]:connect(s,m,cfg['kp'])
    for s,m in [('Ia_positive','MN_E'),('Ia_negative','MN_F')]:connect(s,m,cfg['kd'])
    connect('MN_F','IN_F',.03);connect('MN_E','IN_E',.03)
    connect('IN_F','MN_E',.03,-5.);connect('IN_E','MN_F',.03,-5.)
    compiled=net.compile(dt=DT*1000.,backend='numpy')
    assert not any((compiled.spiking,compiled.gated,compiled.electrical,compiled.electrical_rectified))
    arrays=tuple(np.asarray(getattr(compiled,x)) for x in ('g_max_non','e_lo','e_hi','del_e','time_factor_membrane','g_m','V_rest','i_b','input_connectivity'))
    return compiled,arrays,edges

def run(cfg,plant,feedback=True,backend='numba',duration=8.,save=None):
    compiled,arrays,edges=network(cfg)
    model=mujoco.MjModel.from_xml_string(XML);data=mujoco.MjData(model)
    data.qpos[0]=CENTER;mujoco.mj_forward(model,data)
    p=packed_parameters();v=compiled.V.copy();scratch=np.empty_like(v);inputs=np.zeros(6)
    pressure=np.zeros(2);history=np.zeros((2,128));h_state=np.zeros(2)
    traces=np.empty((round(duration/DT),36));delay_steps=round(plant['delay']/DT)
    assert delay_steps<128
    start=time.perf_counter();neural_wall=0.
    for k in range(traces.shape[0]):
        t=(k+1)*DT;phase=2*np.pi*FREQUENCY*t
        # Descending current is prescribed. No independent force/torque controller.
        wave=cfg['amp']*math.sin(phase+cfg['phase'])
        inputs[0]=max(0.,cfg['bias']+wave+cfg.get('offset',0.))
        inputs[1]=max(0.,cfg['bias']-wave-cfg.get('offset',0.))
        q=data.qpos[0];vel=data.qvel[0]
        pos=np.clip((q-CENTER)/.2,-2.,2.);speed=np.clip(vel/1.,-2.,2.)
        inputs[2]=max(pos,0.) if feedback else 0.;inputs[3]=max(-pos,0.) if feedback else 0.
        inputs[4]=max(speed,0.) if feedback else 0.;inputs[5]=max(-speed,0.) if feedback else 0.
        tick=time.perf_counter()
        if backend=='numba':neural_step(v,scratch,inputs,*arrays)
        else:v[:]=compiled.forward(inputs)
        neural_wall+=time.perf_counter()-tick
        target=np.clip(v[8:10]/5.,0.,1.)*300.
        lengths=np.array([.33-.03*(q-CENTER),.33+.03*(q-CENTER)])
        speeds=np.array([-.03*vel,.03*vel]);commands=np.zeros(2);forces=np.zeros(2);flags=np.zeros(2)
        for j in range(2):
            commands[j],flags[j]=inverse_pressure(target[j],lengths[j],p)
            history[j,k%128]=commands[j]
            delayed=history[j,(k-delay_steps)%128] if k>=delay_steps else 0.
            pressure[j]=pressure_step(pressure[j],delayed,DT,plant['fill'],plant['vent'])
            h_state[j]+=(math.tanh(speeds[j]/.004)-h_state[j])*(1.-math.exp(-DT/.08))
            forces[j]=force_at_pressure(lengths[j],pressure[j],p)*plant['force_gain']*(1.+plant['direction_gain']*h_state[j])
        torque=.03*(forces[0]-forces[1])
        # Known held-out external pulse; no force allocation or actuator bypass.
        disturbance=plant.get('disturbance',0.) if 4.0<=t<4.12 else 0.
        data.qfrc_applied[0]=torque+disturbance
        mujoco.mj_step(model,data)
        traces[k]=np.r_[t,data.qpos[0],data.qvel[0],CENTER+AMPLITUDE*math.sin(phase),v,inputs,commands,pressure,target,forces,lengths,speeds,flags,torque,disturbance]
        if not np.all(np.isfinite(traces[k])):raise RuntimeError('Nonfinite fixture state')
    wall=time.perf_counter()-start
    use=traces[:,0]>=2.
    error=traces[use,1]-traces[use,3]
    result=dict(configuration=cfg,plant=plant,feedback=feedback,backend=backend,duration_s=duration,
        wall_s=wall,neural_wall_s=neural_wall,rms_angle_error_deg=float(np.degrees(np.sqrt(np.mean(error**2)))),
        min_angle_deg=float(np.degrees(traces[:,1].min())),max_angle_deg=float(np.degrees(traces[:,1].max())),
        no_hyperextension=bool(traces[:,1].min()>=0.),max_pressure_kpa=float(traces[:,22:24].max()),
        saturation_steps=int(np.count_nonzero(traces[:,32:34])),
        force_tracking_rms_N=float(np.sqrt(np.mean((traces[use,24:26]-traces[use,26:28])**2))),
        perturbation_peak_error_deg=float(np.degrees(np.max(np.abs(traces[(traces[:,0]>=4)&(traces[:,0]<=5),1]-traces[(traces[:,0]>=4)&(traces[:,0]<=5),3])))))
    if save:
        np.savez_compressed(save.with_suffix('.npz'),trace=traces,names=np.array(['time_s','q_rad','qdot_rad_s','target_q_rad']+NAMES+['PF_F_input','PF_E_input','II_positive_input','II_negative_input','Ia_positive_input','Ia_negative_input','command_F_kpa','command_E_kpa','actual_F_kpa','actual_E_kpa','target_F_N','target_E_N','actual_F_N','actual_E_N','length_F_m','length_E_m','velocity_F_m_s','velocity_E_m_s','inverse_F_flag','inverse_E_flag','actuator_torque_Nm','external_torque_Nm']))
        save.with_suffix('.json').write_text(json.dumps(result,indent=2))
    return result,traces,edges

def checks():
    sys.path.insert(0,str(HERE.parents[1]));from bpa_muscle import BPAMuscle,FESTO_COEFFS
    exported=json.loads((HERE.parents[2]/'Matlab/Functions/festo_lookup_coeffs.json').read_text())
    coefficient_error=max(abs(exported[name]['coeffs'][key]-FESTO_COEFFS[dia][i]) for dia,name in [(10,'f_10'),(20,'f20'),(40,'f40')] for i,key in enumerate(('a0','a1','a3')))
    errors=[];inverse_errors=[];count=0
    for dia in (10,20,40):
        p=packed_parameters(dia)
        reference=BPAMuscle('fixture',dia,.3,.225,fitting_length=.01,tendon_length=.04,series_stiffness=2e5)
        for length in np.linspace(.29,.355,8):
            for pressure in np.linspace(0,620,9):
                reference.set_activation(pressure/620.);actual=force_at_pressure(length,pressure,p)
                errors.append(abs(actual-reference.force(length)))
                if actual>1e-6:
                    recovered,flag=inverse_pressure(actual,length,p);inverse_errors.append(abs(recovered-pressure));assert flag==0
                count+=1
    assert max(errors)<2e-5 and max(inverse_errors)<2e-5 and coefficient_error==0
    _,flag=inverse_pressure(1e7,.33,packed_parameters());assert flag==2
    cfg=dict(kp=.1,kd=.1,amp=.4,bias=1.6,phase=0.)
    native,arrays,_=network(cfg);v=native.V.copy();scratch=np.empty_like(v);maximum=0.;rng=np.random.default_rng(3101)
    inputs=rng.uniform(0.,3.,(2000,6));tick=time.perf_counter()
    for row in inputs:
        neural_step(v,scratch,row,*arrays);expected=native.forward(row);maximum=max(maximum,float(np.max(np.abs(v-expected))))
    assert maximum<1e-12
    # Separate warmed kernel timings, identical input replay and initial state.
    v[:]=0.;tick=time.perf_counter()
    for row in inputs:neural_step(v,scratch,row,*arrays)
    fast_time=time.perf_counter()-tick
    native.reset();tick=time.perf_counter()
    for row in inputs:native.forward(row)
    reference_time=time.perf_counter()-tick
    return dict(static_cases=count,max_static_port_force_error_N=max(errors),max_inverse_pressure_error_kpa=max(inverse_errors),
        source_coefficient_max_error=coefficient_error,infeasible_force_saturates=True,neural_replay_max_error_mV=maximum,
        neural_kernel_native_s=reference_time,neural_kernel_numba_s=fast_time,neural_kernel_speedup=reference_time/fast_time)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    check=checks();(OUT/'CHECKS.json').write_text(json.dumps(check,indent=2));print('CHECKS',json.dumps(check),flush=True)
    nominal=dict(fill=.04,vent=.06,delay=.004,force_gain=1.,direction_gain=0.)
    train_slow=dict(fill=.065,vent=.08,delay=.008,force_gain=.95,direction_gain=0.)
    # A small declared fit set; held-out profiles below are never used here.
    rng=np.random.default_rng(3102);candidates=[]
    for i in range(18):
        cfg=dict(kp=float(rng.uniform(.035,.18)),kd=float(rng.uniform(.02,.12)),amp=float(rng.uniform(.15,.65)),bias=1.6,phase=float(rng.uniform(.05,.65)),offset=float(rng.uniform(.0,.08)))
        scores=[]
        for plant in (nominal,train_slow):
            r,_,_=run(cfg,plant);scores.append(r['rms_angle_error_deg']+(0 if r['no_hyperextension'] else 1000))
        candidates.append(dict(configuration=cfg,training_mean_error_deg=float(np.mean(scores))))
        (OUT/'TRAINING.json').write_text(json.dumps(candidates,indent=2));print('FIT',i+1,round(candidates[-1]['training_mean_error_deg'],4),flush=True)
    winner=min(candidates,key=lambda a:a['training_mean_error_deg']);(OUT/'FROZEN_CONFIGURATION.json').write_text(json.dumps(winner,indent=2))
    heldout={
        'nominal':nominal,
        'heldout_slow_weak':dict(fill=.10,vent=.12,delay=.016,force_gain=.82,direction_gain=.08),
        'heldout_fast_strong':dict(fill=.025,vent=.035,delay=.002,force_gain=1.15,direction_gain=-.08),
        'heldout_push':dict(fill=.075,vent=.10,delay=.012,force_gain=.9,direction_gain=.06,disturbance=.18)}
    results=[]
    for label,plant in heldout.items():
        for feedback in (True,False):
            r,traces,edges=run(winner['configuration'],plant,feedback,save=OUT/(label+('_feedback' if feedback else '_ablated')))
            r['case']=label;results.append(r);print('TEST',label,feedback,round(r['rms_angle_error_deg'],4),flush=True)
    r_native,tr_native,_=run(winner['configuration'],nominal,backend='native',save=OUT/'nominal_native')
    r_fast,tr_fast,_=run(winner['configuration'],nominal,save=OUT/'nominal_repeat')
    r_first=np.load(OUT/'nominal_feedback.npz')['trace']
    verification=dict(max_complete_trace_native_numba_difference=float(np.max(np.abs(tr_native-tr_fast))),
        max_fresh_model_repeat_difference=float(np.max(np.abs(r_first-tr_fast))),
        native_complete_wall_s=r_native['wall_s'],numba_complete_wall_s=r_fast['wall_s'],
        complete_speedup=r_native['wall_s']/r_fast['wall_s'],compilation_included=False,
        note='Single-joint fixture only; complete timing includes logging and per-step timing instrumentation.')
    assert verification['max_complete_trace_native_numba_difference']<1e-9
    assert verification['max_fresh_model_repeat_difference']==0.
    (OUT/'NETWORK.json').write_text(json.dumps(dict(names=NAMES,connections=edges,description='Prescribed alternating PF current; no RG/CPG; genuine ordinary MN synaptic control'),indent=2))
    (OUT/'RESULTS.json').write_text(json.dumps(dict(checks=check,training_winner=winner,results=results,verification=verification,
        source_hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),HERE/'actuator_interface.py',HERE.parents[1]/'bpa_muscle.py']},
        scope='One fixed-base joint; not CAD-robot standing/stepping or validated BPA dynamics'),indent=2))
    (OUT/'fixture.xml').write_text(XML)
    print('FINAL',json.dumps(verification),flush=True)

if __name__=='__main__':main()
