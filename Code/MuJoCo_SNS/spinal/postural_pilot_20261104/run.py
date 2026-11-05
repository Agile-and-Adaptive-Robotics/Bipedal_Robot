"""Bounded neural equivalence and MuJoCo ankle-bench verification, not walking."""
import argparse,json,os,time
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import sys
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np
import mujoco
from model import build,XML,NAMES,INPUTS


def electrical_check():
    result={}
    for architecture in ('premotor','direct_mn'):
        a=build(architecture);b=build(architecture,backend='native')
        rng=np.random.default_rng(61104);error=0.
        for i in range(1500):
            u=rng.uniform(0.,.5,len(INPUTS));u[:2]*=2.
            a.compiled.forward(u);b.compiled.forward(u)
            error=max(error,float(np.max(abs(a.compiled.V-b.compiled.V))))
        result[architecture]=dict(maximum_voltage_error_mV=error,pass_=error<1e-11)
    dc=[]
    for activation in (.05,.2,.5,.8):
        row={}
        for architecture in ('premotor','direct_mn'):
            n=build(architecture,ib_policy='off');u=n.u.copy();u[0]=n.command(activation)
            for _ in range(1500):n.compiled.forward(u)
            row[architecture]=float(np.clip(n.compiled.V[n.idx['MN_POS']]/5.,0.,1.))
        row['command']=activation;row['difference']=abs(row['premotor']-row['direct_mn']);dc.append(row)
    result['gain_match']=dc
    result['passed']=all(result[a]['pass_']for a in ('premotor','direct_mn')) and max(r['difference']for r in dc)<1e-7
    return result


def trial(architecture,ib_policy,kp,kd,initial=.04,noise=0.,seed=61104,push=8.,parameters=None,backend='compact',delay_s=0.):
    n=build(architecture,ib_policy,backend=backend,parameters=parameters);m=mujoco.MjModel.from_xml_string(XML);d=mujoco.MjData(m)
    d.qpos[0]=0.;d.act[:]=.05;mujoco.mj_forward(m,d)
    arms=-d.actuator_moment[:,0];order=np.array([int(np.argmax(arms)),int(np.argmin(arms))])
    if arms[order[0]]<=0 or arms[order[1]]>=0:raise ValueError('No physical antagonistic torque pair')
    length0=d.actuator_length.copy();capacity=300.*float(np.min(abs(arms)))
    d.qpos[0]=initial;mujoco.mj_forward(m,d)
    count=8000;log=np.empty((count//5,9));neur=np.empty((count//5,len(NAMES)))
    dt=.001;rng=np.random.default_rng(seed);eta=np.zeros(2);tick=time.perf_counter();used=0;guard=False
    if not 0<=delay_s<=.12:raise ValueError('Delay test must be between 0 and 120 ms')
    delay_steps=round(delay_s/dt);history=np.empty((delay_steps+1,8));history[:]=np.r_[d.qpos[0],d.qvel[0],d.actuator_length,d.actuator_velocity,d.actuator_force]
    for i in range(count):
        t=i*dt;d.qfrc_applied[0]=push if 2.<=t<2.05 else 0.
        row=history[i%len(history)];row[0]=d.qpos[0];row[1]=d.qvel[0]
        row[2:4]=d.actuator_length;row[4:6]=d.actuator_velocity;row[6:8]=d.actuator_force
        sensor=history[(i-delay_steps)%len(history)]
        torque=-kp*sensor[0]-kd*sensor[1]
        commands=np.clip(np.array([max(torque,0.),max(-torque,0.)])/capacity,0.,1.)
        n.u[:]=0.;n.u[0:2]=n.command(commands)
        eta*=np.exp(-dt/.01);eta+=noise*np.sqrt(1.-np.exp(-2*dt/.01))*rng.normal(size=2)
        n.u[:2]+=eta
        for side,a in zip(('POS','NEG'),order):
            n.u[n.ports['Ia_'+side]]=1.6*np.clip(sensor[4+a]/.1,0.,1.)
            n.u[n.ports['II_'+side]]=1.2*np.clip((sensor[2+a]-length0[a])/.02,0.,1.)
            n.u[n.ports['Ib_'+side]]=1.5*np.clip(-sensor[6+a]/300.,0.,1.)
        n.compiled.forward(n.u)
        out=np.clip(n.compiled.V[[n.idx['MN_POS'],n.idx['MN_NEG']]]/5.,0.,1.)
        d.ctrl[order]=.05+.95*out
        mujoco.mj_step(m,d)
        if i%5==0:
            log[used]=[t,d.qpos[0],d.qvel[0],*d.ctrl[order],*-d.actuator_force[order],torque,d.qfrc_applied[0]]
            neur[used]=n.compiled.V;used+=1
        if not np.isfinite(d.qpos).all() or abs(d.qpos[0])>.4:
            guard=True;break
    log=log[:used];neur=neur[:used];tail=log[:,0]>max(3.,log[-1,0]-2.)
    rms=float(np.degrees(np.sqrt(np.mean(log[tail,1]**2)))) if tail.any()else 1e6
    settled=np.flatnonzero((log[:,0]>2.05)&(abs(log[:,1])<np.radians(1.))&(abs(log[:,2])<np.radians(5.)))
    settle=None
    for k in settled:
        window=log[k:k+100]
        if len(window)==100 and np.all(abs(window[:,1])<np.radians(1.))and np.all(abs(window[:,2])<np.radians(5.)):
            settle=float(log[k,0]-2.05);break
    result=dict(architecture=architecture,ib_policy=ib_policy,kp=kp,kd=kd,initial_rad=initial,noise_nA=noise,seed=seed,push_Nm=push,parameters=parameters or {},backend=backend,added_sensor_delay_s=delay_s,
        guard=guard,simulated_s=float(d.time),wall_s=time.perf_counter()-tick,
        maximum_angle_deg=float(np.degrees(abs(log[:,1]).max())),final_rms_deg=rms,
        settling_after_push_s=settle,mean_squared_excitation=float(np.mean(log[:,3:5]**2)),accepted=bool(not guard and rms<1. and settle is not None and settle<3.),
        scope='Fixed-base 10 kg single-ankle MuJoCo bench with synthetic tension actuators; not Gait2392, BPA pressure control, or walking')
    return result,log,neur


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--kp',type=float,default=80.);p.add_argument('--kd',type=float,default=10.);p.add_argument('--noise',type=float,default=0.);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False)
    check=electrical_check();(a.output/'ELECTRICAL_CHECK.json').write_text(json.dumps(check,indent=2))
    if not check['passed']:raise RuntimeError('Reference/gain gate failed')
    rows=[]
    for architecture in ('direct_mn','premotor'):
        for ib in ('off','inhibitory','support_positive'):
            result,log,neur=trial(architecture,ib,a.kp,a.kd,noise=a.noise)
            tag=architecture+'_'+ib
            np.savez_compressed(a.output/(tag+'.npz'),mechanics=log,neural=neur,names=np.array(NAMES))
            rows.append(result);(a.output/'RESULTS.json').write_text(json.dumps(rows,indent=2))
            print(json.dumps(result),flush=True)
    (a.output/'model.xml').write_text(XML)
    (a.output/'CONFIG.json').write_text(json.dumps(dict(frame='x anterior,y superior,z right',neural_dt_s=.001,physics_dt_s=.001,initial_rad=.04,perturbation='8 Nm for50ms at2s',activation_floor=.05,kp=a.kp,kd=a.kd,noise_nA=a.noise),indent=2))


if __name__=='__main__':main()
