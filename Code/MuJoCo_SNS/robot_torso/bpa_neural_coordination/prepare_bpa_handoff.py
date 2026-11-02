"""Adopt source-specific evidence, preserve superseded probes, and draw data plots."""
import json,xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from gait_standing_objective import score

ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa')
MODEL=Path(r'D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\mjc\gait2392_simbody\gait2392_simbody_cvt3.xml')

def comparison(a,b,fields):
 x=np.load(a);y=np.load(b);result={}
 for key in fields:
  p,q=x[key],y[key];result[key]=dict(shape_equal=p.shape==q.shape,bit_exact=p.shape==q.shape and np.array_equal(p,q),max_abs_difference=float(np.max(abs(p-q))) if p.shape==q.shape else None)
 return result

def main():
 nominal=ROOT/'pressure_standing/syn6_nominal_fresh_v1';fast=ROOT/'pressure_standing/syn6_nominal_fastio_v1'
 checks={}
 for suffix,fields in (('.npz',('qfull','ctrl','neuro')),('.apparatus.npz',('samples','neuro')),('.pressure.npz',('samples',))):
  checks.update({suffix+':'+k:v for k,v in comparison(nominal.with_suffix(suffix),fast.with_suffix(suffix),fields).items()})
 assert all(x['bit_exact'] for x in checks.values()),'Do not adopt fast motor IO until full trajectory parity passes'
 selections=[('nominal','pressure_standing/syn6_nominal_fastio_v1'),('hold20s','pressure_standing/syn6_nominal20s_exact_v1'),
 ('map_live','pressure_forward_recheck/syn6_maperror_live_equilibrium_v3'),('map_replay','pressure_forward_recheck/syn6_maperror_replay_equilibrium_v3'),
 ('delay_live','pressure_series_recheck/syn6_delay_live_series_v2'),('delay_replay','pressure_series_recheck/syn6_delay_replay_series_v2'),
 ('combined_live','pressure_forward_recheck/syn6_slowweak_live_equilibrium_v3'),('combined_replay','pressure_forward_recheck/syn6_slowweak_replay_equilibrium_v3'),
 ('push_live','pressure_heldouts/syn6_torso_push_live'),('push_replay','pressure_heldouts/syn6_torso_push_replay')]
 adopted=[];base_neuro=np.load(nominal.with_suffix('.apparatus.npz'))['neuro']
 for label,relative in selections:
  prefix=ROOT/relative;r=json.loads(prefix.with_suffix('.result.json').read_text());a=score(r,prefix)[1]
  z=np.load(prefix.with_suffix('.apparatus.npz'));n=len(z['neuro']);prefix_error=float(np.max(abs(z['neuro']-base_neuro[:n]))) if n<=len(base_neuro) else None
  extended=prefix.with_suffix('.apparatus_extended.npz');rig=None
  if extended.exists():
   e=np.load(extended);columns=e['columns'].tolist();values=e['samples']
   rig=dict(peak_absolute={c:float(np.max(abs(values[:,i]))) for i,c in enumerate(columns) if i>0},joint_parameters=json.loads(str(e['joint_parameters_json'])))
  adopted.append(dict(condition=label,result=str(prefix.with_suffix('.result.json')),assessment=a,completed_physics_samples=r['completed_samples'],failure=r['failure'],
   pressure_peak_kpa=r.get('pressure_peak_kpa'),inverse_saturation_fraction=r.get('inverse_saturation_fraction'),native_force_error_N=r.get('native_pressure_force_max_error_N'),
   upward_rig_bw=r.get('mean_upward_rig_support_bw'),downward_rig_bw=r.get('mean_downward_rig_force_bw'),net_rig_bw=r.get('mean_net_rig_force_bw'),mean_foot_load_bw=r.get('mean_ground_load_bw'),
   mesh_height_range_m=r['actual_visual_foot_mesh_clearance_range_m'],knee_max_deg=r['max_knee_angle_deg'],neural_prefix_change_from_nominal_max_mV=prefix_error,apparatus=rig))
 fresh=json.loads((ROOT/'BPA_HELDOUT_FRESH_CHECKS.json').read_text())
 report=dict(machine='easteregg2',controller_parameters=json.loads((ROOT/'rustuna_syn6_v1/BEST.json').read_text())['trial']['params'],
  fast_motor_io_full_trajectory_bit_exact=True,fast_motor_io_comparisons=checks,adopted_results=adopted,fresh_heldout_checks=fresh,
  superseded=['Initial uncertainty runs incorrectly checked unloaded total route length rather than sleeve length after series stretch.','Series-v2 force-error/composite exploratory rows applied force gain after equilibrium; their v3 rows solve the perturbed forward curve inside equilibrium. Delay-v2 has unit force gain and is unaffected.'],
  scope='Virtual Gait2392 standing-mode SNS/BPA comparison with declared geometry and pressure assumptions. Not mounted CAD hardware, independent balance, a Hunt curriculum pass, or locomotion. Live feedback uses existing myostatic/load pathways; BAL_* task channels are disabled by reflex-stand.')
 out=ROOT/'BPA_COMPARISON_ADOPTED.json';out.write_text(json.dumps(report,indent=2),encoding='utf-8')
 figures=ROOT/'figures';figures.mkdir(exist_ok=True)
 plt.rcParams.update({'font.family':'Arial','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':240})
 config=json.loads((ROOT/'virtual_bpa_design/VIRTUAL_GAIT_BPA_CONFIGURATION.json').read_text())
 actors=[e.attrib['name'] for e in ET.parse(MODEL).getroot().find('actuator')]
 z=np.load(nominal.with_suffix('.npz'));p=np.load(nominal.with_suffix('.pressure.npz'));a=np.load(nominal.with_suffix('.apparatus.npz'))
 channels=p['channels'].tolist();pid=p['names'].tolist().index('vas_med_l');time=p['time'];offset=len(z['ctrl'])-len(time);u=z['ctrl'][offset:]
 fig,axes=plt.subplots(4,1,figsize=(6.5,8.8),sharex=True,constrained_layout=True)
 axes[0].plot(time,u[:,actors.index('vas_med_l')],color='#0077BB',marker='o',markevery=45,ms=3,label='Vastus medialis, left')
 axes[0].plot(time,u[:,actors.index('semimem_l')],color='#EE7733',ls='--',marker='^',markevery=50,ms=3,label='Semimembranosus, left');axes[0].set_ylabel('MN command, u');axes[0].legend(loc='best')
 for name,label,color,style,mark in [('reference_pressure_kpa','Reference','#555555','--','s'),('actual_pressure_kpa','Actual','#0077BB','-','o')]:axes[1].plot(time,p['samples'][:,channels.index(name),pid],label=label,color=color,ls=style,marker=mark,markevery=45,ms=3)
 axes[1].set_ylabel('Pressure (kPa)');axes[1].legend(loc='best')
 axes[2].plot(time,u[:,actors.index('vas_med_l')]*config['actuators']['vas_med_l']['force_command_scale_N'],color='#555555',ls='--',marker='s',markevery=45,ms=3,label='Mapped neural reference')
 axes[2].plot(time,p['samples'][:,channels.index('tension_N'),pid],color='#0077BB',marker='o',markevery=50,ms=3,label='Actual BPA tension');axes[2].set_ylabel('Tension (N)');axes[2].legend(loc='best')
 for i,label,color,style,mark in [(1,'Right','#0077BB','-','o'),(2,'Left','#EE7733','--','^')]:axes[3].plot(time,np.degrees(a['samples'][:,i]),label=label,color=color,ls=style,marker=mark,markevery=45,ms=3)
 axes[3].axhline(0,color='#888888',ls=':',lw=1);axes[3].set_ylabel('Knee angle (°)');axes[3].legend(loc='best');axes[3].set_xlabel('Physics time after release (s)')
 for ax,title in zip(axes,['(A) Knee antagonist motor commands','(B) Left vastus medialis: pressure response','(C) Left vastus medialis: tension response','(D) Knee kinematics; negative is flexion']):ax.set_title(title,loc='left',fontsize=11);ax.set_xlim(0,.8);ax.grid(alpha=.18)
 for ext in ('png','pdf','svg'):fig.savefig(figures/('virtual_bpa_startup.'+ext))
 plt.close(fig)
 plt.rcParams.update({'font.size':17})
 fig,axes=plt.subplots(2,2,figsize=(11,6.8),sharex=True,constrained_layout=True)
 for col,(label,live,replay) in enumerate([('20% forward force-map loss',selections[2][1],selections[3][1]),('16 ms pressure delay',selections[4][1],selections[5][1])]):
  for relative,name,color,style,mark in [(selections[0][1],'Nominal','#999999',':','s'),(live,'Live SNS feedback','#0077BB','-','o'),(replay,'Nominal sensory replay','#EE7733','--','^')]:
   prefix=ROOT/relative;app=np.load(prefix.with_suffix('.apparatus.npz'));press=np.load(prefix.with_suffix('.pressure.npz'));tt=app['samples'][:,0]
   axes[0,col].plot(tt,np.degrees(app['samples'][:,2]),label=name,color=color,ls=style,marker=mark,markevery=25,ms=5)
   idx=press['names'].tolist().index('vas_med_l');ch=press['channels'].tolist().index('actual_pressure_kpa');axes[1,col].plot(press['time'],press['samples'][:,ch,idx],color=color,ls=style,marker=mark,markevery=25,ms=5)
   if name=='Nominal sensory replay':
    endpoint=float(np.degrees(app['samples'][-1,2]));axes[0,col].plot(tt[-1],endpoint,'kx',ms=7,mew=1.5)
    axes[0,col].annotate('Domain limit',xy=(tt[-1],endpoint),xytext=(tt[-1]+(.08 if col==0 else -.38),endpoint+2.5),arrowprops=dict(arrowstyle='->',lw=.8),fontsize=17)
  axes[0,col].set_title(f'({chr(65+col)}) '+label,loc='left');axes[0,col].set_ylabel('Left knee angle (°)');axes[1,col].set_ylabel('Pressure (kPa)');axes[1,col].set_xlabel('Physics time after release (s)')
  for row in range(2):axes[row,col].set_xlim(0,1.3);axes[row,col].grid(alpha=.18)
 axes[0,0].legend(fontsize=17,loc='upper right')
 for ext in ('png','pdf','svg'):fig.savefig(figures/('virtual_bpa_feedback_heldouts.'+ext))
 plt.close(fig)
 (figures/'virtual_bpa_startup.alt.txt').write_text('Four synchronized plots over0.8seconds after release show antagonist motor commands, reference and actual left-vastus pressure, mapped force reference and actual BPA tension, and bilateral knee flexion. Pressure rises from zero; the knees initially flex to about16degrees and recover near the protected2degree flexion margin. This is a supported virtual Gait2392 actuator-replacement simulation, not mounted CAD hardware.',encoding='utf-8')
 (figures/'virtual_bpa_feedback_heldouts.alt.txt').write_text('Two held-out conditions compare live SNS feedback with nominal sensory replay during virtual pressure-actuated standing. With20percent force-map loss or16millisecond pressure delay, live feedback returns the left knee toward the protected extension margin while replay drifts into the characterized equilibrium boundary. Black crosses mark the last recorded steps before termination, not observed falls. The lower plots show corresponding left-vastus pressures.',encoding='utf-8')
 (figures/'FIGURE_NOTES.md').write_text('# Virtual BPA figure provenance\n\nMachine: easteregg2. These figures show a separately fitted syn6 standing controller on Gait2392 with assumed virtual BPA dimensions and pressure dynamics, retained support springs, and measured pressure/tension. They are not calibrated CAD robot results.\n\nThe startup figure shows 0.8 s after release; the preceding 0.6 s pinned neural initialization is excluded. MN output is normalized excitation, the tension reference is its explicit force-scale conversion, and pressure dynamics create the realized force. Negative knee angle denotes flexion.\n\nThe held-out figure shows the initial loading transient. Black crosses mark the last valid samples before nominal sensory replay leaves the characterized pressure-supported series-equilibrium domain; they do not denote an observed fall. The replay preserves command biases and length-aware pressure inversion. Both live-feedback successes repeat bit-exactly in fresh processes.\n\nAll source paths, full apparatus loads, qualifications, and acceptance metrics: `'+str(out)+'`.\n',encoding='utf-8')
 print(json.dumps(dict(adopted_report=str(out),figures=str(figures),fast_motor_io_bit_exact=True,conditions=[dict(condition=x['condition'],accepted=x['assessment']['accepted']) for x in adopted]),indent=2))

if __name__=='__main__':main()
