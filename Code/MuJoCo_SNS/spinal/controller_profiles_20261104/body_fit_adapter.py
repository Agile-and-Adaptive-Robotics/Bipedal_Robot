"""Persistent Rustuna worker coordinating fresh-process body diagnostics."""
import hashlib,json,os,subprocess,sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
from dominant_cycles import measure,self_check
CONFIG=None

def handle_request(request):
    global CONFIG
    if request['op']=='preflight':
        CONFIG=request['manifest']
        for source,expected in CONFIG['source_sha256'].items():
            if hashlib.sha256(Path(source).read_bytes()).hexdigest()!=expected:raise ValueError('Source changed: '+source)
        prior=json.loads((Path(CONFIG['passed_air'])/'COMPLETE.json').read_text());neural=json.loads((Path(CONFIG['reflex_fit'])/'FINAL.json').read_text())
        ready=bool(prior[0]['accepted']and neural['accepted']and neural['fresh_process_matched']and self_check()['passed'])
        return dict(ok=True,ready_for_optimization=ready,scope='Prior deafferented physical and ordered neural gates; new body evaluations are diagnostic until both stages repeat')
    if CONFIG is None or request['op']!='evaluate':raise ValueError('Preflight required')
    params=request['params'];folder=Path(CONFIG['output'])/('trial_'+str(request['trial']));log=folder.with_suffix('.log')
    if folder.exists():raise ValueError('Preserve existing body evaluation')
    for source,expected in CONFIG['source_sha256'].items():
        if hashlib.sha256(Path(source).read_bytes()).hexdigest()!=expected:raise ValueError('Source changed during study: '+source)
    with log.open('a',encoding='utf-8')as stream:
        run=subprocess.run([sys.executable,str(HERE/'body_probe.py'),'--reflex-fit',CONFIG['reflex_fit'],'--stage','2','--expression',str(params['expression']),'--local-ii-scale',str(params['local_II_scale']),'--output',str(folder)],stdout=stream,stderr=subprocess.STDOUT,timeout=180,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
    path=folder/'body.result.json'
    if not path.exists():return dict(ok=True,value=10000.,diagnostics=dict(accepted=False,body_folder=str(folder),stop=json.loads((folder/'STOP.json').read_text())if(folder/'STOP.json').exists()else dict(exit_code=run.returncode)))
    result=json.loads(path.read_text());z=np.load(folder/'capture/aligned_trace.npz');t=z['t_state'];keep=(t>=5.)&(t<14.);names=z['joint_names'].tolist();q=z['qpos'][:,z['joint_qpos_addresses']]
    principal=measure(t[keep],-np.degrees(q[keep][:,[names.index('knee_angle_'+s)for s in('r','l')]]),np.degrees(q[keep][:,[names.index('hip_flexion_'+s)for s in('r','l')]]));old=result['assessment']
    error=abs(principal.get('period_s',0.)/1.233-1.)
    accepted=bool(principal['accepted']and old['finite']and old['complete']and old['frontal_gate']and old['knee_extension_gate']and error<=.15)
    knee=old.get('knee_ranges_deg',[0.,0.]);hip=old.get('hip_ranges_deg',[0.,0.])
    value=5.*error+5.*principal.get('bilateral_phase_error',.5)+.02*sum(abs(v-35.)for v in knee)
    if not accepted:value+=100.+2.*sum(max(20.-v,0.)for v in knee)+2.*sum(max(v-80.,0.)for v in hip)
    if principal.get('reason')=='dominant_frequency_outside_walking_window':value+=500.+50.*max(principal['dominant_Hz'])
    response=None
    if CONFIG.get('objective_profile')=='database_kinematics_v1':
        response=json.loads((folder/'body.kinematics.json').read_text())
        if response['schema']!='database_kinematics_v1':raise ValueError('Historical v1 objective requires its archived source snapshot')
        value=response['objective'];accepted=bool(accepted and response['accepted'])
    elif CONFIG.get('objective_profile')=='stage_specific_kinematics_v2':
        response=json.loads((folder/'body.kinematics.json').read_text())
        value=response['objective'];accepted=response['accepted']
    return dict(ok=True,value=float(value),diagnostics=dict(accepted=accepted,body_folder=str(folder),principal=principal,target_period_relative_error=error,legacy=old,kinematic_response=response,measured=result['measured'],wall_s=result['wall_s']))
