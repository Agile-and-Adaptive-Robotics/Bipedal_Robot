"""Report the measured command-to-torque chain, without inferring steps from bursts."""
from pathlib import Path
import json,numpy as np
H=Path(__file__).resolve().parent;O=H/'runs/swing_release_20261030_v1'
report={}
for case in ('air','ground'):
    z=np.load(O/(case+'.mechanics.npz'));names=z['names'].tolist();acts=z['acts'].tolist();t=z['t'];v=z['voltage'];late=(t>=5)&(t<13)
    rows={}
    for side,offset in [('r',0),('l',7)]:
        swing=(v[:,names.index('RG_F_'+side)]>v[:,names.index('RG_E_'+side)])&late
        stance=(~swing)&late
        out={}
        for state,keep in [('swing',swing),('stance',stance)]:
            muscle_torque=z['rho'][keep,:, :]*(z['force'][keep]-z['passive_force'][keep])[:,None,:]
            contributors={}
            for axis,row in [('hip',offset),('knee',offset+3),('ankle',offset+4)]:
                tau=muscle_torque[:,row,:].mean(axis=0);order=np.argsort(-abs(tau))[:8]
                contributors[axis]=[{ 'muscle':acts[i],'active_Nm':float(tau[i]),'activation':float(z['act'][keep,i].mean()),
                    'PF_current':float(z['currents'][keep,0,i].mean()),'reflex_current':float(z['currents'][keep,1:3,i].sum(axis=1).mean())} for i in order]
            out[state]=dict(samples=int(keep.sum()),active_knee_Nm=float(z['terms'][keep,0,offset+3].mean()),
                muscle_passive_knee_Nm=float(z['terms'][keep,1,offset+3].mean()),other_passive_knee_Nm=float(z['terms'][keep,2,offset+3].mean()),
                constraint_knee_Nm=float(z['terms'][keep,3,offset+3].mean()),bias_knee_Nm=float(z['terms'][keep,4,offset+3].mean()),contributors=contributors)
        rows[side]=out
    report[case]=rows
target=np.load(H/'runs/targets/targets.npz');strength=target['rho']*(target['unit_force']*target['activation'])[:,None,:]
target_summary={}
for row,label in [(0,'hip_r'),(3,'knee_r'),(7,'hip_l'),(10,'knee_l')]:
    pos=np.maximum(strength[:,row],0).sum(axis=1);neg=np.maximum(-strength[:,row],0).sum(axis=1)
    target_summary[label]=dict(positive_mean_Nm=float(pos.mean()),negative_mean_Nm=float(neg.mean()),
        opposing_fraction=float(np.mean(np.minimum(pos,neg)/np.maximum(np.maximum(pos,neg),1e-6))),net_mean_abs_Nm=float(np.mean(abs(pos-neg))))
report['training_torques']=target_summary
(O/'force_chain.json').write_text(json.dumps(report,indent=2))
for case in ('air','ground'):
    for side in ('r','l'):
        row=report[case][side]['swing'];print(case,side,'swing knee Nm',row['active_knee_Nm'],'passive',row['muscle_passive_knee_Nm'],row['other_passive_knee_Nm'])
        for a in row['contributors']['knee'][:5]:print(a)
print('TARGETS',json.dumps(target_summary))
