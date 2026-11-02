"""Read-only hierarchy/current comparison of same-controller air and ground."""
from pathlib import Path
import argparse,json
import numpy as np

def stats(x):
    return dict(minimum=float(x.min()),maximum=float(x.max()),range=float(np.ptp(x)),mean=float(x.mean()),std=float(x.std()))

def review(root,label):
    z=np.load(root/'configurations'/f'{label}.mechanics.npz',allow_pickle=False)
    result=json.loads((root/'configurations'/f'{label}.result.json').read_text())
    names=z['names'].tolist();t=z['t'];v=z['voltage'];use=(t>=5)&(t<=13)
    report={'label':label,'window_s':[float(t[use][0]),float(t[use][-1])],
            'full_voltage_samples':len(t),'full_neurons':len(names),'voltage_scope':'Actual pre-neural-update SwingDiagnostics samples, every25steps; not interpolated into post-physics capture.',
            'hierarchy':{},'HUNT_AF_to_RG_current':{},'HUNT_AF_relays':{}}
    select=[n for n in names if n.startswith(('RG_','PF_HIP','PF_KNEE','PF_BI_JOINT','PF_IN_BI_PAIR','HUNT_BANK_E_','HUNT_BANK_F_','HUNT_AF_'))]
    for name in select:report['hierarchy'][name]=stats(v[use,names.index(name)])
    timing=result['compiled']['integrated_fsa_design']['timing']
    currents={}
    for edge in timing['edges']:
        pre=v[:,names.index(edge['source'])];post=v[:,names.index(edge['target'])]
        activation=np.clip((pre-edge['new_threshold_mV'])/(edge['new_saturation_mV']-edge['new_threshold_mV']),0,1)
        current=edge['new_conductance_uS']*activation*(edge['reversal_mV']-post)
        currents.setdefault(edge['target'],np.zeros(len(t)));currents[edge['target']]+=current
        report['HUNT_AF_relays'][edge['source']]={**stats(pre[use]),'target':edge['target'],'conductance_uS':edge['new_conductance_uS'],'current_nA':stats(current[use])}
    report['HUNT_AF_to_RG_current']={key:stats(values[use])for key,values in currents.items()}
    report['HUNT_AF_current_scope']='Algebraic source-attributed conductance current from actual saved pre-update full voltages and this run\'s installed timing ledger. Separate from directly logged aggregate current; no inferred missing applied input.'
    a=np.load(root/'figures'/label/'aligned_trace.npz',allow_pickle=False);an=a['neural_names'].tolist();mask=(a['t_state']>=5)&(a['t_state']<=13)
    report['selected_voltage_and_current']={}
    for i,name in enumerate(an):
        row={'voltage_mV':stats(a['neural_voltage'][mask,i])}
        for part in('applied','chemical','leak','bias','ionic'):row[part+'_nA']=stats(a['neural_'+part+'_current_nA'][mask,i])
        row['maximum_update_error_mV']=float(abs(a['neural_update_residual_mV'][mask,i]).max())
        report['selected_voltage_and_current'][name]=row
    report['recorded_current_availability']={'all_samples_valid':bool(a['neural_current_available'].all()),'available':int(a['neural_current_available'].sum()),'total':len(a['t_state'])}
    # The legacy recorder already saved the complete pre-update MN-current
    # decomposition and route-arm/force arrays. Keep its native phase/time.
    import mujoco
    model=mujoco.MjModel.from_binary_path(str(root/'figures'/label/'model.mjb'))
    jointnames=[model.joint(int(j)).name for j in z['joints']];acts=z['acts'].tolist()
    groups=z['current_groups'].tolist();active_force=z['force']-z['passive_force']
    report['phase_current_and_torque']={}
    for side in('r','l'):
        fphase=v[:,names.index('RG_F_'+side)]>v[:,names.index('RG_E_'+side)]
        for phase,phase_mask in(('F',fphase),('E',~fphase)):
            mask=use&phase_mask;entry={'samples':int(mask.sum()),'MN_current_groups':{},'knee_torque_Nm':{}}
            for act in('vas_int_'+side,'bifemsh_'+side,'bifemlh_'+side,'soleus_'+side,'tib_ant_'+side):
                mi=acts.index(act);entry['MN_current_groups'][act]={g:stats(z['currents'][mask,gi,mi])for gi,g in enumerate(groups)}
                entry['MN_current_groups'][act]['actual_force_N']=stats(z['force'][mask,mi])
                entry['MN_current_groups'][act]['control']=stats(z['ctrl'][mask,mi])
                entry['MN_current_groups'][act]['actual_activation']=stats(z['act'][mask,mi])
            row=jointnames.index('knee_angle_'+side);contributions=z['rho'][:,row,:]*active_force
            for term in(0,1,2,3,4):entry['knee_torque_Nm'][z['torque_terms'].tolist()[term]]=stats(z['terms'][mask,term,row])
            ranked=np.argsort(contributions[mask].mean(axis=0));take=np.r_[ranked[:5],ranked[-5:]]
            entry['knee_route_active_torque_Nm']={acts[ai]:stats(contributions[mask,ai])for ai in take}
            entry['scope']='F/E assigned from actual ipsilateral RG contrast; positive native knee torque extends. Muscle active/passive split is the recorder zero-activation probe definition, not a new constitutive fit.'
            report['phase_current_and_torque'][side+'_'+phase]=entry
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--campaign',required=True,type=Path);a=p.parse_args()
    out=a.campaign/'figures'/'ground_freeze_diagnosis';out.mkdir(exist_ok=True)
    reports={name:review(a.campaign,name)for name in('balanced_hunt23_s2','balanced_hunt23_s4')}
    (out/'HIERARCHY_CURRENTS.json').write_text(json.dumps(reports,indent=2),encoding='utf-8')
    for label,report in reports.items():
        print(label)
        for name in('RG_E_r','RG_F_r','PF_HIP-E_r','PF_HIP-F_r','PF_KNEE-E_r','PF_KNEE-F_r','PF_BI_JOINT_RECT_r','PF_BI_JOINT_HAM_r'):
            print(name,report['hierarchy'].get(name))
        print('Timing currents',report['HUNT_AF_to_RG_current'])
        print('Native knee torque by actual RG phase',{phase:row['knee_torque_Nm']for phase,row in report['phase_current_and_torque'].items()})
        for name in('MN_vas_int_r','MN_bifemsh_r','MN_soleus_r','MN_tib_ant_r'):
            print(name,report['selected_voltage_and_current'][name])
if __name__=='__main__':main()
