"""Four declared bilateral contact-path ablations; all other pathways unchanged."""
import json
import numpy as np
from scipy.signal import find_peaks
from import_virtual_walker import HERE
from run_virtual_walker import run
from assess_trial import assess
base=HERE/'runs/W2L_VirtualWalker_Newer';rows=[]
cases=[('contact_RG_off',dict(contact_rg_scale=0.)),('contact_PF_off',dict(contact_pf_scale=0.)),('contact_MN_off',dict(contact_mn_scale=0.)),('all_contact_outputs_off',dict(contact_rg_scale=0.,contact_pf_scale=0.,contact_mn_scale=0.))]
for label,params in cases:
 tag='canonical108_ground_'+label+'_v1';folder=base/tag;r=run(newer=True,duration=10.,tag=tag,**params);a=assess(folder);z=np.load(folder/'trajectory.npz');mask=z['t']>=1.;names=z['neuron_names'].tolist();rg={}
 for name in ['L RG ext','L RG flx','R RG ext','R RG flx']:
  v=z['neuro'][mask,names.index(name)];pk,_=find_peaks(v,prominence=1.,distance=round(.2/r['physics_dt']));rg[name]=dict(peaks=len(pk),range_mV=float(np.ptp(v)))
 contacts=[]
 for i,name in enumerate(['left_heel','left_toe','right_heel','right_toe']):
  signal=z['contact_sensor_binary'][mask,i];contacts.append(dict(sensor=name,on_fraction=float(np.mean(signal)),rising_edges=int(np.sum(np.diff(signal)>0))))
 out=dict(label=label,result=r,assessment=a,RG=rg,actual_contact_signals=contacts,scope='Exact16 direct contact edges classified by actual target layer; only selected conductances set to0. Contact sensors, Ia, legacyII,12homonymousII,Ib,c1,V3 and all other weights retained. Absolute homologous L/R gains asserted equal.');(folder/'CONTACT_ABLATION.json').write_text(json.dumps(out,indent=2));rows.append(out);(HERE/'runs/CANONICAL108_CONTACT_ABLATIONS.json').write_text(json.dumps(rows,indent=2));print(json.dumps(dict(label=label,swings=[a['legs'][s]['unloaded_clear_swings'] for s in ('L','R')],alternation=a['alternating_transition_fraction'],support=r['mean_net_vertical_assistance_BW'],failure=r['failure']),indent=2),flush=True)
