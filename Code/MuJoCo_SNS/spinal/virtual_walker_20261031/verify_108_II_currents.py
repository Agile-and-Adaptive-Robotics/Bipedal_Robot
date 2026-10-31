"""Actual sampled homonymous II→MN currents and physical input binding checks."""
import json
import numpy as np
from import_virtual_walker import HERE,REPO
from native_kernel import NativeController
base=HERE/'runs/W2L_VirtualWalker_Newer';run=base/'canonical108_air_baseline_v2';z=np.load(run/'trajectory.npz');meta=json.loads((base/'IMPORT_AUDIT.json').read_text());ledger=json.loads((HERE/'runs/canonical_12II/BUILD_LEDGER.json').read_text());c=NativeController(meta['source_snapshot'],meta['muscles']);names=z['neuron_names'].tolist();rows=[];mask=z['t']>=1.
for mapping in ledger['mapping']:
 pre=c.ids[mapping['relay_id']];post=c.ids[mapping['MN_id']];edge=next(e for e in c.edges if int(e[0])==pre and int(e[1])==post);input=next(q for q in c.inputs if q[0]==pre);assert input[1]>=0 and input[2]=='II' and meta['muscles'][input[1]]['id']==mapping['sensor_id'];g=edge[2]*np.clip((z['neuro'][mask,pre]-edge[3])/(edge[4]-edge[3]),0,1);current=g*(edge[5]-z['neuro'][mask,post]);assert np.max(current)>0.;rows.append(dict(relay=mapping['relay_name'],sensor=mapping['sensor'],actuator=mapping['actuator'],MN=mapping['MN_name'],physical_bank_index=input[1],source_data_type=input[2],current_mean_nA=float(np.mean(current)),current_max_nA=float(np.max(current)),relay_voltage_range_mV=[float(np.min(z['neuro'][mask,pre])),float(np.max(z['neuro'][mask,pre]))],gmax_uS=float(edge[2])))
out=dict(actual_neurons=len(names),channels_verified=len(rows),source_sha256=json.loads((run/'RESULT.json').read_text())['source_sha256'],scope='Sampled postsynaptic driving currents g(Vpre)*(Erev−VMN), Size1 native-equivalent nA, not source-type names alone.',rows=rows);(run/'HOMONYMOUS_II_CURRENT_VERIFICATION.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
