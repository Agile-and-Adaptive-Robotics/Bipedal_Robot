"""Read actual CAD262 electrical command/reflex topology; no simulations."""
import json,numpy as np
import non_gait_contracts as N

if __name__=='__main__':
    net=N.cad_network();c=net.compiled;names=[p['name'] for p in net.net.populations]
    post,pre=np.nonzero(c.g_max_non)
    rows=[dict(source=names[a],target=names[b],gmax=float(c.g_max_non[b,a]),Erev=float(c.del_e[b,a])) for b,a in zip(post,pre)]
    prefixes=['RC_','IaIN_','IIIN_','IIX_','IBIN_','IBEXC_']
    result=dict(machine='easteregg2',source_neurons=len(names),population_counts={p:sum(name.startswith(p) for name in names) for p in prefixes},active_outgoing_counts={p:sum(row['source'].startswith(p) for row in rows) for p in prefixes},active_POSTURE_to_MN=[row for row in rows if row['source']=='POSTURE' and row['target'].startswith('MN_')],input_ports=[name for name in net.inputs if name.startswith('POST_')],source_files=[N.evidence(N.ROBOT/'run_neural_robot.py'),N.evidence(N.SPINAL/'build_network.py')],physics_steps=0,training_steps=0)
    path=N.RAW_RECORDS/'CAD_DESCENDING_TOPOLOGY_v1.json';assert not path.exists();path.write_text(json.dumps(result,indent=2));print(json.dumps(dict(population_counts=result['population_counts'],outgoing_counts=result['active_outgoing_counts'],posture_MN_edges=len(result['active_POSTURE_to_MN']),direct_POST_inputs=len(result['input_ports']))))
