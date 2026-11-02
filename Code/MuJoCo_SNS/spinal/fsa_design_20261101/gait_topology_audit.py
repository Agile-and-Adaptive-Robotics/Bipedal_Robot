"""Read-only topology construction/audit: no neural or physics time steps."""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import argparse,hashlib,json,sys
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
HERE=Path(__file__).resolve().parent;S=HERE.parent;REPO=S.parents[2]
sys.path[:0]=[str(S),str(S/'w2l_cpg')]
EXPECTED={'s3k':'0e232bd04001e342a774d3de08ba5d0ed378907673bc2c5ccfa81daae22ae743','w2lvar':'49b946c627222e030aa7bc741d7de20faeb804646d791b51548e6d3d704d4e2f','syn6':'0843e32e4e5cf16de974720620c0101e05bf1424018c7815d21b8b1182d6e5b0'}


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def summarize(names,edges):
    kinds=('IaIN_','RC_','IIX_','IIIN_','IBIN_','IBEXC_')
    counts={k:sum(n.startswith(k) for n in names) for k in kinds}
    direct=[e for e in edges if e['source'].startswith(('Ia_','II_','Ib_')) and e['target'].startswith('MN_')]
    central=[e for e in edges if e['source'].startswith(('Ia_','II_','Ib_')) and e['target'].startswith(('RG_','PF_','InE_','InF_'))]
    direct_counts={kind:{'excitatory':sum(e['source'].startswith(kind+'_') and e['reversal_mV']>0 for e in direct),'inhibitory':sum(e['source'].startswith(kind+'_') and e['reversal_mV']<0 for e in direct)} for kind in ('Ia','II','Ib')}
    central_counts={kind:{'edges':sum(e['source'].startswith(kind+'_') for e in central),'unique_g_uS':sorted({e['gmax_uS'] for e in central if e['source'].startswith(kind+'_')})} for kind in ('Ia','II','Ib')}
    heel=[e for e in edges if e['source'].startswith('HEEL_')];two=[]
    for a in heel:
        if not a['target'].startswith(('InE_','InF_')):continue
        for b in edges:
            if b['source']==a['target'] and b['target'].startswith('RG_'):
                two.append({'source':a['source'],'via':a['target'],'target':b['target'],'first_g_uS':a['gmax_uS'],'first_reversal_mV':a['reversal_mV'],'second_g_uS':b['gmax_uS'],'second_reversal_mV':b['reversal_mV'],'two_hop_sign_at0to5mV':int(np.sign(a['reversal_mV'])*np.sign(b['reversal_mV'])),'not_a_product_gain':'Conductance shunting,thresholds,and phase state matter'})
    return {'nodes':len(names),'nonzero_pairs':len(edges),'local_interneuron_counts':counts,'direct_afferent_to_MN':direct_counts,'direct_central':central_counts,'heel_edges':heel,'heel_two_hop':two}


def raw(n):
    c=n.compiled;names=[p['name'] for p in n.net.populations];edges=[{'source':names[j],'target':names[i],'gmax_uS':float(c.g_max_non[i,j]),'reversal_mV':float(c.del_e[i,j])} for i,j in zip(*np.nonzero(c.g_max_non))]
    r=summarize(names,edges);r['compiled_chemical_sha256']={k:hashlib.sha256(np.ascontiguousarray(getattr(c,k)).tobytes()).hexdigest() for k in ('g_max_non','del_e','e_lo','e_hi')};return r


def audit(family):
    graph_path=REPO/'Neuromechanical_Models/Mujoco_SNS_models/Gait2392_standing_families'/(family+'_graph_ledger.json');graph=json.loads(graph_path.read_text())
    if sha(graph_path)!=EXPECTED[family]:raise ValueError('Selected source ledger differs from verified EB475WS4 ledger')
    os.environ['AARL_NET']=family
    import build_network as B,params as P
    inputs=[e.attrib['name'] for e in ET.parse(graph['actuation']['model']).getroot().find('actuator')]
    default=B.build(inputs);default_G=P.G.copy();first=raw(default)
    P.G['full_rules']=1.;full=B.build(inputs);effective=float(P.G['full_rules']);second=raw(full)
    selected=summarize([n['name'] for n in graph['nodes']],graph['chemical_edges'])
    return {'role':'Gait neural pathway composition reviewer','family':family,'selected_source_ledger':str(graph_path),'selected_source_sha256':sha(graph_path),'EB475WS4_exact_hash_match':True,
            'raw_default':first,'requested_full_rules1':second,'requested_full_rules_value':1.,'effective_full_rules_after_variant_build':effective,
            'default_vs_requested_full_rules1_chemical_matrices_identical':first['compiled_chemical_sha256']==second['compiled_chemical_sha256'],
            'default_G_relevant':{k:default_G[k] for k in ('full_rules','ia_in','renshaw','heel_rge','toe_rge','ib_rge','heel_in_f_exc','ia_f_central','ii_f_central','ii_e_central','ib_e_central')},
            'selected_standing':selected,'saved_masks':graph['settings'],'deafferentation_scope':'Runner use_aff=False/AARL_DEAFF modifies encoded currents,not these graph matrices. No-afferents alone is not asserted to mute all separately encoded contact/HUNT_AF paths.',
            'source_sha256':{str(p):sha(p) for p in (S/'build_network.py',S/'build_network_w2lvar.py',S/'build_network_syn6.py',S/'params.py',S/'runner.py')},
            'neural_steps':0,'physics_steps':0,'training':False,'wiring_files_edited':False,'global_release':False}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--family',choices=tuple(EXPECTED),required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();r=audit(a.family);a.out.mkdir(parents=True,exist_ok=True);f=a.out/(a.family+'_TOPOLOGY.json');f.write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps({'family':a.family,'raw':r['raw_default']['local_interneuron_counts'],'full_rules1':r['requested_full_rules1']['local_interneuron_counts'],'effective_full_rules':r['effective_full_rules_after_variant_build'],'neural_steps':0,'output':str(f)}))
