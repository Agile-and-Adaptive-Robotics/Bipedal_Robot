"""Assemble and statically verify revised drawing outputs. No circuit execution."""
import json
import xml.etree.ElementTree as E
import pymupdf
from audit_saved_graphs import OUT,ROOT,digest

b=OUT/'drawings'
w2l=['W2L108_L_central_feedback','W2L108_L_central_feedback_key','W2L108_R_central_feedback','W2L108_R_central_feedback_key']
panels=w2l+['w2lvar_central_feedback_RG','w2lvar_central_feedback_PF','w2lvar_muscle_bundle_members','syn6_central_feedback_RG','syn6_central_feedback_PF','s3k_central_feedback_RG_PF']
for name,items in [('W2L108_RG_PF_Feedback',w2l),('RG_PF_Feedback_revised',panels)]:
    d=pymupdf.open()
    for slug in items:
        with pymupdf.open(b/(slug+'.pdf')) as p:d.insert_pdf(p)
    d.save(b/(name+'.pdf'));d.close()
checks=[]
for slug in panels:
    E.parse(b/(slug+'.svg'))
    with pymupdf.open(b/(slug+'.pdf')) as p:
        assert len(p)==1 and p[0].get_text().strip()
    checks.append({'drawing':slug,'svg_xml_valid':True,'pdf_single_page_with_text':True,'svg_sha256':digest(b/(slug+'.svg'))})
represented=0
for family,scope in [('w2lvar','RG'),('w2lvar','PF'),('syn6','RG'),('syn6','PF'),('s3k','RG_PF')]:
    key=json.loads((b/(family+'_central_feedback_'+scope+'_key.json')).read_text())
    original=json.loads((ROOT/('Neuromechanical_Models/Mujoco_SNS_models/Gait2392_standing_families/'+family+'_graph_ledger.json')).read_text())
    edge_map={e['id']:e for e in original['chemical_edges']}
    for bundle in key['bundles_and_relays']:
        for edge in bundle['full_edge_fields']:
            assert edge==edge_map[edge['id']]
            represented+=1
        if 'connection_count' in bundle:assert bundle['connection_count']==len(bundle['actual_edge_ids'])
record={'pages':len(panels),'checks':checks,'represented_bundle_or_relay_edges_verified':represented,'triangle_convention':'Point on incoming line; broad base on receiving object','feedback_categories':['Ia','II','Ib','heel','toe'],'W2L108_INs':'Both RG/PF interneurons explicitly drawn in each module; source V3→contra RG ext IN included','missing_routes':'Shown as absent; not invented','bundles':'Exact source channel counts on lines; named members in w2lvar roster and per-page JSON keys','circuits_executed':0}
(b/'FEEDBACK_REVISION_VERIFICATION.json').write_text(json.dumps(record,indent=2))
print(json.dumps({'pages':len(panels),'verified_edges':represented,'circuits_run':0}))
