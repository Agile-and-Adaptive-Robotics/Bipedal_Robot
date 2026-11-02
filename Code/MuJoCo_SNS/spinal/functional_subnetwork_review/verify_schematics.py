"""Check drawing data against saved sources; never import a circuit."""
import json
from pathlib import Path
import xml.etree.ElementTree as E
import pymupdf
from audit_saved_graphs import ROOT,OUT

DIR=OUT/'drawings/v4'
BASE=ROOT/'Neuromechanical_Models/Mujoco_SNS_models'
def native(path):
    r=E.parse(path).getroot();nodes={n.findtext('ID'):n for n in r.iter('Neuron')};types={t.findtext('ID'):t for t in r.iter('SynapseType')};rows=[]
    for e in r.iter('Connexion'):
        a,b=e.findtext('SourceID'),e.findtext('TargetID');t=types[e.findtext('SynapseTypeID')]
        rows.append({'id':e.findtext('ID'),'source':nodes[a].findtext('Name'),'target':nodes[b].findtext('Name'),'gmax_uS':float(t.findtext('SynAmp')),'reversal_mV':float(t.findtext('Equil')),'threshold_mV':float(t.findtext('ThreshV')),'saturation_mV':float(t.findtext('SaturateV'))})
    return {e['id']:e for e in rows}
indexes={f:{e['id']:e for e in json.loads((BASE/('Gait2392_standing_families/'+f+'_graph_ledger.json')).read_text())['chemical_edges']} for f in ('s3k','w2lvar','syn6')}
indexes['W2L108']=native(ROOT/'tmp/functional_subnetwork_remote/W2L108_verified/editable/W2L_VirtualWalker_Newer_Standalone.asim')
indexes['Rat134']=native(ROOT/'tmp/functional_subnetwork_remote/named_conditions/rat134/Rat134_ReceptorGeometry_Standalone.asim')
indexes['W2L704']={e['id']:e for e in json.loads((BASE/'W2L_704_80_baseline/W2L_704_80_graph_ledger.json').read_text())['chemical_edges']}
rf=json.loads((BASE/'Hunt_fourcell_session1_left_RF_HAM_connectome.json').read_text())['edges']
indexes['RF_HAM']={f'RFHAM:left:{i}':{'source':e['from'],'target':e['to'],'gmax_uS':e['gain']} for i,e in enumerate(rf)}
results=[]
for p in sorted(DIR.glob('*_connections.json')):
    d=json.loads(p.read_text());slug=p.name.removesuffix('_connections.json')
    family='syn6' if slug.startswith('Gait_') else next(f for f in indexes if slug.startswith(f))
    count=0
    for edge in d['source_connections']:
        source=indexes[family][edge['id']]
        for key in ('source','target','gmax_uS','reversal_mV','threshold_mV','saturation_mV'):
            if key in source:assert edge[key]==source[key],(slug,edge['id'],key)
        count+=1
    geometry=json.loads((DIR/(slug+'_geometry_check.json')).read_text())
    assert not geometry['dangling_endpoints'] and not geometry['wires_through_node_bodies'],slug
    E.parse(DIR/(slug+'.svg'))
    with pymupdf.open(DIR/(slug+'.pdf')) as pdf:assert len(pdf)==1 and pdf[0].get_text().strip()
    results.append({'drawing':slug,'source_family':family,'source_edge_records_checked':count,'source_endpoint_and_parameter_parity':True,'geometry_pass':True,'svg_xml_pass':True,'pdf_text_pass':True})
with pymupdf.open(DIR/'RG_PF_Schematics.pdf') as combined:pages=len(combined)
report={'drawings':results,'pdf_pages':pages,'represented_edge_records_checked':sum(r['source_edge_records_checked'] for r in results),'circuits_executed':0,'scope':'Checks of displayed saved-source connections and drawing geometry only; not circuit behavior or completeness of every source model.'}
(DIR/'VERIFICATION.json').write_text(json.dumps(report,indent=2))
print(json.dumps({'drawings':len(results),'pages':pages,'edge_records_checked':report['represented_edge_records_checked'],'circuits_executed':0}))
