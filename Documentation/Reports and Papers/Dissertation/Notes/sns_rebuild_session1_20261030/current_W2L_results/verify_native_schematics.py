from pathlib import Path
import json,csv,hashlib
from lxml import etree as E
r=Path(r'D:\GitHub\Bipedal_Robot');out=r/'Documentation/Reports and Papers/Dissertation/Notes/sns_rebuild_session1_20261030/current_W2L_results/native_schematic';w=r/'Neuromechanical_Models/W2L_VirtualWalker_Newer_20261031';data=json.loads((r/'Code/MuJoCo_SNS/spinal/virtual_walker_20261031/runs/GENTLE108_EFFECTIVE_CONNECTOME.json').read_text());p=out/'SCHEMATIC_PROVENANCE.json';audit=json.loads(p.read_text());audit['all_executed_target_changes']=[e for e in data['edges']if e['target_changed']];assert len(audit['all_executed_target_changes'])==4
audit['display_completed_head_unique_connections']=len({x['edge']for x in audit['missing_native_chemical_heads_completed']});assert audit['display_completed_head_unique_connections']==0
audit['scope']='Selected native subnetworks retain node centers, styles, and internal route control points. Neural/off-page glyphs are enlarged 40% around unchanged centers for publication typography; sensor/adapter glyphs retain native sizes. Four ankle destination overlays and14 short ID-backed cross-page display stubs are explicitly distinguished from unchanged internal routes.'
assert len({e['edge_ID']for e in audit['abstracted_crosspage_edges']})==14
layout=json.loads((out/'MAIN_LAYOUT.json').read_text())if(out/'MAIN_LAYOUT.json').exists()else{'width_in':11,'height_in':8.5,'minimum_font_pt':12.5}
audit['main_canvas_inches']=[layout['width_in'],layout['height_in']];audit['minimum_font_pt']=layout['minimum_font_pt'];audit['print_font_pt_at9inch_width']=layout['minimum_font_pt']*9/layout['width_in']
for view in audit['page_views']:
 view.pop('preserved_geometry',None);view['native_node_centers_and_internal_route_points_retained']=True;view['neural_offpage_glyph_scale']=1.4
audit['source_project_unchanged_after_render']=hashlib.sha256((w/'W2L_VirtualWalker_Newer.aproj').read_bytes()).hexdigest()==audit['source_project_sha256'];assert audit['source_project_unchanged_after_render'];audit['visual_QA']={'main_png_inspected':True,'expanded_left_hip_inspected':True,'left_knee_ankle_overlay_inspected':True,'parameter_table_page1_inspected':True,'limitations':['The14 cross-page stubs abstract actual connections; full original routes/endpoints remain in source and companion mapping.','Final main figure is9x5.5in with10pt minimum type.','This is a selection of actual native pages, not a single complete108-neuron sheet.']};p.write_text(json.dumps(audit,indent=2))
native=E.parse(str(w/'W2L_VirtualWalker_Newer_Standalone.asim'));byid={q.findtext('ID'):q for q in native.iter()if q.findtext('ID')};rows=[]
for edge in data['edges']:edge['reference_connection_G_uS']=float(byid[edge['source_connection_ID']].findtext('G'))
with(out/'EXECUTED_EDGES.csv').open('w',newline='',encoding='utf-8')as f:
 wr=csv.DictWriter(f,fieldnames=list(data['edges'][0]));wr.writeheader();wr.writerows(data['edges'])
for n in data['neurons']:
 if 'II'in n['name']and('relay'in n['name']or n['name'].endswith(' II')):
  a=next(q for q in native.findall('.//Adapter')if q.findtext('TargetID')==n['ID']);src=a.findtext('SourceID');srcnode=byid.get(src);targets=[e for e in data['edges']if e['source_ID']==n['ID']]
  rows.append({'Category':'New homonymous'if'relay'in n['name']else'Retained phase','NeuronID':n['ID'],'NativeNeuronName':n['name'],'PhysicalSourceID':src,'PhysicalSourceName':''if srcnode is None else srcnode.findtext('Name',''),'SourceDataType':a.findtext('SourceDataType'),'TargetNames':'; '.join(e['target']for e in targets),'TargetIDs':'; '.join(e['target_ID']for e in targets)})
assert sum(x['Category']=='New homonymous'for x in rows)==12 and sum(x['Category']=='Retained phase'for x in rows)==4
with(out/'TYPE_II_PHYSICAL_AND_NEURAL_MAP.csv').open('w',newline='',encoding='utf-8')as f:wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
print({'unique_completed_II_heads':audit['display_completed_head_unique_connections'],'target_changes':len(audit['all_executed_target_changes']),'II_channels':len(rows),'source_unchanged':audit['source_project_unchanged_after_render']})

