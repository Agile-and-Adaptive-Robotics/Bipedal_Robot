"""Local path/visual repairs plus a bounded, instrumented source-model run."""
from pathlib import Path
import xml.etree.ElementTree as E
import json,hashlib,uuid,copy,re
HERE=Path(__file__).resolve().parent
OUT=HERE/'runs/native_source_smoke';OUT.mkdir(parents=True,exist_ok=True)
original=(HERE/'Source/Original_Biarticular_net.aproj').read_bytes().decode('utf-8')
text=re.sub(r'<ProjectPath>.*?</ProjectPath>','<ProjectPath>'+str(HERE).replace('\\','/')+'/</ProjectPath>',original,count=1)
text=re.sub(r'<Texture>.*?</Texture>','<Texture/>',text)
text=text.replace('<SetSimEnd>False</SetSimEnd>','<SetSimEnd>True</SetSimEnd>',1)
text=re.sub(r'<SimEndTime\b[^>]*/>','<SimEndTime Value="5.01" Scale="None" Actual="5.01"/>',text,count=1)
(HERE/'Biarticular_net.aproj').write_bytes(text.encode('utf-8'))
ap=E.fromstring(text)
pages=0
for d in ap.iter('DiagramXml'):
 if d.text and d.text.strip():E.fromstring(d.text);pages+=1
root=E.parse(HERE/'Biarticular_net_Standalone.asim').getroot()
for name,value in [('ProjectPath',str(OUT).replace('\\','/')+'/'),('SetSimEnd','True'),('SimEndTime','5.01'),('EnableSimRecording','False')]:root.find(name).text=value
for e in root.iter('Texture'):e.text=''
charts=root.find('DataCharts');chart=copy.deepcopy(charts.find('DataChart'))
chart.find('ID').text=str(uuid.uuid5(uuid.NAMESPACE_URL,'deng-biarticular-source-baseline'))
chart.find('Name').text='Source baseline RG PF MN joints root'
chart.find('OutputFilename').text='source_baseline.txt'
for name,value in [('StartTime','0'),('EndTime','5'),('CollectInterval','.001')]:
 e=chart.find(name)
 if e is not None:e.text=value
columns=chart.find('DataColumns');template=copy.deepcopy(columns.find('DataColumn'));columns.clear();inventory=[]
def add(name,ident,kind):
 e=copy.deepcopy(template)
 for key,value in [('ID',str(uuid.uuid5(uuid.NAMESPACE_URL,'deng-baseline-'+name))),('ColumnName',name),('TargetID',ident),('DataType',kind)]:e.find(key).text=value
 columns.append(e);inventory.append(dict(name=name,target=ident,type=kind))
for n in root.iter('Neuron'):
 name=n.findtext('Name') or '';clean=name.replace('_',' ')
 if (' RG ' in clean or ' PF ' in clean) and ' IN ' not in clean or ' MN ' in clean:
  add(name,n.findtext('ID'),'MembraneVoltage')
for j in root.iter('Joint'):
 name=j.findtext('Name') or ''
 if any(k in name.lower() for k in ('hip','knee','ankle')):add(name+' angle',j.findtext('ID'),'JointRotation')
rb=root.find('.//Organism/RigidBody');add('Root height',rb.findtext('ID'),'WorldPositionY')
charts.append(chart)
E.indent(root);E.ElementTree(root).write(OUT/'Biarticular_source_smoke.asim',encoding='utf-8',xml_declaration=True)
# Keep the project-named runtime companion local and usable as well.
runtime=copy.deepcopy(root);runtime.find('ProjectPath').text=str(HERE).replace('\\','/')+'/'
E.ElementTree(runtime).write(HERE/'Biarticular_net.asim',encoding='utf-8',xml_declaration=True)
result=dict(project=str(HERE/'Biarticular_net.aproj'),runtime=str(OUT/'Biarticular_source_smoke.asim'),
 original_source_preserved=True,neural_or_muscle_parameters_changed=False,diagram_pages_validated=pages,
 adaptations=['Local ProjectPath repaired','Unavailable grass texture cleared (visual only)','Finite 5.01-second run; source charts end at 5 seconds','Added recording of existing RG/PF/MN/joint/root signals','Replay recording disabled; chart recording retained'],
 columns=inventory,source_runtime_sha256=hashlib.sha256((HERE/'Biarticular_net_Standalone.asim').read_bytes()).hexdigest())
(OUT/'PREPARATION.json').write_text(json.dumps(result,indent=2))
print(json.dumps(dict(project=result['project'],runtime=result['runtime'],columns=len(inventory),diagram_pages=pages),indent=2))
