"""Read-only source-copy probe: physical contact -> actual adapter output -> neuron current."""
from pathlib import Path
from lxml import etree as E
import copy,json,uuid,subprocess,numpy as np,hashlib
R=Path(__file__).resolve().parents[3];S=R/'Neuromechanical_Models/W2L_VirtualWalker_Newer_20261031/W2L_VirtualWalker_Newer_Standalone.asim';W=R/'Neuromechanical_Models/W2L_VirtualWalker_20261031/build/native_contact_unit_probe_20261031';W.mkdir(exist_ok=True)
tree=E.parse(str(S));root=tree.getroot();root.find('ProjectPath').text=str(W)+'\\';root.find('SimEndTime').text='.51';charts=root.find('DataCharts');base=copy.deepcopy(charts[0]);charts.clear();charts.append(base)
def put(n,k,v):
 z=n.find(k)
 if z is None:z=E.SubElement(n,k)
 z.text=str(v)
def uid(s):return str(uuid.uuid5(uuid.NAMESPACE_URL,'native-contact-unit-probe-20261031/'+s))
for k,v in [('ID',uid('chart')),('Name','Native contact unit probe'),('OutputFilename','Native contact unit probe.txt'),('EndTime','.5'),('CollectInterval','.001')]:put(base,k,v)
cols=base.find('DataColumns');cols.clear();labels=['TimeSlice','Time'];entries=[]
for a in root.findall('.//Adapters/Adapter'):
 if a.findtext('SourceDataType')!='ContactCount':continue
 name=a.findtext('Name');entry={'adapter_id':a.findtext('ID'),'source_id':a.findtext('SourceID'),'target_id':a.findtext('TargetID'),'RobotIOScale':a.findtext('RobotIOScale'),'DelayBufferMode':a.findtext('DelayBufferMode'),'gain':E.tostring(a.find('Gain'),encoding='unicode')};entries.append(entry)
 for suffix,target,kind in [('contact_count',a.findtext('SourceID'),'ContactCount'),('calculated_output',a.findtext('ID'),'CalculatedVal'),('next_output',a.findtext('ID'),'NextVal'),('neuron_adapter_current_A',a.findtext('TargetID'),'AdapterCurrent'),('neuron_vm_V',a.findtext('TargetID'),'MembraneVoltage')]:
  label=name+' '+suffix;c=E.SubElement(cols,'DataColumn')
  for k,v in [('ID',uid(label)),('Type','DataColumn'),('ColumnName',label),('TargetID',target),('DataType',kind)]:put(c,k,v)
  labels.append(label)
model=W/'native_probe.asim';tree.write(str(model),encoding='utf-8',xml_declaration=True);exe=Path(r'D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin\AnimatSimulator.exe');q=subprocess.run([str(exe),str(model)],cwd=exe.parent,capture_output=True,text=True,timeout=20);(W/'run.log').write_text(q.stdout+'\n'+q.stderr);print('exit',q.returncode,q.stdout[-1000:])
report={'source':str(S),'source_sha256':hashlib.sha256(S.read_bytes()).hexdigest(),'source_coefficients_unchanged':True,'exit':q.returncode,'adapter_xml':entries,'labels':labels}
file=W/'Native contact unit probe.txt'
if file.exists() and file.stat().st_size:
 a=np.loadtxt(file,skiprows=1);a=a[(a[:,1]<=.5)&np.any(a[:,2:]!=0,axis=1)];np.savez_compressed(W/'probe.npz',data=a);report['samples']=len(a);report['cases']=[]
 for j,entry in enumerate(entries):
  offset=2+j*5;counts=a[:,offset];case={'adapter':entry['adapter_id'],'contact_values':np.unique(counts).tolist()}
  for c in (0.,1.):
   sel=np.flatnonzero((counts==c)&(a[:,1]>=.005))
   if len(sel):
    rows=a[sel,offset:offset+5];case[str(c)]={'first_time':float(a[sel[0],1]),'first_row':rows[0].tolist(),'calc_minmax':[float(rows[:,1].min()),float(rows[:,1].max())],'current_minmax_A':[float(rows[:,3].min()),float(rows[:,3].max())]}
  report['cases'].append(case)
(W/'PROBE_REPORT.json').write_text(json.dumps(report,indent=2));print(json.dumps(report.get('cases',{}),indent=2))
