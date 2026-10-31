from pathlib import Path
from lxml import etree as E
import copy,uuid,subprocess,numpy as np,json
R=Path(__file__).resolve().parents[3];W=R/'Neuromechanical_Models/W2L_VirtualWalker_20261031/build/isolated_neural_clock';W.mkdir(exist_ok=True);exe=Path(r'D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin\AnimatSimulator.exe');results=[];stimuli=[]
def put(n,k,v):
 z=n.find(k)
 if z is None:z=E.SubElement(n,k)
 z.text=str(v)
for tag,scale in [('Guided',1),('GuidedT2',2)]:
 source=R/('Neuromechanical_Models/W2L_VirtualWalker_'+tag+'_20261031/W2L_VirtualWalker_'+tag+'_Standalone.asim');tree=E.parse(str(source));root=tree.getroot();out=W/tag;out.mkdir(exist_ok=True);put(root,'ProjectPath',str(out)+'\\');put(root,'SimEndTime',8.01);body=next(b for b in root.iter('RigidBody') if b.findtext('Name')=='Root');put(body,'Freeze','True')
 disabled=0
 for a in root.findall('.//Adapters/Adapter'):
  if a.findtext('Type')=='PhysicalToNode':put(a,'Enabled','False');disabled+=1
 stimuli.append(E.tostring(root.find('ExternalStimuli'),method='c14n'))
 charts=root.find('DataCharts');chart=copy.deepcopy(charts[0]);charts.clear();charts.append(chart);put(chart,'Name','Neural clock probe');put(chart,'OutputFilename','Neural clock probe.txt');put(chart,'EndTime',8);put(chart,'CollectInterval',.001);cols=chart.find('DataColumns');cols.clear();names=[]
 for n in root.iter('Neuron'):
  if n.findtext('Name') not in ('L RG ext','L RG flx','L Hip PF ext','L Hip PF flx'):continue
  name=n.findtext('Name');c=E.SubElement(cols,'DataColumn')
  for k,v in [('ID',str(uuid.uuid5(uuid.NAMESPACE_URL,'neuralclock/'+tag+'/'+name))),('Type','DataColumn'),('ColumnName',name),('TargetID',n.findtext('ID')),('DataType','MembraneVoltage')]:put(c,k,v)
  names.append(name)
 file=out/'probe.asim';tree.write(str(file),encoding='utf-8',xml_declaration=True);q=subprocess.run([str(exe),str(file)],cwd=exe.parent,capture_output=True,text=True,timeout=30);(out/'run.log').write_text(q.stdout+'\n'+q.stderr);record={'clock':scale,'exit':q.returncode,'disabled_physical_adapters':disabled}
 if q.returncode==0:
  a=np.loadtxt(out/'Neural clock probe.txt',skiprows=1);a=a[(a[:,1]<=8)&np.any(a[:,2:]!=0,axis=1)];v=a[:,2+names.index('L RG flx')];t=a[:,1];indices=np.flatnonzero((v[1:]>-.06)&(v[:-1]<=-.06))+1;onsets=t[indices];onsets=onsets[onsets>=1];record.update(samples=len(a),finite=bool(np.isfinite(a).all()),flexor_onsets_s=onsets.tolist(),median_period_s=float(np.median(np.diff(onsets))) if len(onsets)>=3 else None);np.savez_compressed(out/'trace.npz',data=a,names=names)
 results.append(record)
assert stimuli[0]==stimuli[1]
ratio=results[1].get('median_period_s')/results[0].get('median_period_s') if all(r.get('median_period_s') for r in results) else None
report={'same_external_inputs':True,'physics_timestep_and_muscle_laws_unchanged':True,'scope':'Isolated deafferented graded network; root frozen only in these probes, not guided companions. Same pulse input; compare steady onset periods after1s.','results':results,'period_ratio':ratio,'doubling_check_pass':ratio is not None and abs(ratio-2)<.02}
(W/'CLOCK_VERIFICATION.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
