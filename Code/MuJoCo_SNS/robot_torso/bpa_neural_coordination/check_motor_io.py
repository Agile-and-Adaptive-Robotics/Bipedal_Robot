"""Verify preallocated motor mapping against all recorded actual neural states."""
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np

ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa')
MODEL=Path(r'D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\mjc\gait2392_simbody\gait2392_simbody_cvt3.xml')

def main():
 actors=[e.attrib['name'] for e in ET.parse(MODEL).getroot().find('actuator')]
 cases=[]
 for family,relative in [('syn6','pressure_standing/syn6_nominal_fresh_v1'),('s3k','standing/s3k_baseline_v2'),('w2lvar','standing/w2lvar_baseline_v2')]:
  source=ROOT/relative;z=np.load(source.with_suffix('.npz'));a=np.load(source.with_suffix('.apparatus.npz'));names=np.load(ROOT/('family_parity/'+family+'_parity.npz'))['neuron_names'].tolist()
  assert len(actors)==z['ctrl'].shape[1]
  ids=np.array([names.index('MN_'+n) for n in actors]);v=a['neuro'];offset=len(z['ctrl'])-len(v)
  expected=z['ctrl'][offset:];observed=np.empty_like(expected);buffer=np.empty(len(actors))
  for i,state in enumerate(v):
   np.take(state,ids,out=buffer);np.divide(buffer,5.,out=buffer);np.clip(buffer,0.,1.,out=buffer);observed[i]=buffer
  error=float(np.max(abs(observed-expected)));assert np.array_equal(observed,expected)
  cases.append(dict(family=family,physical_frames=len(v),actuators=len(actors),max_control_error=error,bit_exact=True))
 report=dict(machine='easteregg2',passed=True,cases=cases,scope='Preallocated vector motor mapping versus actual recorded scalar mapping, all 92 actuators in three families. E_HI=5 mV, no phase-machine override.')
 (ROOT/'motor_io_parity.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
