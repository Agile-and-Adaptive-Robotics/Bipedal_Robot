from source_inputs import historical_source
"""Read-only original/native/newer GUI RG differences; no fitting."""
import json,xml.etree.ElementTree as E
from import_virtual_walker import HERE,REPO
from native_kernel import NativeController
original=REPO/'Neuromechanical_Models/W2L_VirtualWalker_20261031/W2L_VirtualWalker_Standalone.asim';newer=historical_source('W2L_VirtualWalker_Newer_GUI1nA');bank=json.loads((HERE/'runs/W2L_VirtualWalker/IMPORT_AUDIT.json').read_text())['muscles'];a=NativeController(original,bank,neural_time_scale=2.);b=NativeController(newer,bank,symmetric_rg_ext_tonic=0.,neural_time_scale=2.);fields=['rest_mV','membrane_tau_ms','size','gmaxCa','Ca_m_midpoint','Ca_m_slope','Ca_m_tau_ms','Ca_h_midpoint','Ca_h_slope','Ca_h_tau_ms','Ca_reversal_mV'];diff=[];edges=[]
for name in ['L RG ext','L RG flx','L RG ext IN','L RG flx IN']:
 i=a.names.index(name);j=b.names.index(name);node=dict(name=name,parameters={k:dict(original=float(x),newer=float(y)) for k,x,y in zip(fields,a.params[i],b.params[j]) if x!=y},initial_voltage=dict(original=float(a.v[i]),newer=float(b.v[j])),initial_m=dict(original=float(a.m[i]),newer=float(b.m[j])),initial_h=dict(original=float(a.h[i]),newer=float(b.h[j])));diff.append(node)
 old={(a.names[int(e[0])],a.names[int(e[1])]):e[2:].tolist() for e in a.edges if a.names[int(e[1])]==name};new={(b.names[int(e[0])],b.names[int(e[1])]):e[2:].tolist() for e in b.edges if b.names[int(e[1])]==name}
 for pair in sorted(set(old)|set(new)):
  if old.get(pair)!=new.get(pair):edges.append(dict(source=pair[0],target=pair[1],original=old.get(pair),newer=new.get(pair),fields=['gmax','threshold_mV','saturation_mV','reversal_mV']))
out=dict(original=str(original),newer=str(newer),parameters=diff,incoming_edge_differences=edges,source_startup_original=[dict(target=a.names[i],start=s,end=e,nA=v) for i,s,e,v in a.stimuli],source_startup_newer=[dict(target=b.names[i],start=s,end=e,nA=v) for i,s,e,v in b.stimuli]);(HERE/'runs/newer_rg_isolation_v1/ORIGINAL_NEWER_RG_DIFF.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
