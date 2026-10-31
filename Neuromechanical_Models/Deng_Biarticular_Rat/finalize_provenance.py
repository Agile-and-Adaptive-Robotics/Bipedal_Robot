"""Refresh hashes and explicit local adaptation inventory without touching models."""
from pathlib import Path
import hashlib,json
HERE=Path(__file__).resolve().parent
p=HERE/'SOURCE_PROVENANCE.json';d=json.loads(p.read_text())
d.pop('standalone_alias',None)
d['preserved_originals']=[str(HERE/'Source/Original_Biarticular_net.aproj'),str(HERE/'Biarticular_net_Standalone.asim')]
d['standalone_alias']='Biarticular_net.asim is an instrumented, finite, local-path companion, not a byte-identical original.'
d['local_adaptations']=['Local paths repaired; visual grass texture cleared; finite 5.01-second run; existing-state recording; replay disabled. No neural or muscle parameter edits.']
d['local_files']=[dict(file=p.name,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in [HERE/'Biarticular_net.aproj',HERE/'Biarticular_net.asim',HERE/'Biarticular_net_Standalone.asim',HERE/'Source/Original_Biarticular_net.aproj']]
expected={'Biarticular_net_Standalone.asim':'a7f67185c2c9ce11223f253bb1e0620554dd98c69d49f6d03d6afbbae0e97c23','Original_Biarticular_net.aproj':'126648ee594b0c409b4b0b78e8be6cd3664d8fb5119f43d5524fac535e7d2427'}
for x in d['local_files']:
 if x['file'] in expected:assert x['sha256']==expected[x['file']],x
p.write_text(json.dumps(d,indent=2));print('Both preserved originals match their verified public source hashes.')
