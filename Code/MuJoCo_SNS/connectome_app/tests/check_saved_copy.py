"""Read-only verification of the native UI's working-copy output."""
import json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from project import normalize
from runtime import build_project

saved=Path(sys.argv[1]);fixture=Path(__file__).resolve().parents[1]/'examples/paired_mn.sns.json'
actual=normalize(json.loads(saved.read_text(encoding='utf-8')))
source=normalize(json.loads(fixture.read_text(encoding='utf-8')))
assert actual['source']==source['source']
assert actual['interface']==source['interface']
assert {n['id'] for n in actual['nodes']}=={n['id'] for n in source['nodes']}
assert {e['id'] for e in actual['edges']}=={e['id'] for e in source['edges']}
changes=[]
for key in ('nodes','edges'):
    original={x['id']:x for x in source[key]}
    for row in actual[key]:
        assert row.get('source_record')==original[row['id']].get('source_record')
        if key=='edges':
            for field in ('fromId','toId','map'):
                assert row.get(field)==original[row['id']].get(field),(row['id'],field)
            if row.get('synapse')!=original[row['id']].get('synapse'):
                changes.append({'id':row['id'],'synapse':row.get('synapse')})
model=build_project(actual)
assert model['steps_executed']==0
for changed in changes:
    assert changed['synapse']['gmax_uS'] in [e['params']['max_conductance'] for e in model['network'].connections]
print('Edited chemical edges:',json.dumps(changes))
print(json.dumps({k:v for k,v in model.items() if k not in ('network','compiled')},indent=2))
