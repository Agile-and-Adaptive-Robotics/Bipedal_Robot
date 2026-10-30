"""Register requested Hunt drawings in the original editor, preserving its library."""
import json,re,hashlib,os
from pathlib import Path
H=Path(__file__).resolve().parent;S=H.parent;repo=S.parents[2];dest=repo/'Neuromechanical_Models/Mujoco_SNS_models'
path=S/'connectome_block_editor.html';text=path.read_text(encoding='utf-8');old=text
def embedded_blob(value):
    start=re.search(r'window\.TPL_DATA_EMBED\s*=\s*',value).end();data,end=json.JSONDecoder().raw_decode(value[start:]);return value[start:start+end]
blob=embedded_blob(text);before=hashlib.sha256(blob.encode()).hexdigest()
canonical=S/'ben_connectome_canonical.json';canonical_hash=hashlib.sha256(canonical.read_bytes()).hexdigest()
root=H/'runs/swing_release_20261030_v1/editor_registration';root.mkdir(exist_ok=True)
backup=root/'connectome_block_editor_before_registration.html'
if not backup.exists():backup.write_text(old,encoding='utf-8')
specs={}
for suffix,title in [('subsystems','Hunt Session1 - RF HAM modules'),('full','Hunt Session1 - exact full graph'),('right_RF_HAM','Hunt Session1 - right RF HAM'),('left_RF_HAM','Hunt Session1 - left RF HAM')]:
    spec=json.loads((dest/('Hunt_fourcell_session1_'+suffix+'_connectome.json')).read_text());spec['_note']='Actual trained four-cell candidate; standing/gait-expression conductances at one pose. Ground stepping failed. Boundary bundles in the modules view are display references; exact graph is authoritative.'
    specs['hunt_session1_'+suffix]=dict(name=title,spec=spec)
begin='// HUNT_SESSION1_TEMPLATES_BEGIN';end='// HUNT_SESSION1_TEMPLATES_END'
insertion=begin+'\nconst HUNT_SESSION1_SPECS = '+json.dumps(specs).replace('</','<\\/')+';\nObject.entries(HUNT_SESSION1_SPECS).forEach(([key, value]) => {\n  TEMPLATES[key] = {grp:"walker", name:value.name, build:() => JSON.parse(JSON.stringify(value.spec))};\n});\n'+end+'\n'
if begin in text:text=re.sub(re.escape(begin)+r'[\s\S]*?'+re.escape(end)+r'\n?',lambda m:insertion,text,count=1)
else:
    anchor='(function renderTplOptions() {';assert text.count(anchor)==1;text=text.replace(anchor,insertion+'\n'+anchor,1)
# Retain the outer model when exporting/persisting from inside a SUB.
found=re.search(r'function modelSpec\(\) \{[\s\S]*?\n\}',text);assert found
fixed='''function modelSpec() {
  let ns = nodes, es = edges, gn = grpNames, hid = hiddenG;
  for (let i = ctxStack.length - 1; i >= 0; i--) {
    const f = ctxStack[i];
    f.hostNode.sub = specOf(ns, es, gn, hid);
    ns = f.nodes; es = f.edges; gn = f.grpNames; hid = f.hiddenG;
  }
  return specOf(ns, es, gn, hid);
}'''
text=text.replace(found.group(0),fixed,1)
after_blob=embedded_blob(text);assert hashlib.sha256(after_blob.encode()).hexdigest()==before
tmp=root/'registered_editor.tmp';tmp.write_text(text,encoding='utf-8');os.replace(tmp,path)
assert hashlib.sha256(canonical.read_bytes()).hexdigest()==canonical_hash
report=dict(original_editor=str(path),template_keys=list(specs),old_embedded_library_unchanged=True,canonical_connectome_unchanged=True,storage_key_unchanged='cme_models_v1',live_tab_verified=False,backup=str(backup))
(root/'registration.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
