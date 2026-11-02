"""Curate template visibility only; never rebuild graphs or touch browser storage."""
from pathlib import Path
import hashlib,json,re

HERE=Path(__file__).resolve().parent
REFERENCE={
 'blank':('Blank editable model','workspace'),
 'rules':('Ben circuit rules — reference (2026-09-24)','reference'),
 'synergy6':('Six-component NMF — reduced data reference','reference'),
 'rybak':('Rybak RG/PF — literature-rule draft','reference'),
 'shevtsova':('Shevtsova RG — literature-rule reference','reference'),
 'shinohara':('Shinohara feedback — literature-rule reference','reference'),
}
BUILTINS=('deng','lit','bilateral','w2l','w2lfile','blank')

def write_if_changed(path,text):
 # Avoid rewriting an open editor unnecessarily, and preserve its original
 # bytes if Windows temporarily locks the destination during replacement.
 if path.exists() and path.read_text(encoding='utf-8')==text:return
 staged=path.with_name(path.name+'.catalog-tmp')
 staged.write_text(text,encoding='utf-8')
 try:staged.replace(path)
 finally:
  if staged.exists():staged.unlink()

def build_catalog(data,extra_keys=()):
 entries={}
 for key in sorted(set(data)|set(BUILTINS)|set(extra_keys)):
  label,state=REFERENCE.get(key,(key.replace('_',' ')+' — archived snapshot','legacy'))
  if key.startswith(('walker_v9','walker_v10')):reason='Superseded v9/v10 parameters do not walk under current physics (AGENTS.md).'
  elif '_flat' in key:reason='Pre-schematic duplicate preserved for historical inspection.'
  elif '_pruned' in key:reason='Historical pruning experiment; not a current runtime graph.'
  elif key.startswith('hunt_session1'):reason='Earlier Hunt/RF-HAM snapshot; not the current executed graph.'
  elif state=='legacy':reason='Dated or illustrative template; not a verified current executed configuration.'
  else:reason='Reusable reference only; not a claim of current controller execution.'
  entries[key]={'label':label,'state':state,'default_visible':state!='legacy','reason':reason}
 return {'schema_version':1,'default_policy':'Only blank and explicitly labelled reusable references. Archived graph data remain available with Legacy enabled.','storage_policy':'No localStorage reads, writes, migration, deletion, or template auto-loading.','current_runtime_policy':'Current graph identity must come from its executed export and parameter/source manifest, never an older template label.','entries':entries,'runtime_import_required':['s3k current saved configuration','w2lvar current saved configuration','syn6 current saved configuration','biw2lvar-Hunt-PT / RF-HAM current saved configuration','W2L108 graded / W2L704 spiking hybrid','rat134 engineering condition','CAD54-actuator /262-neuron configuration','Gait central spiking current saved configuration']}

def refresh_catalog():
 sidecar=HERE/'connectome_templates.json';before=hashlib.sha256(sidecar.read_bytes()).hexdigest();data=json.loads(sidecar.read_text(encoding='utf-8'));htmlpath=HERE/'connectome_block_editor.html';html=htmlpath.read_text(encoding='utf-8');extra=re.findall(r'"(hunt_session1[^"\s]+)"\s*:',html);catalog=build_catalog(data,extra)
 write_if_changed(HERE/'connectome_template_catalog.json',json.dumps(catalog,indent=2,ensure_ascii=False)+'\n')
 blob='/*TPL-CATALOG-START*/\nwindow.TPL_CATALOG_EMBED = '+json.dumps(catalog,ensure_ascii=False,separators=(',',':'))+';\n/*TPL-CATALOG-END*/'
 pattern=r'/\*TPL-CATALOG-START\*/.*?/\*TPL-CATALOG-END\*/'
 if re.search(pattern,html,re.S):html=re.sub(pattern,lambda _:blob,html,flags=re.S)
 else:html=html.replace('/*TPL-EMBED-END*/','/*TPL-EMBED-END*/\n'+blob,1)
 write_if_changed(htmlpath,html);assert hashlib.sha256(sidecar.read_bytes()).hexdigest()==before
 return catalog

if __name__=='__main__':
 c=refresh_catalog();print(json.dumps({'templates':len(c['entries']),'default_visible':[k for k,v in c['entries'].items()if v['default_visible']],'legacy':sum(v['state']=='legacy'for v in c['entries'].values()),'graph_data_rebuilt':False,'browser_storage_touched':False},indent=2))
