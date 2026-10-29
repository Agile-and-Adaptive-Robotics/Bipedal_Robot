"""Fetch small, pinned public source files; keep data in ignored runs/."""
from pathlib import Path
import urllib.request,json,hashlib,concurrent.futures
O=Path(__file__).resolve().parent/'runs/ik_sources';O.mkdir(parents=True,exist_ok=True)
def get(url):
    req=urllib.request.Request(url,headers={'User-Agent':'Gait2392-IK-audit'})
    return urllib.request.urlopen(req,timeout=40).read()
repo='antoinefalisse/predictsim_mtp'
info=json.loads(get('https://api.github.com/repos/'+repo))
branch=info['default_branch'];tree=json.loads(get('https://api.github.com/repos/'+repo+'/git/trees/'+branch+'?recursive=1'))
(O/'falisse_tree.json').write_text(json.dumps(tree,indent=2))
sha=tree['sha'];paths=[x['path'] for x in tree['tree'] if x['type']=='blob' and (
    x['path'].endswith(('motion.mot','GRF.mot')) or
    x['path'] in ('settings.py','README.md') or
    x['path'].endswith('new_model_scaled_FK_contacts.osim'))]
def fetch(path):
    url='https://raw.githubusercontent.com/'+repo+'/'+sha+'/'+urllib.parse.quote(path)
    data=get(url);target=O/'falisse'/path;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
    return dict(path=str(target),url=url,sha256=hashlib.sha256(data).hexdigest(),bytes=len(data))
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:records=list(pool.map(fetch,paths))
(O/'falisse_manifest.json').write_text(json.dumps(dict(repo=repo,commit=sha,files=records),indent=2))
for row in records:print(row['bytes'],row['path'])
