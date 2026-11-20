"""Transactional native model packages and separate editor presentation files."""
import copy,json,os,tempfile
from pathlib import Path
import native_w2l as native
from project import ProjectError,fingerprint

def atomic_package(outputs):
    """Stage every file first; restore previous bytes if replacement fails."""
    staged={};prior={}
    try:
        for path,data in outputs.items():
            path.parent.mkdir(parents=True,exist_ok=True);prior[path]=path.read_bytes() if path.exists() else None
            fd,tmp=tempfile.mkstemp(prefix='.'+path.name+'.',suffix='.tmp',dir=path.parent)
            staged[path]=tmp
            with os.fdopen(fd,'wb') as f:f.write(data);f.flush();os.fsync(f.fileno())
        written=[]
        try:
            for path,tmp in staged.items():os.replace(tmp,path);written.append(path)
        except Exception:
            for path in reversed(written):
                if prior[path] is None:path.unlink(missing_ok=True)
                else:path.write_bytes(prior[path])
            raise
    finally:
        for tmp in staged.values():
            if os.path.exists(tmp):os.unlink(tmp)

def save_model(doc,destination,source=None):
    doc=copy.deepcopy(doc);metadata=doc['native'];asim,project=native.build_bytes(doc)
    source=Path(source) if source else None
    same=source==destination
    if same:
        apath=destination.parent/metadata['asim_name'];ppath=destination.parent/metadata['project_name'] if project else None
    else:
        apath=destination.with_suffix('.asim');ppath=destination.with_suffix('.aproj') if project else None
        for target in (apath,ppath):
            if target and target!=destination and target.exists():raise ProjectError('A companion model already exists at '+str(target)+'. Choose another name to keep both copies.')
    outputs={apath:asim}
    if project:
        root=native.parse(project);before=native.E.tostring(root)
        if not same:
            native.field(root,'ProjectName',destination.stem);native.field(root,'SimulationFile',apath.name)
            for item in root.iter('ProjectPath'):item.text=str(destination.parent)+os.sep
        if root.findtext('SimulationFile')!=apath.name:native.field(root,'SimulationFile',apath.name)
        if native.E.tostring(root)!=before:project=native.xml(root)
        outputs[ppath]=project
    # The model may reference local chart forms and textures. Copy those files
    # with Save As, without replacing other assets already in the target folder.
    if not same:
        for old in metadata.get('source_files',{}):
            for asset in Path(old).parent.iterdir():
                if asset.suffix.lower() not in ('.aform','.bmp','.png','.jpg'):continue
                target=destination.parent/asset.name;data=asset.read_bytes()
                if target.exists() and target.read_bytes()!=data:raise ProjectError('A different asset already exists: '+str(target))
                outputs[target]=data
        assets={p.name for p in outputs if p.suffix.lower() in ('.aform','.bmp','.png','.jpg')}
        for path in (apath,ppath):
            if path is None:continue
            root=native.parse(outputs[path]);before=native.E.tostring(root)
            for item in root.iter():
                if len(item)==0 and item.text and Path(item.text).name in assets:item.text=str(destination.parent/Path(item.text).name)
            if native.E.tostring(root)!=before:outputs[path]=native.xml(root)
        asim=outputs[apath];project=outputs[ppath] if ppath else None
    adapter_ids={x.findtext('ID') for x in native.parse(asim).iter('Adapter')};metadata['adapters']=[a for a in metadata['adapters'] if a['id'] in adapter_ids]
    metadata.update(asim_name=apath.name,asim_b64=native.enc(asim),asim_sha256=native.sha(asim),project_name=ppath.name if ppath else None,project_b64=native.enc(project) if project else None,project_sha256=native.sha(project) if project else None,primary_suffix=destination.suffix.lower())
    metadata['source_files']={str(p):native.sha(data) for p,data in outputs.items() if p.suffix.lower() in ('.asim','.aproj')}
    doc['source']['file']=str(destination)
    sidecar=destination.with_suffix('.editor.json');outputs[sidecar]=json.dumps({'schema':'sns-native-presentation/1','document':doc},ensure_ascii=False,indent=2,allow_nan=False).encode('utf-8')
    atomic_package(outputs)
    return doc,{str(p):fingerprint(p) for p in outputs if p.suffix.lower() in ('.asim','.aproj') or p==sidecar}
