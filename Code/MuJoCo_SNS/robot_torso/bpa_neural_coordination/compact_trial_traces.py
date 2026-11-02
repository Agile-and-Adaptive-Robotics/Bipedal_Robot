"""RAM-only full telemetry for gate calculations, sparse/windowed saved witnesses.

Only intercepts newly owned output paths beneath an explicit case directory.
Source/model archives elsewhere retain their normal NumPy IO behavior.
"""
from pathlib import Path
import hashlib,json
from collections.abc import Mapping
import numpy as np

class Archive(Mapping):
 def __init__(self,data):self.data=data;self.files=list(data)
 def __getitem__(self,key):return self.data[key].copy()
 def __iter__(self):return iter(self.data)
 def __len__(self):return len(self.data)
 def __enter__(self):return self
 def __exit__(self,*args):return False
 def close(self):pass

class CompactTraces:
 def __init__(self,folder):self.folder=Path(folder).resolve();self.cache={};self.paths=set()
 def _path(self,file):
  try:return Path(file if isinstance(file,(str,bytes,Path)) else file.name).resolve()
  except (TypeError,ValueError):return None
 def __enter__(self):
  self.original_save=np.savez_compressed;self.original_save_plain=np.savez;self.original_load=np.load
  def save(file,*args,**kwargs):
   path=self._path(file)
   if path is None or not path.is_relative_to(self.folder):return self.original_save(file,*args,**kwargs)
   if args:raise ValueError('Compact case telemetry requires named arrays')
   data={k:np.asarray(v).copy() for k,v in kwargs.items()};key=path
   if path.name.endswith('.tmp.npz') and '.npz.' in path.name:key=path.with_name(path.name.split('.npz.',1)[0]+'.npz')
   self.cache[key]=data;self.paths.add(path);self.paths.add(key)
  def load(file,*args,**kwargs):
   path=self._path(file)
   if path in self.cache:return Archive(self.cache[path])
   return self.original_load(file,*args,**kwargs)
  np.savez_compressed=save;np.savez=save;np.load=load;return self
 def __exit__(self,*args):
  np.savez_compressed=self.original_save;np.savez=self.original_save_plain;np.load=self.original_load
  # runner.py creates empty temporary files before atomic replacement. They
  # contain no scientific data; remove only verified owned zero-byte outputs.
  for path in self.paths:
   target=path.resolve()
   if not target.is_relative_to(self.folder):raise ValueError('Refuse cleanup outside owned case')
   if target.is_file() and target.stat().st_size==0:target.unlink()
  return False
 def materialize(self,full=False):
  manifest=[];compact={}
  for index,(path,data) in enumerate(sorted(self.cache.items(),key=lambda kv:str(kv[0]))):
   entries={}
   for key,a in data.items():
    digest=hashlib.sha256(a.tobytes()).hexdigest();entries[key]=dict(shape=list(a.shape),dtype=str(a.dtype),sha256=digest)
    label=f'a{index}_{key}'
    # Short state snapshots are preserved; long history/logger arrays are
    # sampled uniformly plus startup/end windows. No hidden pass recomputation.
    if a.ndim and len(a)>512:
     take=np.unique(np.r_[np.arange(0,len(a),100),np.arange(min(16,len(a))),np.arange(max(0,len(a)-128),len(a))]);value=a[take]
     # Do not duplicate all888neural states in exploratory witnesses.
     if key=='neuro' and value.ndim==2 and value.shape[1]>100:value=value[:,:32]
     compact[label]=value;compact[label+'_sample_indices']=take
    elif a.size<=20000:compact[label]=a
   if full:self.original_save(path,**data)
   manifest.append(dict(original_path=str(path),full_trace_retained=full,arrays=entries))
  self.original_save(self.folder/'COMPACT_TELEMETRY.npz',**compact)
  (self.folder/'TELEMETRY_MANIFEST.json').write_text(json.dumps(dict(policy='Full gate calculations used RAM arrays. Exploratory witnesses are sparse/windowed; source/model files never intercepted.',full_retained=full,archives=manifest),indent=2),encoding='utf-8')
