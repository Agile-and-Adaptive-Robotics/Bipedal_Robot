"""Persistent scientific worker. Adapter handle_request owns model-specific gates."""
import sys,json,argparse,traceback,contextlib,importlib.util
from pathlib import Path

def serve(adapter_path):
 # Keep one scientific interpreter and its imports/JIT cache alive across trials.
 sys.path.insert(0,str(adapter_path.parent))
 if not adapter_path.stem.isidentifier():raise ValueError('Adapter filename must be a valid Python module name for persistent JIT caching')
 spec=importlib.util.spec_from_file_location(adapter_path.stem,adapter_path)
 adapter=importlib.util.module_from_spec(spec)
 sys.modules[spec.name]=adapter  # Required for Numba's cached environment reload.
 with contextlib.redirect_stdout(sys.stderr):spec.loader.exec_module(adapter)
 for raw in sys.stdin:
  try:
   request=json.loads(raw)
   if request.get('op')=='shutdown':print(json.dumps({'ok':True,'shutdown':True}),flush=True);return
   with contextlib.redirect_stdout(sys.stderr):answer=adapter.handle_request(request)
   if not isinstance(answer,dict):raise TypeError('Adapter must return a JSON object')
   print(json.dumps(answer,allow_nan=False),flush=True)
  except Exception as e:
   traceback.print_exc(file=sys.stderr);print(json.dumps(dict(ok=False,error_type=type(e).__name__,error=str(e))),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--adapter',required=True);a=p.parse_args();serve(Path(a.adapter).resolve())
