"""Cheap cross-interpreter plumbing check only; deliberately refuses optimization."""
import sys,os
def handle_request(request):
 if request['op']=='preflight':return dict(ok=True,ready_for_optimization=False,reason='Protocol probe only; no accepted rat baseline or actuator mapping attached',worker_python=sys.executable,worker_version=sys.version)
 if request['op']=='evaluate':
  x=request['params']['x'];return dict(ok=True,value=(x-2)**2,diagnostics={'toy_only':True,'pid':os.getpid()})
 raise ValueError('Unknown request')
