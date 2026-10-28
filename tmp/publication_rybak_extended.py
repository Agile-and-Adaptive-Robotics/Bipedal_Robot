from pathlib import Path
import subprocess,sys,json,time,os
R=Path('D:/GitHub/Bipedal_Robot');S=R/'Code/MuJoCo_SNS/spinal';O=S/'reports_20261028_publication_flow';os.environ['PYTHONIOENCODING']='utf-8';os.environ['OPENBLAS_NUM_THREADS']='1';t=time.time()
with (O/'rybak_extended.log').open('w',encoding='utf-8') as f:
 try:rc=subprocess.run([sys.executable,'-u',str(S/'replication/verify_rybak2006.py')],cwd=O,stdout=f,stderr=subprocess.STDOUT,timeout=900).returncode
 except subprocess.TimeoutExpired:rc=124
(O/'rybak_extended_status.json').write_text(json.dumps({'returncode':rc,'seconds':time.time()-t,'scope':'fresh existing reduced-circuit test; no parameter changes','prior_timeout_seconds':240},indent=2));print(rc,flush=True)
