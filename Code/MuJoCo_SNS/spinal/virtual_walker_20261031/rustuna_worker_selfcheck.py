"""Protocol-only persistent-worker test; no optimizer or scientific model."""
import json,time
from rustuna_transfer import Worker,HERE
out=HERE/'runs/rustuna_protocol_probe';out.mkdir(parents=True,exist_ok=True)
w=Worker(HERE/'rustuna_protocol_probe.py',out)
try:
 preflight=w.request({'op':'preflight','manifest':{}},10.);assert preflight['ready_for_optimization'] is False
 rows=[];started=time.perf_counter()
 for i in range(100):
  x=i*.05;a=w.request(dict(op='evaluate',params={'x':x}),5.);assert a['value']==(x-2)**2;rows.append(a)
 elapsed=time.perf_counter()-started;assert len({x['diagnostics']['pid'] for x in rows})==1
 result=dict(probe='Toy protocol only; no optimization or MuJoCo simulation',preflight=preflight,requests=100,total_request_wall_s=elapsed,ms_per_round_trip=elapsed/100*1000,same_worker_pid=True,known_objective_exact=True)
 (out/'PROTOCOL_SELFCHECK.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
finally:w.close()
