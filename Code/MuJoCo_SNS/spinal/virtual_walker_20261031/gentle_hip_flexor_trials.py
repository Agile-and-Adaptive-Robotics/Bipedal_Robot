"""Two bilateral anatomically verified hip-flexor-only reductions atop hip0.6."""
import hashlib,json
from import_virtual_walker import HERE
from run_virtual_walker import run
from assess_trial import assess
cfg=json.loads((HERE/'GENTLE_HIP_0P6.json').read_text());expected=cfg.pop('expected_source_sha256');assert hashlib.sha256(open(cfg['source'],'rb').read()).hexdigest()==expected
ledger=json.loads((HERE/'runs/canonical_12II/BUILD_LEDGER.json').read_text());bindings=[x for x in ledger['mapping'] if x['joint']=='hip' and x['action']=='flexor'];assert len(bindings)==2 and all(x['moment_arm_m']>0 and x['MN_name'].endswith(' Hip MN flx') for x in bindings)
rows=[]
for factor in (.67,.5):
 tag='canonical108_ground_hip0p6_flex_'+str(factor).replace('.','p')+'_v1';r=run(newer=True,duration=10.,tag=tag,hip_flexor_pf_mn_scale=factor,**cfg);folder=HERE/'runs/W2L_VirtualWalker_Newer'/tag;a=assess(folder);out=dict(tag=tag,hip_factor=.6,hip_flexor_factor=factor,actual_flexor_bindings=bindings,result=r,assessment=a);rows.append(out);(HERE/'runs/GENTLE_HIP_FLEXOR_TRIALS.json').write_text(json.dumps(rows,indent=2));print(json.dumps(dict(factor=factor,gate=a['alternating_supported_stepping_gate'],swings=[a['legs'][s]['unloaded_clear_swings'] for s in 'LR'],alternation=a['alternating_transition_fraction'],support=r['mean_net_vertical_assistance_BW'],relative_progress=r['root_forward_relative_to_belt_m'],failure=r['failure']),indent=2),flush=True)
