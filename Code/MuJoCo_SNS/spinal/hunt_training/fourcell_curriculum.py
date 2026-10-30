"""Physical checks of the freshly refitted four-cell candidate; no silent graduation."""
import json
from pathlib import Path
import balance_campaign as B
H=Path(__file__).resolve().parent;root=H/'runs/fourcell_hunt_20261029_v1'
assert (root/'TRAINING_COMPLETE.json').exists()
B.O=root/'curriculum';B.O.mkdir(exist_ok=True)
base=dict(weights_directory=str(root/'weights'),biarticular_pair_gain=.3,commissural_gain=.3,learn_reciprocal=True)
results={}
def run(name,changes):
    r=B.run(name,{**base,**changes});results[name]=r
    if r:print(name,json.dumps(r['assessment']),flush=True)
    return r
run('stage1',dict(stage=1))
plain=run('stage2_plain',dict(stage=2))
axis=run('stage2_axis',dict(stage=2,axis_feedback=1.5,hip_abduction_goal=.08))
scaled=run('stage2_scaled_probe',dict(stage=2,axis_feedback=1.5,hip_abduction_goal=.08,motor_scales={'knee_ext_r':.3,'knee_ext_l':.3,'knee_flex_r':2.,'knee_flex_l':2.},calibration_probe=True))
accepted=[r for r in (plain,axis,scaled) if r and r['assessment']['accepted']]
candidates=[r for r in (plain,axis,scaled) if r]
winner=max(accepted or candidates,key=lambda r:r['assessment']['score'])
stand={**winner['configuration'],'stage':3,'restriction':1.,'height_offset':-.02,'start_hip_deg':-7.,'start_knee_deg':-10.,'start_ankle_deg':0.,'balance':dict(kp=.5,kd=.2,ka=0.,delay_s=.1,replace_inherited=False)}
run('stage3',stand)
for height in (-.02,0.,-.04):run('stage4_height_'+str(height),{**stand,'stage':4,'height_offset':height})
summary=dict(training=str(root/'TRAINING_COMPLETE.json'),results={k:dict(configuration=v['configuration'],assessment=v['assessment'],measured=v['measured']) for k,v in results.items() if v},
    note='Scaled output case is an explicit physical calibration probe. No later stage graduates if prerequisites fail. Stage 5 is withheld until loaded alternating ground steps pass.')
(B.O/'summary.json').write_text(json.dumps(summary,indent=2))
