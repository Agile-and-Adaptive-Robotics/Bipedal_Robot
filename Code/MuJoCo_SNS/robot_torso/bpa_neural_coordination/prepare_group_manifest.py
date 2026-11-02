import json
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa')
cfg=dict(adapter_path=str(HERE/'group_standing_objective.py'),output_dir=str(ROOT/'rustuna_w2lvar_groups_v1'),campaign_root=str(ROOT),seed=100032,max_trials=12,n_startup_trials=4,wall_budget_s=600,trial_timeout_s=150,search_space={k:dict(low=0.,high=1.) for k in ('knee_ext_iiin','hip_ext_iiin')},fixed_parameters={},initial_candidates=[dict(knee_ext_iiin=a,hip_ext_iiin=b) for a,b in ((1.,1.),(0.,1.),(1.,0.),(0.,0.))],scope='Causal two-group IIIN standing calibration after source-current audit and single-group ablations. Preserves source graph; overrides two target-group conductance classes only, tied bilaterally. Global PF/Ib/bias unchanged. No walking/BPA claim.')
p=ROOT/'GAIT_W2LVAR_GROUP_MANIFEST.json'
if p.exists() and json.loads(p.read_text())!=cfg:raise FileExistsError(p)
p.write_text(json.dumps(cfg,indent=2),encoding='utf-8');print(p)
