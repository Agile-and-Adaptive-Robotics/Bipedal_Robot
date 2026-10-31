"""Write consumed-parameter bounded study manifest; does not launch a run."""
import json
from import_virtual_walker import HERE
manifest=dict(adapter_path=str(HERE/'w2l_ground_optimization.py'),output_dir=str(HERE/'runs/rustuna108_ground_home_v1'),search_space=dict(pf_mn_recruitment_scale=dict(low=.8,high=1.5,log=False),toe_mn_gmax=dict(low=0.,high=.05,log=False)),fixed_parameters={},seed=100031,max_trials=24,wall_budget_s=600,trial_timeout_s=90,initial_candidates=[dict(pf_mn_recruitment_scale=1.,toe_mn_gmax=.01),dict(pf_mn_recruitment_scale=.8,toe_mn_gmax=.01),dict(pf_mn_recruitment_scale=1.5,toe_mn_gmax=.01),dict(pf_mn_recruitment_scale=1.,toe_mn_gmax=0.),dict(pf_mn_recruitment_scale=1.,toe_mn_gmax=.05)])
(HERE/'RUSTUNA108_GROUND_HOME_MANIFEST.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(manifest,indent=2))
