"""Strict complete loader for the normalized standing mapping checkpoint."""
import argparse,json,hashlib
from pathlib import Path
import gait_standing_family as F

def main():
 p=argparse.ArgumentParser();p.add_argument('--checkpoint',required=True);p.add_argument('--out',required=True);p.add_argument('--tag',required=True);p.add_argument('--duration',type=float,default=8.);p.add_argument('--validate-only',action='store_true');a=p.parse_args();config=json.loads(Path(a.checkpoint).read_text());params=config['global_parameters'];scales=config['pathway_group_scales']
 if set(params)!={'pf_mn_gain','ib_group_gain','standing_bias_scale'} or set(scales)!={'IIIN:'+g for g in config['source_fanin']}:raise ValueError('Incomplete or unconsumed candidate parameter set')
 for name,expected in config['expected_source_sha256'].items():
  source=Path(F.__file__) if name==Path(F.__file__).name else F.SPINAL/name
  if hashlib.sha256(source.read_bytes()).hexdigest()!=expected:raise ValueError('Source changed since validated checkpoint; revalidate explicitly: '+str(source))
 for path,expected in config['expected_dependency_sha256'].items():
  if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=expected:raise ValueError('Dependency changed since validation: '+path)
 if a.validate_only:print('Complete saved parameters and source dependencies validated; no simulation started');return
 F.run('w2lvar',Path(a.out),duration=a.duration,tag=a.tag,gain=params['pf_mn_gain'],ib_group_gain=params['ib_group_gain'],standing_bias_scale=params['standing_bias_scale'],pathway_group_scales=scales,capture_neural_inputs=True)

if __name__=='__main__':main()
