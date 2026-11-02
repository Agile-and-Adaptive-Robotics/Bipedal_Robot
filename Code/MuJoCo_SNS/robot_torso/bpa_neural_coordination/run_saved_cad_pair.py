"""Run the frozen2s literal-SNS pressurepair in a new outputfolder safely."""
import argparse,json,hashlib,platform,shutil
from pathlib import Path
BASE=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\cad_bpa\pair_neural_recruitment\normalized_unclamped');CAD_ROOT=BASE.parents[1]
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True,help='New empty absolute folder under the authorized external cad_bpa directory');p.add_argument('--mode',choices=['reference','mixed','both'],default='both');args=p.parse_args()
 import cad_pair_unclamped as U
 if not args.output.is_absolute():raise ValueError('Supply a full absolute output path')
 dest=args.output.resolve()
 if platform.node().lower()!='easteregg2':raise ValueError('This exact saved configuration is prepared for easteregg2; do not substitute another machine runtime silently')
 if not dest.is_relative_to(CAD_ROOT.resolve()) or dest==CAD_ROOT.resolve():raise ValueError('Output must be an owned subfolder under '+str(CAD_ROOT))
 if dest.exists() and any(dest.iterdir()):raise ValueError('Choose a new empty output folder; retained evidence is never overwritten')
 protocol=json.loads((BASE/'RECRUITMENT_PREREGISTRATION.json').read_text())
 for path,expected in protocol['source_sha256'].items():
  if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=expected:raise ValueError('Frozen scientific source changed: '+path)
 core=json.loads((BASE/'mixed/cad.json').read_text())
 for relative,expected in core['source_hashes'].items():
  path=U.C.ROBOT.parent/relative
  if hashlib.sha256(path.read_bytes()).hexdigest()!=expected:raise ValueError('Original model/controller dependency changed: '+str(path))
 if (BASE/'ASSUMED_PAIR_CONFIGURATION.json').read_bytes()!=(U.A.OLD/'ASSUMED_PAIR_CONFIGURATION.json').read_bytes():raise ValueError('Pressurefixture differs from frozen original')
 dest.mkdir(parents=True,exist_ok=True);U.OUT=dest
 if args.mode=='mixed':
  folder=dest/'reference';folder.mkdir();reference=json.loads((BASE/'reference/BRIDGE_SUMMARY.json').read_text());reference.update(reference_fixture_only_not_new_run=True,original_reference_summary=str(BASE/'reference/BRIDGE_SUMMARY.json'));(folder/'BRIDGE_SUMMARY.json').write_text(json.dumps(reference,indent=2))
 modes=('reference','mixed') if args.mode=='both' else (args.mode,)
 for mode in modes:U.run(mode)
 print('Saved new known-nominal2s study: '+str(dest))
if __name__=='__main__':main()
