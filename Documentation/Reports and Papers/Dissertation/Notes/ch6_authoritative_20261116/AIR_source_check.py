"""Check archived Git blobs against the air recording's declared launch hash."""
from pathlib import Path
import subprocess
import hashlib
import json

repo=Path(r'D:\GitHub\Bipedal_Robot')
out=Path(__file__).resolve().parent
source_path='Code/MuJoCo_SNS/spinal/controller_profiles_20261104/body_probe.py'
target='d73deeebf3813b428e8833b33b075cb15d9e38f8853bcc990b18547d2d324d5a'
revisions=subprocess.check_output(['git','log','--all','--format=%H','--',source_path],cwd=repo).decode().splitlines()
records=[]
for revision in revisions:
    raw=subprocess.check_output(['git','show',f'{revision}:{source_path}'],cwd=repo)
    digest=hashlib.sha256(raw).hexdigest()
    records.append({'revision':revision,'path':source_path,'sha256':digest,'matches_capture':digest==target})
    if digest==target:
        (out/'AIR_launch_body_probe.py').write_bytes(raw)
(out/'AIR_SOURCE_CHECK.json').write_text(json.dumps({'target_sha256':target,'git_blobs_checked':records},indent=2),encoding='utf-8')
print(json.dumps(records,indent=2))
