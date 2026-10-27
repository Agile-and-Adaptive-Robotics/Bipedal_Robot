"""Verify the published OpenSim model, independently of the MuJoCo port."""
import os
import sys
from pathlib import Path
import json
import re

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[3]
SOURCE=REPO/'tmp/codex_dirusso_original'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import opensim as osim

path=SOURCE/'import/models/gait9dof18musc_Thelen_20170320.osim'
credits=re.search(r'<credits>(.*?)</credits>',path.read_text(),re.S).group(1)
model=osim.Model(str(path))
model.setUseVisualizer(False)
state=model.initSystem()
coords=model.getCoordinateSet()
muscles=model.getMuscles()
report=dict(path=str(path),opensim_version=osim.GetVersionAndDate(),
            coordinates=[coords.get(i).getName() for i in range(coords.getSize())],
            locked_coordinates=[coords.get(i).getName() for i in range(coords.getSize())
                                if coords.get(i).getLocked(state)],
            muscles=[muscles.get(i).getName() for i in range(muscles.getSize())],
            muscle_classes=sorted(set(muscles.get(i).getConcreteClassName()
                                       for i in range(muscles.getSize()))),
            credits=credits,
            runtime_system_initialized=True,
            scope='Model-load and system-initialization check; not a gait simulation')
(HERE/'author_model_inventory.json').write_text(json.dumps(report,indent=2))
print(json.dumps({k:v for k,v in report.items() if k!='credits'},indent=2))
