import json
from pathlib import Path
import dirusso_recovery as R
HERE=Path(__file__).resolve().parent
best=json.loads((HERE/'dirusso_pso_best.json').read_text())
print(json.dumps(R.simulate(best['params'],support=.9,duration=10.,
                           save='dirusso_pso_warning_checked'),indent=2))
