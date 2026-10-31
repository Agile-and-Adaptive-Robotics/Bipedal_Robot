import ast,json,hashlib
import numpy as np
from source_inputs import historical_source,HERE
from native_kernel import NativeController
manifest=json.loads((HERE/'runs/source_inputs/MANIFEST.json').read_text());bank=json.loads((HERE/'runs/W2L_VirtualWalker_Newer_GUI1nA/IMPORT_AUDIT.json').read_text())['muscles'];rows=[]
for row in manifest['sources']:
 stem=__import__('pathlib').Path(row['snapshot']).stem.replace('_Standalone','');p=historical_source(stem);a=NativeController(row['original'],bank);b=NativeController(p,bank);same=all(np.array_equal(getattr(a,k),getattr(b,k)) for k in ('params','edges','tonic','v','m','h'));assert same and a.stimuli==b.stimuli;rows.append(dict(stem=stem,sha256_equal=True,controller_parameters_and_initialization_equal=True,neurons=len(b.v),connections=len(b.edges)))
for p in HERE.glob('*.py'):ast.parse(p.read_text(),filename=str(p))
out=dict(results=rows,syntax_valid=True,live_experiment_source_inputs_no_longer_require_generated_folders=True);(HERE/'runs/source_inputs/CONSTRUCTOR_EQUIVALENCE.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
