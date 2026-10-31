"""Immutable inputs for historical experiments, not active walker families."""
from pathlib import Path
import hashlib,json
HERE=Path(__file__).resolve().parent
def historical_source(stem):
 path=HERE/'runs/source_inputs'/f'{stem}_Standalone.asim'
 if not path.exists():raise FileNotFoundError(f'Historical input is archived, not an active model. Restore {path.name} using runs/source_inputs/MANIFEST.json and the retained archive ZIP.')
 manifest=json.loads((path.parent/'MANIFEST.json').read_text());expected=next(x['sha256'] for x in manifest['sources'] if x['snapshot']==str(path));assert hashlib.sha256(path.read_bytes()).hexdigest()==expected,'Historical source snapshot changed'
 return path
