"""Repackage only the HTML after a UI change; retain the built scene data."""
from pathlib import Path
import base64,gzip
p=Path(__file__).resolve();output=p.parents[1]/'packaging/output'
data=base64.b64encode(gzip.compress((output/'scene_data.json').read_bytes())).decode()
template=(p.parent/'packaging_view.html').read_text(encoding='utf-8')
(output/'AH001_Packaging.html').write_text(template.replace('__DATA__',data),encoding='utf-8')
