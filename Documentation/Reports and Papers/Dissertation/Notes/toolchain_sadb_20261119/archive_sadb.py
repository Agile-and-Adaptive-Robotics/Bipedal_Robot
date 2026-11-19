"""Preserve the four earlier SADb graphics with their hash and original chapter path."""
from pathlib import Path
import hashlib
import json
from zipfile import ZipFile, ZIP_DEFLATED

here=Path(__file__).resolve().parent
diss=here.parents[1]
root=(diss/"Figures/50-database").resolve()
archive=(root/"Deprecated/20261119_toolchain_sadb").resolve()
assert archive.is_relative_to(root/"Deprecated")
archive.mkdir(exist_ok=True,parents=True)
active_sources="\n".join(p.read_text(encoding="utf8") for p in (diss/"Final/chapters").glob("*.tex"))
old=["sadb_clusters_app.png","sadb_pivot_app.png","sadb_table_detail.jpg","sadb_table_record.jpg"]
manifest={}
for name in old:
    src=(root/name).resolve()
    dst=(archive/name).resolve()
    assert src.parent==root and dst.parent==archive and name not in active_sources
    data=(src if src.exists() else dst).read_bytes()
    if src.exists():
        assert not dst.exists() or dst.read_bytes()==data
        src.replace(dst)
    assert dst.read_bytes()==data
    manifest[name]={"former":src.relative_to(diss).as_posix(), "archived":dst.relative_to(diss).as_posix(),
                    "sha256":hashlib.sha256(data).hexdigest()}
description=json.dumps(manifest,indent=2)+"\n"
(archive/"MANIFEST.json").write_text(description,encoding="utf-8")
out=here/"SADB_FIGURE_SOURCES.zip"
with ZipFile(out,"w",compression=ZIP_DEFLATED) as z:
    for path in [*sorted(archive.glob("*.*")),*sorted((root/"components/sadb_app_20261119").glob("*.*")),here/"sadb_figures.py"]:
        z.write(path,path.relative_to(diss).as_posix())
with ZipFile(out) as z:
    assert z.testzip() is None
    for rel in z.namelist():
        assert (diss/rel).read_bytes()==z.read(rel)
print("Preserved", len(old), "prior chapter graphics and the source app captures.")
