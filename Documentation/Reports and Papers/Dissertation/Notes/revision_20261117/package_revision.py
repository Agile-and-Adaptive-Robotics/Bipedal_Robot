"""Package the final working sources, compiled proof, and original figure files after verification."""
import argparse
import hashlib
import json
import shutil
import zipfile
from pathlib import Path

notes=Path(__file__).resolve().parent
dis=notes.parents[1]
source=dis/"Final"
figures=dis/"Figures"
parser=argparse.ArgumentParser()
parser.add_argument("build",type=Path)
parser.add_argument("live_text_directory",type=Path)
args=parser.parse_args()
build=args.build.resolve()
def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
def pack(path,members):
    with zipfile.ZipFile(path,"w",compression=zipfile.ZIP_DEFLATED,compresslevel=6) as out:
        for name,p in sorted(members.items()):
            info=zipfile.ZipInfo(name,(2026,11,17,0,0,0))
            info.compress_type=zipfile.ZIP_DEFLATED
            out.writestr(info,p.read_bytes())
    with zipfile.ZipFile(path) as check:
        assert check.testzip() is None
        assert set(check.namelist())==set(members)
        for name,p in members.items():
            assert hashlib.sha256(check.read(name)).hexdigest()==digest(p)
    return {"sha256":digest(path),"members":len(members),"bytes":path.stat().st_size}

members={}
for p in source.iterdir():
    if p.is_file(): members[p.relative_to(dis).as_posix()]=p
for p in (source/"chapters").glob("*.tex"):
    members[p.relative_to(dis).as_posix()]=p
for line in (build/"main.fls").read_text(encoding="utf-8",errors="replace").splitlines():
    if not line.startswith("INPUT "):continue
    path=(source/line[6:].strip().strip('"')).resolve()
    if path.is_relative_to(figures) and path.is_file():
        members[path.relative_to(dis).as_posix()]=path
for chapter_dir in ["15-background","20-muscles","30-torque","50-database","60-sns","70-overall","94-AppendixC","95-AppendixD"]:
    assert any(k.startswith("Figures/"+chapter_dir+"/") for k in members), chapter_dir
assert "Figures/70-overall/bpa_construction.pdf" in members
for chapter_dir,pattern in [("70-overall","bpa_*.alt.txt"),("94-AppendixC","extensor_cam_*.alt.txt")]:
    for p in (figures/chapter_dir).glob(pattern):
        members[p.relative_to(dis).as_posix()]=p
result={}
name="Bolen_Dissertation_revised_sources.zip"
result[name]=pack(notes/name,members)
name="Bolen_Dissertation_revised.pdf"
shutil.copy2(build/"main.pdf",notes/name)
assert digest(notes/name)==digest(build/"main.pdf")
result[name]={"sha256":digest(notes/name),"bytes":(notes/name).stat().st_size}
captured={p.name:p for p in args.live_text_directory.iterdir() if p.is_file()}
assert len(captured)==23
name="LIVE_SOURCE_TEXT.zip"
result[name]=pack(notes/name,captured)
old=figures/"70-overall"/"Deprecated"
preserved={p.relative_to(dis).as_posix():p for p in old.rglob("*") if p.is_file() and
           (p.name.lower().startswith("bpa_construction") or p.name.lower().startswith("bpa_exploded_cutaway"))}
assert len(preserved)>=13
manifest=old/"20261117/manifest.json"
preserved[manifest.relative_to(dis).as_posix()]=manifest
name="SUPERSEDED_FIGURE_SOURCES.zip"
result[name]=pack(notes/name,preserved)
cam_archive=notes/"CAM_FIGURE_ORIGINALS.zip"
if cam_archive.exists():
    with zipfile.ZipFile(cam_archive) as cam_sources:
        assert cam_sources.testzip() is None
        result[cam_archive.name]={"sha256":digest(cam_archive),"members":len(cam_sources.namelist()),"bytes":cam_archive.stat().st_size}
(notes/"DELIVERY_HASHES.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
print(json.dumps(result,indent=2))
