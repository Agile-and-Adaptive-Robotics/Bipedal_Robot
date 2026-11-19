"""Package the manuscript and changed artwork using their destinations in Overleaf."""
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
import argparse
from difflib import unified_diff
import hashlib
import json
import re
from capture_baseline import fingerprint

here=Path(__file__).resolve().parent
diss=here.parents[1]
repo=diss.parents[2]
scratch=repo/"tmp/toolchain_sadb_20261119"
figs=diss/"Figures"
final=diss/"Final"
chapters=["50-database","60-sns","70-overall","95-AppendixD","96-AppendixE"]
all_folders=chapters+["92-AppendixA","93-AppendixB","94-AppendixC"]
previous=repo/"tmp/ch6_advisor_20261119/overleaf_source_r9"
capture=json.loads((here/"LIVE_CAPTURE.json").read_text())
with ZipFile(diss/"Bolen_Dissertation.zip") as export:
    raw=export.read("thesis.bib")
lf=raw.replace(b"\r\n",b"\n")
local_bib=(final/"thesis.bib").read_bytes()
entry=(here/"Chomsky_2023.bib").read_bytes()
entry_lf=entry.replace(b"\r\n",b"\n")
additions=(entry_lf,entry_lf.replace(b"\n",b"\r\n"))
prior=next(local_bib.removesuffix(b"\n"+a) for a in additions
           if local_bib.endswith(b"\n"+a) and hashlib.sha256(local_bib.removesuffix(b"\n"+a)).hexdigest()==capture["previous_local_bibliography_sha256"])
assert prior.replace(b"\r\n",b"\n")==lf
assert list(fingerprint(prior.decode("utf8",errors="replace").replace("\r\n","\n")))==capture["bibliography"]["fingerprint"]
assert hashlib.sha256(prior).hexdigest()==capture["previous_local_bibliography_sha256"]
assert local_bib.startswith(prior) and local_bib.count(b"@misc{chomsky_false_promise_2023,")==1
assert local_bib[len(prior):].strip().replace(b"\r\n",b"\n")==entry.strip().replace(b"\r\n",b"\n")
members={"chapters/"+c+".tex":(final/"chapters"/(c+".tex")).read_bytes() for c in chapters}
if (here/"FRONT_MATTER_HASHES.json").exists():
    for name,sha in json.loads((here/"FRONT_MATTER_HASHES.json").read_text()).items():
        contents=(final/name).read_bytes()
        assert hashlib.sha256(contents).hexdigest()==sha
        members[name]=contents
members["Chomsky_2023.bib"]=entry
assets={}
include={}
with ZipFile(here.parent/"ch6_advisor_20261119/Revised_chapters_and_appendices.zip") as z:
    prior_transfer=set(z.namelist())
with ZipFile(diss/"Bolen_Dissertation.zip") as export:
    exported=set(export.namelist())
    for chapter in all_folders:
        folder=figs/chapter
        used=set(re.findall(r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}",(final/"chapters"/(chapter+".tex")).read_text(encoding="utf8")))
        active={p.name for p in folder.iterdir() if p.is_file() and p.suffix.lower() in {".png",".pdf",".jpg",".jpeg",".eps",".svg"}}
        assert active==used,(chapter,active^used)
        include[chapter]=sorted(active)
        changed_images={a for a in active if (key:="figs/"+chapter+"/"+a) not in exported or (folder/a).read_bytes()!=export.read(key)}
        for p in sorted(folder.iterdir()):
            if not p.is_file() or (p.name not in active and not(p.suffix==".txt" and any(p.name in {Path(a).stem+".alt.txt",Path(a).stem+"_alt.txt"} for a in active))):
                continue
            name="figs/"+chapter+"/"+p.name
            if p.suffix==".txt" and name not in prior_transfer and not any(p.name in {Path(a).stem+".alt.txt",Path(a).stem+"_alt.txt"} for a in changed_images):
                continue
            data=p.read_bytes()
            if name not in exported or data!=export.read(name):
                members[name]=data
                assets[name]="updated or new in this pass" if not (previous/name).exists() or (previous/name).read_bytes()!=data else "already in the preceding transfer"
    changes=[]
    for name,data in members.items():
        if not name.endswith((".tex",".cls")):
            continue
        old=export.read(name).decode("utf8").replace("\r\n","\n")
        new=data.decode("utf8").replace("\r\n","\n")
        changes.extend(unified_diff(old.splitlines(keepends=True),new.splitlines(keepends=True),
                                    fromfile="author_export/"+name,tofile="local_revision/"+name))
    (here/"OVERLEAF_SOURCE_CHANGES.diff").write_text("".join(changes),encoding="utf8")
font_repairs=json.loads((here/"FONT_EMBEDDING_FILES.json").read_text())
for relative,hashes in font_repairs.items():
    data=(figs/relative).read_bytes()
    assert hashlib.sha256(data).hexdigest()==hashes["after_sha256"],relative
    name="figs/"+relative
    if name in members:
        assert members[name]==data
    members[name]=data
    assets[name]="font embedding verified against archived artwork"
sources="\n".join("  "+name for name in sorted(members) if name.endswith((".tex",".cls")))
guide=(
    "Updating the existing Bolen-Dissertation project in Overleaf\n\n"
    "1. Apply the supplied source files to the matching existing paths in your project.\n"
    "   The review folder also contains OVERLEAF_SOURCE_CHANGES.diff for use with tracked edits or comments.\n"
    +sources+"\n\n"
    "2. Merge the contents of figs/ into the matching folders beneath figs/ in Overleaf.\n"
    "   The PDF graphics must be uploaded as files. The .txt companions describe the figures.\n"
    "   The package includes font repairs outside Chapters 5–7.\n\n"
    "3. Append the entry in Chomsky_2023.bib to your existing thesis.bib if the key\n"
    "   chomsky_false_promise_2023 is absent. The supplied project already includes ong_predicting_2019.\n\n"
    "4. Recompile main.tex. The Abstract should read i. In Chapter 7, Figure 7.2C is the\n"
    "   exploded BPA and Figure 7.5 shows how to repeat a controlled simulation experiment.\n\n"
    "Locally, four superseded Chapter 5 snapshots are stored in Deprecated. If you keep the same\n"
    "archive in Overleaf, move the following old figs/50-database/ files after updating Chapter 5:\n"
    "  sadb_clusters_app.png\n  sadb_pivot_app.png\n  sadb_table_detail.jpg\n  sadb_table_record.jpg\n")
members["READ_ME_FIRST.txt"]=guide.encode("utf8")
zip_path=here/"Dissertation_Overleaf_updates.zip"
with ZipFile(zip_path,"w",compression=ZIP_DEFLATED) as z:
    for key,value in members.items(): z.writestr(key,value)
with ZipFile(zip_path) as z:
    assert z.testzip() is None and set(z.namelist())==set(members)
    for name,value in members.items(): assert z.read(name)==value
ap=argparse.ArgumentParser()
ap.add_argument("--stage",type=Path)
args=ap.parse_args()
if args.stage:
    stage=args.stage.resolve()
    assert stage.is_relative_to(scratch) and not stage.exists()
    stage.mkdir(parents=True)
    with ZipFile(diss/"Bolen_Dissertation.zip") as z:
        for entry_info in z.infolist():
            if entry_info.is_dir():
                continue
            path=(stage/entry_info.filename).resolve()
            assert path.is_relative_to(stage)
            path.parent.mkdir(parents=True,exist_ok=True)
            path.write_bytes(z.read(entry_info))
    original_bib=(stage/"thesis.bib").read_bytes()
    assert original_bib.replace(b"\r\n",b"\n")==prior.replace(b"\r\n",b"\n")
    with ZipFile(zip_path) as z:
        for name in z.namelist():
            if name in {"main.tex","psutex.cls"} or name.startswith(("chapters/","figs/")):
                path=(stage/name).resolve()
                assert path.is_relative_to(stage)
                path.parent.mkdir(parents=True,exist_ok=True)
                path.write_bytes(z.read(name))
        additive=z.read("Chomsky_2023.bib")
        if prior.endswith(b"\r\n"):
            additive=additive.replace(b"\r\n",b"\n").replace(b"\n",b"\r\n")
        (stage/"thesis.bib").write_bytes(prior+b"\n"+additive)
        assert (stage/"thesis.bib").read_bytes()==local_bib
    for chapter,active in include.items():
        for name in active:
            assert (figs/chapter/name).read_bytes()==(stage/"figs"/chapter/name).read_bytes(),(chapter,name)
    for relative,hashes in font_repairs.items():
        assert hashlib.sha256((stage/"figs"/relative).read_bytes()).hexdigest()==hashes["after_sha256"],relative
item={"source_export_sha256":hashlib.sha256((diss/"Bolen_Dissertation.zip").read_bytes()).hexdigest(),"active_images_by_chapter":include,"total_graphics_audited":sum(map(len,include.values())),"members":{name:hashlib.sha256(data).hexdigest() for name,data in members.items()},"graphic_changes":assets,
      "font_repairs":font_repairs,"current_live_bib_fingerprint":capture["bibliography"]["fingerprint"],"local_bib_preserves_previous_bytes":True}
(here/"TRANSFER.json").write_text(json.dumps(item,indent=2)+"\n",encoding="utf-8")
print("Transfer verified:",len(members),"members,",item["total_graphics_audited"],"directly included images across eight folders.")
