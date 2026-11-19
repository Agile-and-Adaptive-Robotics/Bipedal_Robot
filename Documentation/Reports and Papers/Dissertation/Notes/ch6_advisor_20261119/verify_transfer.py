"""Apply the transfer to an isolated copy of the existing Overleaf export, then compile it."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import sys
from zipfile import ZipFile

here=Path(__file__).resolve().parent
diss=here.parents[1]
repo=diss.parents[2]
p=argparse.ArgumentParser()
p.add_argument("target",type=Path)
p.add_argument("--compare-builds",nargs=2,type=Path,metavar=("LOCAL","TRANSFER"))
args=p.parse_args()
if args.compare_builds:
    import pymupdf
    local,transfer=args.compare_builds
    log=(transfer/"main.log").read_text(errors="replace")
    blg=(transfer/"main.blg").read_text(errors="replace")
    assert not re.search(r"^!|LaTeX Warning: .*undefined|Package natbib Warning|multiply defined|Overfull \\[hv]box",log,re.I|re.M)
    assert not re.search(r"^Warning--|error message",blg,re.M)
    with pymupdf.open(local/"main.pdf") as a, pymupdf.open(transfer/"main.pdf") as b:
        compared=dict(local_pages=len(a),transfer_pages=len(b),different_page_text=[])
        start=next(p for level,t,p in a.get_toc() if level==1 and "Tools for Constructing" in t)
        end=next(p for level,t,p in a.get_toc() if level==1 and p>start)
        sb=next(p for level,t,p in b.get_toc() if level==1 and "Tools for Constructing" in t)
        eb=next(p for level,t,p in b.get_toc() if level==1 and p>sb)
        compared["chapter6_matching_pages"]=sum(a[i].get_text()==b[j].get_text() for i,j in zip(range(start-1,end-1),range(sb-1,eb-1)))
        assert compared["chapter6_matching_pages"]==end-start==eb-sb
        compared["different_page_text"]=[i+1 for i in range(min(len(a),len(b))) if a[i].get_text()!=b[i].get_text()]
    compare_path=here/"TRANSFER_CHECK.json"
    evidence=json.loads(compare_path.read_text())
    evidence["local_transfer_compile"]=dict(**compared,undefined_references=0,undefined_citations=0,
        BibTeX_warnings=0,overfull_boxes=0,transfer_sha256=hashlib.sha256((here/"Revised_chapters_and_appendices.zip").read_bytes()).hexdigest())
    compare_path.write_text(json.dumps(evidence,indent=2)+"\n",encoding="utf-8")
    print(compared)
    assert compared["local_pages"] == compared["transfer_pages"] and not compared["different_page_text"]
    sys.exit(0)
target=args.target.resolve()
assert target.is_relative_to(repo/"tmp") and not target.exists()
target.mkdir(parents=True)
def write(destination, name, data):
    path=(destination/name).resolve()
    assert path.is_relative_to(destination)
    path.parent.mkdir(exist_ok=True,parents=True)
    path.write_bytes(data)
export = diss/"Bolen_Dissertation.zip"
assert hashlib.sha256(export.read_bytes()).hexdigest()=="bdf32431a4b823961bad1fa9d9c93976ab0957cc071915ef6302d41d4039b5fc"
with ZipFile(export) as z:
    for member in z.infolist():
        if not member.is_dir():
            write(target,member.filename,z.read(member))
with ZipFile(here/"LIVE_SOURCES.zip") as z:
    for name in z.namelist():
        write(target,"chapters/"+name,z.read(name))
with ZipFile(here/"Revised_chapters_and_appendices.zip") as z:
    assert z.testzip() is None
    for member in z.namelist():
        if member.startswith(("figs/","chapters/")):
            write(target,member,z.read(member))
    bib=target/"thesis.bib"
    original=bib.read_bytes()
    assert b"ong_predicting_2019" not in original
    bib.write_bytes(original+b"\n"+z.read("Ong_2019.bib"))
assert (target/"main.tex").exists()
inventory={}
issues=[]
for chapter in ["60-sns","70-overall","92-AppendixA","93-AppendixB","94-AppendixC","95-AppendixD","96-AppendixE"]:
    figures=diss/"Figures"/chapter
    src=(target/"chapters"/(chapter+".tex")).read_text(encoding="utf-8")
    current=set(re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}",src))
    root_images={x.name for x in figures.iterdir() if x.is_file() and x.suffix.lower() in {".jpg",".jpeg",".pdf",".png",".eps",".svg"}}
    assert root_images==current, (chapter,root_images^current)
    for name in current:
        assert (figures/name).read_bytes()==(target/"figs"/chapter/name).read_bytes(),(chapter,name)
    inventory[chapter]=sorted(current)
(here/"TRANSFER_CHECK.json").write_text(json.dumps(dict(source_export_sha256=hashlib.sha256(export.read_bytes()).hexdigest(),
    member_paths_verified=True, prior_bibliography_preserved=True,included_active_figure_folders=inventory,
    count_direct_graphic_files=sum(map(len,inventory.values()))),indent=2)+"\n",encoding="utf-8")
print("Checked",sum(map(len,inventory.values())),"directly included figures across seven folders.")
