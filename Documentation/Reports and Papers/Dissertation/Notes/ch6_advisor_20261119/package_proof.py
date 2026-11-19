"""Verify the complete proof, extract Chapter 6, and package the four edited live sources."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
from zipfile import ZipFile, ZipInfo, ZIP_DEFLATED

import pymupdf as fitz
from PIL import Image, ImageDraw

here = Path(__file__).resolve().parent
diss = here.parents[1]
root = diss.parents[2]
parser = argparse.ArgumentParser()
parser.add_argument("build", type=Path)
build = parser.parse_args().build.resolve()
output = diss / "output/pdf"
output.mkdir(exist_ok=True,parents=True)
render = build.parent / "visual"
render.mkdir(exist_ok=True)

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
def png(page, scale):
    pix = page.get_pixmap(matrix=fitz.Matrix(scale,scale), alpha=False)
    return Image.frombytes("RGB", [pix.width,pix.height],pix.samples)
def safe_copy(src, dest):
    tmp = dest.with_suffix(".part")
    shutil.copyfile(src, tmp)
    tmp.replace(dest)

full = output/"Bolen_Dissertation_20261119.pdf"
chapter = output/"Bolen_Dissertation_Chapter_6.pdf"
book = fitz.open(build / "main.pdf")
toc = book.get_toc()
first = next(p for level,title,p in toc if level==1 and "Tools for Constructing" in title)
following = next(p for level,title,p in toc if level==1 and p>first)
rows = [entry for entry in toc if first <= entry[2] < following]
clip = fitz.open()
clip.insert_pdf(book,from_page=first-1,to_page=following-2)
clip.set_toc([[level,title,p-first+1] for level,title,p in rows])
# The extract retains the dissertation's printed folios and matching viewer labels.
chapter_start = int(re.search(r"\}\{(\d+)\}\{chapter\.6", (build/"main.toc").read_text()).group(1))
clip.set_page_labels([dict(startpage=0, style="D", prefix="", firstpagenum=chapter_start)])
clip.set_metadata(dict(title="Chapter 6. Tools for Constructing and Testing Neuromechanical Models",
                       author=book.metadata.get("author") or "Ben Bolen"))
tmp = build/"chapter6.pdf"
clip.save(tmp,garbage=4,deflate=True)
clip.close()
with fitz.open(tmp) as extracted:
    assert len(extracted) == following-first and len(extracted.get_toc()) == len(rows)
    assert all(extracted[i].get_text() == book[first+i-1].get_text() for i in range(len(extracted)))
    assert extracted[0].get_label()==str(chapter_start)
safe_copy(tmp,chapter)
safe_copy(build/"main.pdf",full)
assert digest(full)==digest(build/"main.pdf")

names = ["60-sns.tex","70-overall.tex","93-AppendixB.tex","95-AppendixD.tex"]
label_changes = json.loads((here/"FIGURE_LABELS.json").read_text(encoding="utf-8"))
new_figures = {"60-sns":("syn6_suspended_kinematics","fig:sns_syn6_overlay"),
               "95-AppendixD":("research_toolchain_workflow","fig:app_research_toolchain")}
for folder,(stem,_) in new_figures.items():
    for ending in (".pdf",".alt.txt"):
        f=folder+"/"+stem+ending
        label_changes.append(dict(file=f,active_sha256=digest(diss/"Figures"/f)))
package = here/"Revised_chapters_and_appendices.zip"
with ZipFile(package,"w",compression=ZIP_DEFLATED) as z:
    z.write(here/"Ong_2019.bib","Ong_2019.bib")
    info = ZipInfo("Transfer.txt", (2026, 11, 19, 0, 0, 0))
    z.writestr(info,"Copy the four chapter sources into chapters/ and the figs/ files into their corresponding Overleaf figure folders. Append the entry in Ong_2019.bib to the existing thesis.bib.\n",compress_type=ZIP_DEFLATED)
    for name in names:
        z.write(diss/"Final/chapters"/name,"chapters/"+name)
    for change in label_changes:
        z.write(diss/"Figures"/change["file"],"figs/"+change["file"])
with ZipFile(package) as z:
    assert z.read("Ong_2019.bib")== (here/"Ong_2019.bib").read_bytes()
    assert (diss/"Final/thesis.bib").read_bytes().count(b"@article{ong_predicting_2019,")==1
    for name in names:
        assert z.read("chapters/"+name)==(diss/"Final/chapters"/name).read_bytes()
    for change in label_changes:
        assert hashlib.sha256(z.read("figs/"+change["file"])).hexdigest()==change["active_sha256"]
with ZipFile(here/"LIVE_SOURCES.zip") as z:
    assert (diss/"Final/chapters/10-introduction.tex").read_text(encoding="utf-8")==z.read("10-introduction.tex").decode()
    assert (diss/"Final/chapters/30-torque.tex").read_text(encoding="utf-8")==z.read("30-torque.tex").decode()
    # Audit every directly included image and label for the changed chapters against the live source.
    assets = {}
    for name in names:
        live=z.read(name).decode()
        current=(diss/"Final/chapters"/name).read_text(encoding="utf-8")
        for pat in (r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}",r"\\label\{([^}]+)\}"):
            expected = re.findall(pat,live)
            if name.removesuffix(".tex") in new_figures:
                stem, key = new_figures[name.removesuffix(".tex")]
                expected.append(key if "label" in pat else stem+".pdf")
            assert sorted(expected) == sorted(re.findall(pat, current))
        figs = re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}",current)
        base = diss/"Figures"/name.removesuffix(".tex")
        for f in figs:
            assert "Components" not in f and "Deprecated" not in f and (base/f).is_file(), (name,f)
        assets[name]=figs

full_log=(build/"main.log").read_text(errors="replace")
bib=(build/"main.blg").read_text(errors="replace")
assert not re.search(r"^!|LaTeX Warning: (?:.*undefined|.*multiply defined)|Package natbib Warning|Overfull \\[hv]box",
                     full_log,re.M|re.I)
assert not re.search(r"^Warning--|error message",bib,re.M)
prefix = chapter_start-first
pdf6 = [book[i-1].get_text() for i in range(first,following)]
aux_text="\n".join(path.read_text(encoding="utf-8",errors="replace") for path in (build/"chapters").glob("*.aux"))
def physical(label):
    match=re.search(r"\\newlabel\{"+re.escape(label)+r"\}\{\{[^}]*\}\{(\d+)\}",aux_text)
    assert match,label
    return int(match[1])-prefix
special={key:physical(label) for key,label in {"knee":"fig:sns_reflex_current","arm":"fig:sns_beverage_current",
    "library":"fig:simulink_sns_implementation","model_table":"tab:sns_model_configurations","force_table":"tab:force-r2",
    "future":"sec:fw_pipeline","conversion_plot":"fig:spike_animatlab","spike_plot":"fig:spike_rhythm",
    "syn6_overlay":"fig:sns_syn6_overlay","research_workflow":"fig:app_research_toolchain"}.items()}
print("Chapter",first,following-1,"printed",chapter_start,following-1+prefix,"special",special)
# A complete visual overview and individual page rasters for proofing at publication size.
for begin in range(first,following,6):
    pages=list(range(begin,min(begin+6,following)))
    sheet=Image.new("RGB",(3*360,2*483),"#dddddd")
    pen=ImageDraw.Draw(sheet)
    for j,p in enumerate(pages):
        im=png(book[p-1],.57)
        x=(j%3)*360; y=(j//3)*483
        sheet.paste(im,(x+5,y+17))
        pen.text((x+6,y+2),f"p. {p+prefix}",fill="black")
    sheet.save(render/f"chapter_{begin+prefix:03}.png")
assert all(k in special for k in ("knee","arm","model_table","library","force_table","future"))
key_pages = set([first,following-1,*special.values()])
for p in key_pages:
    png(book[p-1],1.3).save(render/f"page_{p+prefix:03}.png")
terminology = {}
for i in range(first,len(book)+1):
    value=book[i-1].get_text()
    if re.search(r"\bgraded\b",value,re.I):
        terminology[str(i+prefix)]=[line.strip()[:140] for line in value.splitlines() if re.search(r"\bgraded\b",line,re.I)]
results = dict(full_pages=len(book), chapter_pages=following-first, physical_start=first, physical_end=following-1,
    printed_start=chapter_start, printed_end=following-1+prefix, captured_figure_pages={k:p+prefix for k,p in special.items()},
    chapter_bookmarks=rows, sha256={p.relative_to(root).as_posix():digest(p) for p in (full,chapter,package)},
    figures=assets, source_sha256={name:digest(diss/"Final/chapters"/name) for name in [*names,"10-introduction.tex","30-torque.tex"]},
    wording_in_embedded_art_or_other_chapters=terminology, TeX_overfull_boxes=0, undefined_references=0, undefined_citations=0,
    BibTeX_warnings=0, pdfTeX_warnings=len(re.findall(r"pdfTeX warning",full_log)),
    capture_and_Introduction_match=True, capture_and_torque_chapter_match=True)
(here/"VERIFICATION.json").write_text(json.dumps(results,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
print("Terminology matches:",terminology)
