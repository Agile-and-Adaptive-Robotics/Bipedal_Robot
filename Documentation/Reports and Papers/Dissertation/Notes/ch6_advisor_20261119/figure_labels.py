"""Update labels in three original vector figures and verify that only the label areas change."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import shutil
from zipfile import ZipFile, ZIP_DEFLATED

from fontTools import subset
from fontTools.ttLib import TTFont
import numpy as np
import pymupdf as pdf

here = Path(__file__).resolve().parent
diss = here.parents[1]
figures = diss/"Figures"
staging = diss.parents[2]/"tmp/ch6_advisor_20261119/labels"
staging.mkdir(exist_ok=True, parents=True)
changes = {
    "60-sns/sns_library.pdf": {"Graded synapse":"Nonspiking synapse"},
    "93-AppendixB/animatlab_spiking_conversion.pdf": {
        "(A1) Graded baseline (oscillates)":"(A1) Nonspiking baseline (oscillates)",
        "Graded baseline (oscillates)":"Nonspiking baseline (oscillates)",
        "(C1) Thresholds only, graded synapses (oscillates)":"(C1) Spiking thresholds only (oscillates)",
        "Thresholds only, graded synapses (oscillates)":"Spiking thresholds only (oscillates)"},
    "95-AppendixD/spiking_ground_neuromechanics_short.pdf":{"Graded":"Nonspiking"},
}
descriptions = {
    "60-sns/sns_library.alt.txt":{"Graded and spiking neurons":"Nonspiking and spiking neurons"},
    "93-AppendixB/animatlab_spiking_conversion_alt.txt":{"graded non-spiking baseline":"nonspiking baseline","graded synapses":"voltage-dependent synapses"},
    "95-AppendixD/spiking_ground_neuromechanics_short_alt.txt":{"saved graded and untuned spiking":"saved nonspiking and untuned spiking"},
    "95-AppendixD/animatlab_short_response.alt.txt":{"graded membrane voltages":"membrane voltages of nonspiking neurons"},
}

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def plot_pixels(page):
    pix=page.get_pixmap(matrix=pdf.Matrix(2,2),alpha=False)
    return np.frombuffer(pix.samples,dtype=np.uint8).reshape(pix.height,pix.width,pix.n)

def update_labels(base, out, replacements):
    with pdf.open(base) as doc:
        page=doc[0]
        original=plot_pixels(page).copy()
        selected=[s for block in page.get_text("dict")["blocks"] for line in block.get("lines",[])
                  for s in line["spans"] if s["text"] in replacements]
        assert Counter(s["text"] for s in selected)==Counter(replacements.keys())
        for s in selected:
            rect=pdf.Rect(s["bbox"])
            # Only text is removed. Lines, paths, and images are preserved.
            page.add_redact_annot(rect+(0,1,0,-1),fill=False,cross_out=False)
        page.apply_redactions(images=0,graphics=0,text=0)
        font_source=Path("C:/Windows/Fonts")/("arialbd.ttf" if "Bold" in selected[0]["font"] else "arial.ttf")
        fontfile=staging/(out.stem+".ttf")
        with TTFont(font_source) as font:
            sub=subset.Subsetter()
            sub.populate(text=" ".join(replacements.values()))
            sub.subset(font)
            font.save(fontfile)
        font=pdf.Font(fontfile=str(fontfile))
        boxes=[]
        for s in selected:
            box=pdf.Rect(s["bbox"])
            replacement=replacements[s["text"]]
            width=font.text_length(replacement,fontsize=s["size"])
            x=(box.x0+box.x1-width)/2
            at=(x,s["origin"][1])
            color=pdf.sRGB_to_pdf(s["color"])
            page.insert_text(at,replacement,fontsize=s["size"],fontname="ChapterLabels",
                             fontfile=str(fontfile),color=color)
            boxes.append((box|pdf.Rect(x,box.y0,x+width,box.y1))+(-2,-2,2,2))
        doc.save(out,garbage=4,deflate=True)
    with pdf.open(out) as doc:
        new=plot_pixels(doc[0])
        delta=np.max(np.abs(original.astype(int)-new.astype(int)),axis=2)
        for box in boxes:
            x0,y0,x1,y1=(box*pdf.Matrix(2,2)).round()
            delta[max(0,y0):y1,max(0,x0):x1]=0
        assert not np.any(delta>2), (base, int(delta.max()),np.sum(delta>2))
        layer=doc[0].get_text()
        assert "graded" not in layer.lower()
        for s in replacements.values():
            assert s in layer, (base,s)
    return {"changed_spans":len(selected),"max_outside_label_pixel_delta_255":int(delta.max())}

audit=[]
portable_archive=here/"PRIOR_FIGURE_LABELS.zip"
rename_log=diss/"Notes/figure_names_and_overview_20261119/FIGURE_RENAMES.json"
original_names={(r["chapter"],r["new"]):r["old"] for r in json.loads(rename_log.read_text())}
for name,pairs in {**changes,**descriptions}.items():
    active=figures/name
    archive=active.parent/"Deprecated/20261119_terminology"/original_names.get((active.parent.name,active.name),active.name)
    archive.parent.mkdir(exist_ok=True,parents=True)
    if not archive.exists():
        if portable_archive.exists():
            with ZipFile(portable_archive) as previous:
                archive.write_bytes(previous.read(archive.relative_to(figures).as_posix()))
        else:
            shutil.copyfile(active,archive)
    out=staging/active.name
    detail={}
    if active.suffix==".pdf":
        detail=update_labels(archive,out,pairs)
    else:
        data=archive.read_bytes()
        for old,new in pairs.items():
            assert data.count(old.encode())==1
            data=data.replace(old.encode(),new.encode())
        out.write_bytes(data)
    # The originals and any rerun of the same label pass are the permitted inputs.
    if active.read_bytes()!=out.read_bytes():
        tmp=active.with_suffix(active.suffix+".nov19-tmp")
        shutil.copyfile(out,tmp)
        tmp.replace(active)
    audit.append(dict(file=name,previous=archive.relative_to(figures).as_posix(),
                      prior_sha256=sha(archive),active_sha256=sha(active),labels=pairs,**detail))
(here/"FIGURE_LABELS.json").write_text(json.dumps(audit,indent=2)+"\n",encoding="utf-8")
if not portable_archive.exists():
    with ZipFile(portable_archive,"w",compression=ZIP_DEFLATED) as z:
        for row in audit:
            z.write(figures/row["previous"],row["previous"])
for r in audit:
    print(r["file"],r.get("max_outside_label_pixel_delta_255","accessible description"))
