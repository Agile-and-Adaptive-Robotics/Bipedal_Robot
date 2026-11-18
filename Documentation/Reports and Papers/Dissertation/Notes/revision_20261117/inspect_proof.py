"""Render all pages in scope at 144 dpi, create page overviews, and measure placed figure fonts."""
import argparse
import json
import re
from pathlib import Path
import pymupdf
from PIL import Image, ImageDraw

p = argparse.ArgumentParser()
p.add_argument("build", type=Path)
p.add_argument("--dpi", type=int, default=144)
args=p.parse_args()
build=args.build.resolve()
renders=build/"review"
renders.mkdir(exist_ok=True)
pdf=pymupdf.open(build/"main.pdf")
toc=(build/"main.toc").read_text(encoding="utf-8",errors="replace")
starts={name:int(value) for name,value in re.findall(
    r"\\contentsline \{(?:chapter|app)\}\{\\numberline \{([0-9A-E]+)\}.*?\}\{(\d+)\}",toc)}
offset=next(i for i,page in enumerate(pdf) if page.get_label()=="1" and i>0) if any(p.get_label()=="1" and i>0 for i,p in enumerate(pdf)) else None
if offset is None or "Introduction" not in pdf[offset].get_text():
    offset = next(i for i,p in enumerate(pdf) if re.search(r"Chapter\s*1\s*Introduction",p.get_text(),re.I))
def print_page_number(index):
    return index-offset+1
def span_list(page):
    return [s for block in page.get_text("dict")["blocks"] if "lines" in block
            for line in block["lines"] for s in line["spans"] if s["text"].strip()]
records=[]
groups={}
keys=["6","7","A","B","C","D","E"]
ends={"6":starts["7"],"7":starts["8"],"A":starts["B"],"B":starts["C"],
      "C":starts["D"],"D":starts["E"],"E":len(pdf)-offset+2}
for group in keys:
    images=[]
    for i in range(starts[group]+offset-1,min(len(pdf),ends[group]+offset-1)):
        pg=pdf[i]
        spans=span_list(pg)
        pix=pg.get_pixmap(dpi=args.dpi,alpha=False)
        dest=renders/f"p{i+1:03}.png"
        pix.save(dest)
        images.append(dest)
        counts={}
        for s in spans:
            counts[round(s["size"],1)]=counts.get(round(s["size"],1),0)+1
        records.append({"group":group,"pdf_page":i+1,"printed_page":print_page_number(i),
                        "font_size_counts":counts,"span_count":len(spans),
                        "small_spans":[[round(s["size"],1),s["text"]] for s in spans if s["size"]<9.8][:65],
                        "min_x":round(min(s["bbox"][0] for s in spans),1) if spans else None,
                        "max_x":round(max(s["bbox"][2] for s in spans),1) if spans else None})
    groups[group]=[images[0].name,images[-1].name]
    for j in range(0,len(images),6):
        canvas=Image.new("RGB",(1140,1020),"#e5e6e7")
        draw=ImageDraw.Draw(canvas)
        for k,path in enumerate(images[j:j+6]):
            with Image.open(path) as img:
                img.thumbnail((368,470))
                x=(k%3)*380+(380-img.width)//2
                y=k//3*510+11
                canvas.paste(img,(x,y))
                draw.text(((k%3)*380+12,y+470),f"PDF {int(path.stem[1:])}  |  printed {int(path.stem[1:])-offset}",fill="#222222")
        canvas.save(renders/f"group_{group}_{j//6+1:02}.jpg",quality=88)
data={"pdf_pages":len(pdf),"first_arabic_pdf_page":offset+1,"dpi":args.dpi,"ranges":groups,"pages":records}
(renders/"layout.json").write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8")
print(groups,flush=True)
