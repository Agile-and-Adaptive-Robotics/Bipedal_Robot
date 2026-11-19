"""Compare independent builds, publish both PDFs, and check Chapter 6's destinations."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
from zipfile import ZipFile

import pymupdf as pdf

here=Path(__file__).resolve().parent
diss=here.parents[1]
repo=diss.parents[2]
output=diss/"output/pdf"
full=output/"Bolen_Dissertation_Toolchain_SADb_20261119.pdf"
chapter=output/"Bolen_Chapter_6_Toolchain_20261119.pdf"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_build(build):
    log=(build/"main.log").read_text(encoding="utf8",errors="replace")
    bib=(build/"main.blg").read_text(encoding="utf8",errors="replace")
    invalid=re.findall(r"^!.*|.*LaTeX Warning: [^\n]*(?:undefined|multiply defined|Float too large)|.*Package natbib Warning.*|.*Overfull \\[hv]box.*|.*destination with the same identifier.*|.*Rerun to get cross-references right.*",log,re.I|re.M)
    assert not invalid,(build,invalid)
    assert not re.search(r"^Warning--|error message",bib,re.M),build
    fonts=set()
    with pdf.open(build/"main.pdf") as book:
        page_count=len(book)
        for page in book:
            assert tuple(page.mediabox)==tuple(page.cropbox)==(0,0,612,792) and page.rotation==0,page.number
            fonts.update(font[0] for font in page.get_fonts(full=True))
        for xref in fonts:
            assert (book.extract_font(xref)[3] or
                    book.xref_get_key(xref,"Subtype")[1]=="/Type3" and
                    book.xref_get_key(xref,"CharProcs")[0]!="null"),xref
    return {"pages":page_count,"embedded_font_resources":len(fonts),
            "errors":0,"undefined_citations_or_references":0,"overfull_boxes":0,
            "oversized_floats":0,"bibtex_warnings":0,"pdftex_notices":log.count("pdfTeX warning")}


def extract(book, first, following, path):
    """Keep each border, its destination, the folios, and the subsection bookmarks."""
    dests=book.resolve_names()
    with pdf.open() as part:
        part.insert_pdf(book,from_page=first,to_page=following-1,links=False)
        links={"inside_chapter":0,"in_full_dissertation":0,"web":0}
        for page in book.pages(first,following):
            annots=[]
            for link in page.get_links():
                xref=link["xref"]
                obj=book.xref_object(xref)
                assert not re.search(r"\b\d+ \d+ R\b",obj),obj
                ref=part.get_new_xref()
                part.update_object(ref,obj)
                if link["kind"]==pdf.LINK_URI:
                    links["web"]+=1
                else:
                    key=link["nameddest"]
                    dest=dests[key]
                    if first<=dest["page"]<following:
                        x,y=dest["to"]
                        target=part.page_xref(dest["page"]-first)
                        action=f"<</S/GoTo /D [{target} 0 R /XYZ {x} {y} null]>>"
                        links["inside_chapter"]+=1
                    else:
                        action=f"<</S/GoToR /F {pdf.get_pdf_str(full.name)} /D {pdf.get_pdf_str(key)}>>"
                        links["in_full_dissertation"]+=1
                    part.xref_set_key(ref,"A",action)
                annots.append(f"{ref} 0 R")
            if annots:
                part.xref_set_key(part.page_xref(page.number-first),"Annots","["+" ".join(annots)+"]")
        outline=[]
        for level,title,num,dest in book.get_toc(simple=False):
            if first<num<=following:
                location=dests[dest["nameddest"]]
                x,y=location["to"]
                outline.append([level,title,num-first,{"kind":pdf.LINK_GOTO,"page":num-first-1,
                                                        "to":pdf.Point(x,792-y),"zoom":0}])
        part.set_toc(outline)
        part.set_page_labels([dict(startpage=0,style="D",prefix="",firstpagenum=int(book[first].get_label()))])
        part.set_metadata(dict(author="Ben Bolen",title="Chapter 6. Tools for Constructing and Testing Neuromechanical Models"))
        part.save(path,garbage=4,deflate=True)
    with pdf.open(path) as part:
        assert len(part)==following-first and len(part.get_toc())==len(outline)
        for i,p in enumerate(part):
            src=book[first+i]
            assert p.get_label()==src.get_label() and p.get_text()==src.get_text()
            assert p.get_pixmap().digest==src.get_pixmap().digest,p.number
            original,clipped=src.get_links(),p.get_links()
            assert len(clipped)==len(original),(i,len(original),len(clipped))
            for old,new in zip(original,clipped):
                for key in ("Rect","C","Border","H"):
                    assert book.xref_get_key(old["xref"],key)==part.xref_get_key(new["xref"],key),(i,key)
    return links


def publish(source,target):
    tmp=target.with_suffix(".pdf.part")
    shutil.copyfile(source,tmp)
    tmp.replace(target)
    assert sha(source)==sha(target)


def printed_workflow_sizes(build):
    records=json.loads((here/"FIGURE_PRINT_QA.json").read_text())["workflows"]
    log=(build/"main.log").read_text(encoding="utf8",errors="replace")
    sizes={}
    for r in records:
        filename=r["graphic"]
        tail=log.split("/"+filename+">",1)[1]
        placed_tex_pt=float(re.search(r"Requested size: ([\d.]+)pt",tail)[1])
        original=next((diss/"Figures").glob("*/"+filename))
        with pdf.open(original) as artwork:
            printed_scale=(placed_tex_pt*72/72.27)/artwork[0].rect.width
        relative=printed_scale/r["inclusion_scale"]
        new_head=[round(relative*pt,2) for pt in r["after_head_length_width_print_pt"]]
        text_min=round(relative*r["after_minimum_type_pt"],2)
        assert text_min>=10,(filename,text_min)
        sizes[r["reference_figure"]]={"smallest_workflow_type_pdf_pt":text_min,"arrowhead_length_and_width_pdf_pt":new_head}
    return sizes


if __name__=="__main__":
    args=argparse.ArgumentParser()
    args.add_argument("local",type=Path)
    args.add_argument("transfer",type=Path)
    args=args.parse_args()
    local,other=args.local.resolve(),args.transfer.resolve()
    evidence={"local_build":check_build(local),"transfer_build":check_build(other)}
    evidence["printed_workflow_sizes"]=printed_workflow_sizes(local)
    with pdf.open(local/"main.pdf") as a, pdf.open(other/"main.pdf") as b:
        assert len(a)==len(b) and a.get_toc()==b.get_toc()
        for ap,bp in zip(a,b):
            assert ap.get_text()==bp.get_text(),ap.number
            assert ap.get_pixmap().digest==bp.get_pixmap().digest,ap.number
        evidence["full_page_text_and_72dpi_RGB_renders_matched"]=len(a)
        toc=a.get_toc()
        first=next(num-1 for lv,t,num in toc if lv==1 and "Tools for Constructing" in t)
        following=next(num-1 for lv,t,num in toc if lv==1 and num-1>first)
        evidence["chapter_6"]={"printed_pages":[a[first].get_label(),a[following-1].get_label()],
                               "page_count":following-first,"links":extract(a,first,following,local/"chapter6_checked.pdf")}
        intro=next(num-1 for lv,t,num in toc if lv==1 and t=="Introduction")
        assert a[2].get_text().splitlines()[:2]==["i","Abstract"] and a[intro].get_label()=="1"
        assert a[1].get_text().strip()=="© 2026 Ben Bolen"
        evidence["front_openings"]={t:a[num-1].get_label() for lv,t,num in toc if lv==1 and num<=intro+1}
        figure_labels=["fig:sadb_clusters","fig:sns_reflex_current","fig:sns_syn6_overlay","fig:bpa_exploded_wall","fig:sensorium_hierarchy","fig:discussion_neuromechanical_workflow"]
        aux="".join(p.read_text(errors="replace") for p in (local/"chapters").glob("*.aux"))
        evidence["figure_number_and_folio"]={}
        for lab in figure_labels:
            number,folio=re.search(r"\\newlabel\{"+lab+r"\}\{\{([^}]+)\}\{(\d+)\}",aux).groups()
            evidence["figure_number_and_folio"][lab]=[number,folio]
        assert evidence["figure_number_and_folio"]["fig:discussion_neuromechanical_workflow"][0]=="7.5"
        evidence["visible_link_boxes_in_full_pdf"]=sum(a.xref_get_key(link["xref"],"Border")[1].replace(" ","")=="[001]" for p in a for link in p.get_links())
        assert evidence["visible_link_boxes_in_full_pdf"]>1000
    transfer=json.loads((here/"TRANSFER.json").read_text())
    with ZipFile(diss/"Bolen_Dissertation.zip") as z:
        evidence["chapter_sources_matching_author_export"]=[]
        for file in sorted((diss/"Final/chapters").glob("*.tex")):
            key="chapters/"+file.name
            if key not in z.namelist() or key in transfer["members"]:
                continue
            assert file.read_bytes().replace(b"\r\n",b"\n")==z.read(key).replace(b"\r\n",b"\n"),key
            evidence["chapter_sources_matching_author_export"].append(key)
    output.mkdir(parents=True,exist_ok=True)
    publish(local/"main.pdf",full)
    publish(local/"chapter6_checked.pdf",chapter)
    evidence["sha256"]={x.relative_to(diss).as_posix():sha(x) for x in (full,chapter,here/"Dissertation_Overleaf_updates.zip")}
    evidence["build_directories"]=[p.relative_to(repo).as_posix() for p in (local,other)]
    (here/"FINAL_VERIFICATION.json").write_text(json.dumps(evidence,indent=2)+"\n",encoding="utf8")
    print(json.dumps(evidence,indent=2))
