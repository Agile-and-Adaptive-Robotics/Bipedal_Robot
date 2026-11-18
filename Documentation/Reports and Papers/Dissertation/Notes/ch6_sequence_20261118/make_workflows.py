"""Compose the print figures from the four native PowerPoint diagrams.

export_workflows.ps1 creates the full slide PDF from the verified thumb-drive deck.
The four selected pages can be retained as a small, standalone source for future builds.
"""
import argparse
from pathlib import Path

import pymupdf as pdf


SLIDES = (64, 67, 70, 71)
WIDTH = 466
DARK = (0.09, 0.20, 0.23)


def native_piece(target, source, index, bounds, x, y, width):
    """Crop a native diagram, retaining the source lines, text, and arrows as vectors."""
    clip = pdf.Rect(bounds)
    small = pdf.open()
    small.insert_pdf(source, from_page=index, to_page=index)
    small[0].clip_to_rect(clip)
    rect = pdf.Rect(x,y,x+width,y+width*clip.height/clip.width)
    target.show_pdf_page(rect, small, 0, clip=clip)
    small.close()
    return rect


def label(page, x, y, string, size=11, bold=False):
    font = "hebo" if bold else "helv"
    page.insert_text((x,y),string,fontname=font,fontsize=size,color=DARK)


def save(document, path, title):
    document.set_metadata({"title":title,"author":"Benjamin Bolen",
                           "subject":"Adapted from Bolen's PhD defense diagrams"})
    document.subset_fonts()
    document.save(path,garbage=4,deflate=True,clean=True,no_new_id=True)
    document.close()


def descriptions(folder):
    texts = {
        "simulation_workflow_authoritative":
            "Two rows from the defense share the same three-column layout. In the upper row, OpenSim leads to Ben's model repairs and actuator interfaces and then to MuJoCo and SNS. In the lower row, SOLIDWORKS leads to Ben's CAD import and SNS-Simscape library and then to MuJoCo through the MathWorks bridge. The middle column identifies Ben's additions. An executable mechanical model is prepared before runtime signal exchange.",
        "controller_retest_and_optimization":
            "Panel A shows two left-leg paths from the saved W2L circuit, an extensor pattern neuron driving the extensor motor population and a flexor pattern neuron driving the flexor motor population. Beside them, the CAD loop runs from revising the geometry to updating the SNS-Simscape model and interfaces to retesting in MuJoCo through the MathWorks bridge. Panel B follows a separate sequence: choose a function and its target, run and score trials within parameter bounds, then rerun the selected settings.",
    }
    for stem, text in texts.items():
        (folder/(stem+".alt.txt")).write_text(text+"\n",encoding="utf-8")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("slides_pdf",type=Path)
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    source=pdf.open(args.slides_pdf)
    if len(source)==4:
        pages=source
    else:
        pages=pdf.open()
        for n in SLIDES:
            pages.insert_pdf(source,from_page=n-1,to_page=n-1)
        source_path=args.output/"workflow_slides_64_67_70_71.pdf"
        pages.save(source_path,garbage=4,deflate=True,no_new_id=True)
    for n, page, keyword in zip(SLIDES,pages,("OpenSim","SOLIDWORKS","Retest in MuJoCo","Choose the function")):
        if keyword not in page.get_text():
            raise ValueError(f"The source PDF page for slide {n} does not contain {keyword!r}")

    document=pdf.open()
    page=document.new_page(width=WIDTH,height=151)
    native_piece(page,pages,0,(20,134,691,238),0,0,WIDTH)
    native_piece(page,pages,1,(14,149,685,246),0,83,WIDTH)
    save(document,args.output/"simulation_workflow_authoritative.pdf",
         "OpenSim and CAD routes to a coupled neuromechanical model")

    document=pdf.open()
    page=document.new_page(width=WIDTH,height=302)
    label(page,2,15,"A   Revising geometry and retesting the controller",13,True)
    label(page,4,45,"Pattern neuron",11.5,True)
    label(page,132,45,"Motor neuron",11.5,True)
    # The defense uses one raster image for its two circuit rows. Keep the actual
    # exported wiring and glyphs. Set the smaller slide labels anew at print size.
    for name,top,bounds in [
        ("Left knee extension",54,(82,158,330,196)),
        ("Left knee flexion",120,(82,224,330,262)),
    ]:
        clip=pdf.Rect(bounds)
        raster=pages[2].get_pixmap(matrix=pdf.Matrix(4,4),clip=clip,alpha=False)
        target=pdf.Rect(3,top,209,top+206*clip.height/clip.width)
        page.insert_image(target,stream=raster.tobytes("png"))
        label(page,4,target.y1+15,name,11)
    native_piece(page,pages,2,(416,94,692,321),235,30,230)
    label(page,2,237,"B   Staged controller optimization",13,True)
    native_piece(page,pages,3,(29,122,689,206),0,243,WIDTH)
    save(document,args.output/"controller_retest_and_optimization.pdf",
         "Controller retesting and staged parameter optimization")
    descriptions(args.output)
    for path in sorted(args.output.glob("*.pdf")):
        if path.name.startswith("workflow_slides"):
            continue
        with pdf.open(path) as fig:
            spans=[span for b in fig[0].get_text("dict")["blocks"] if "lines" in b
                   for line in b["lines"] for span in line["spans"] if span["text"].strip()]
            if min(s["size"] for s in spans)<10:
                raise ValueError(f"Text smaller than 10 pt in {path}: {min(s['size'] for s in spans)}")
            fig[0].get_pixmap(dpi=144,alpha=False).save(args.output/(path.stem+".png"))
            print(path.name,min(s["size"] for s in spans),len(spans))
    pages.close()
    if pages is not source:
        source.close()


if __name__=="__main__":
    main()
