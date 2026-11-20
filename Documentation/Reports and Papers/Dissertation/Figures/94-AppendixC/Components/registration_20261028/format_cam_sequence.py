"""Lay out the eight validated cam poses for the appendix and paper supplement.

The source PDFs are the vector output of make_cam_label_figures.py, checked
against the eight poses in CAM_LABEL_VERIFICATION.json. Reuse their plotted
geometry, point leaders, symbols, and exact legend strokes at original scale.
"""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import re
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile

import pymupdf as pdf


APPENDIX = Path(__file__).resolve().parents[2]
DISSERTATION = APPENDIX.parents[1]
NOTES = DISSERTATION / "Notes"
sys.path.insert(0, str(NOTES / "toolchain_sadb_20261119"))
from embed_figure_fonts import embed_figure

WIDTH, OFFSET = 512, 21
KEY_TEXT = ("Femur cylinder", "Femur wall", "Condyle envelope", "Tibia lower", "Tibia upper",
            "Bend radius", "Optimized route", "Cam contacts", "Origin p1", "Insertion p9", "Reference bones")
FONTS = Path(os.environ.get("WINDIR", "C:/Windows")) / "Fonts"


def spans(page):
    return [span for block in page.get_text("dict")["blocks"] if "lines" in block
            for line in block["lines"] for span in line["spans"]]


def read_key(source):
    """Copy each published swatch into its own tightly bounded vector page."""
    swatches = pdf.open()
    refs = source[0].get_contents()
    original_key = []
    for index, ref in enumerate(refs):
        text = source.xref_stream(ref)
        label = re.search(rb"<([0-9a-fA-F]+)>\]TJ", text)
        if not label:
            continue
        name = bytes.fromhex(label[1].decode()).decode("latin1")
        if name not in KEY_TEXT:
            continue
        assert name == KEY_TEXT[len(original_key) // 2]
        baseline = float(re.search(rb"348 ([\d.]+) Tm", text)[1])
        sample = source.xref_stream(refs[index-1])
        assert b"BT" not in sample and sample.strip().startswith(b"q")
        page = swatches.new_page(width=24, height=12)
        obj = swatches.get_new_xref()
        swatches.update_object(obj, "<<>>")
        swatches.update_stream(obj, f"q 1 0 0 1 -319 {6-baseline-3.5:g} cm\n".encode() + sample + b"\nQ")
        page.set_contents(obj)
        original_key.extend([refs[index-1], ref])
    assert len(swatches) == len(KEY_TEXT)
    return swatches, original_key


def key(page, swatches, columns, box):
    page.draw_rect(box, color=(.1, .1, .1), width=.85)
    width = pdf.Font(fontfile=str(FONTS / "arialbd.ttf")).text_length("Legend", fontsize=12)
    page.insert_text(((box.x0+box.x1-width)/2, box.y0+15), "Legend", fontname="ArialBold", fontsize=12)
    rows = (len(KEY_TEXT) + columns - 1) // columns
    col_spacing = (box.width-16)/columns
    row_spacing = 18 if columns == 1 else 15
    for index, name in enumerate(KEY_TEXT):
        col, row = divmod(index, rows)
        x = box.x0 + 8 + col*col_spacing
        y = box.y0 + 29 + row*row_spacing
        page.show_pdf_page(pdf.Rect(x,y-6,x+24,y+6), swatches, index)
        page.insert_text((x+28, y+3.7), name, fontname="Arial", fontsize=11)


def arrange(source, part, swatches, old_key):
    page = source[0]
    if part == "GH":
        for obj in old_key:
            source.update_stream(obj, b"")
    letters = {"(A)", "(D)"} if part == "AF" else {"(G)"}
    before = spans(page)
    moved = [entry for entry in before if entry["text"] in letters]
    assert {entry["text"] for entry in moved} == letters
    for entry in moved:
        box = pdf.Rect(entry["bbox"])
        mid_y = (box.y0+box.y1)/2
        page.add_redact_annot((box.x0, mid_y-1, box.x1, mid_y+1), fill=False, cross_out=False)
    page.apply_redactions(images=pdf.PDF_REDACT_IMAGE_NONE, graphics=pdf.PDF_REDACT_LINE_ART_NONE)
    assert Counter(x["text"] for x in spans(page)) == Counter(x["text"] for x in before if x["text"] not in letters)

    result = pdf.open()
    dst = result.new_page(width=WIDTH, height=561 if part == "AF" else 360)
    dst.show_pdf_page(pdf.Rect(OFFSET, 0, OFFSET+page.rect.width, page.rect.height), source, 0)
    dst.insert_font(fontname="Arial", fontfile=str(FONTS / "arial.ttf"))
    dst.insert_font(fontname="ArialBold", fontfile=str(FONTS / "arialbd.ttf"))
    for entry in moved:
        dst.insert_text((1.5, entry["origin"][1]), entry["text"], fontname="ArialBold", fontsize=entry["size"])
    if part == "AF":
        key(dst, swatches, 4, pdf.Rect(4,483,508,557))
    else:
        key(dst, swatches, 1, pdf.Rect(340,5,496,228))
    result.subset_fonts()
    return result


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--source", type=Path, help="Override the two registered-pose inputs in the source archive.")
    cli.add_argument("--output", type=Path, default=APPENDIX)
    args = cli.parse_args()
    checked = json.loads((NOTES/"revision_20261117/CAM_LABEL_VERIFICATION.json").read_text())
    args.output.mkdir(parents=True, exist_ok=True)
    archive_path = NOTES/"toolchain_sadb_20261119/FIGURE_FONT_SOURCES.zip"
    with TemporaryDirectory() as scratch, ZipFile(archive_path) as archive:
        sources = {}
        for part in ("AF","GH"):
            name = f"extensor_cam_{part}_bones.pdf"
            if args.source:
                src = args.source/name
            else:
                member = "Documentation/Reports and Papers/Dissertation/Figures/94-AppendixC/Deprecated/20261119_font_embedding/" + name
                src = Path(scratch)/("original_" + name)
                src.write_bytes(archive.read(member))
            assert hashlib.sha256(src.read_bytes()).hexdigest() == checked[part]["revised_sha256"]
            embedded = Path(scratch)/name
            embed_figure(src, embedded)
            sources[part] = pdf.open(embedded)
        swatches, old_key = read_key(sources["GH"])
        for part, source in sources.items():
            result = arrange(source, part, swatches, old_key)
            name = f"extensor_cam_{part}_bones.pdf"
            result.save(args.output/name, garbage=4, deflate=True)
            result.close()
            source.close()
        swatches.close()


if __name__ == "__main__":
    main()
