"""Rebuild the vector callouts over the preserved SADb screenshot PDF.

Run with PyMuPDF installed. The default output is the active chapter figure;
--output permits a separate proof. No manuscript or application files are edited.
"""

from pathlib import Path
import argparse
import hashlib

import pymupdf


HERE = Path(__file__).resolve().parent
BASE = HERE / "sadb_bubble_cropped_base.pdf"
BASE_SHA256 = "0ae64beec42faec19da753452e607c01bd300df77f43e84c3f7cf4663ca7349f"
CALLOUTS = [
    ("A", "Tabs", (5, 72, 148, 126), [(140, 103), (163, 103)]),
    ("B", "Map controls", (560, 4, 900, 44), [(680, 39), (680, 133)]),
    ("C", "Focus record", (5, 175, 148, 240), [(140, 215), (163, 215)]),
    ("D", "Curation and\ntranslation", (5, 380, 148, 455), [(140, 419), (163, 419)]),
    ("E", "Citation neighborhood", (400, 590, 740, 632), [(560, 589), (560, 557), (665, 305)]),
    ("F", "Connected\npapers", (958, 175, 1096, 242), [(958, 215), (936, 215)]),
    ("G", "Topic color key", (5, 590, 370, 632), [(130, 589), (130, 577), (430, 559)]),
]


def build(output):
    assert hashlib.sha256(BASE.read_bytes()).hexdigest() == BASE_SHA256
    source = pymupdf.open(BASE)
    assert len(source) == 1 and tuple(source[0].rect) == (0, 0, 792, 496)
    result = pymupdf.open()
    page = result.new_page(width=1100, height=640)
    page.show_pdf_page(pymupdf.Rect(154, 70, 946, 566), source, 0)
    blue, white = (0.05, 0.24, 0.43), (1, 1, 1)
    for letter, title, rect, points in CALLOUTS:
        page.draw_rect(pymupdf.Rect(rect), color=None, fill=white, overlay=True)
        page.insert_text((rect[0] + 3, rect[1] + 18), letter,
                         fontsize=17, fontname="hebo", color=blue)
        spare = page.insert_textbox(
            pymupdf.Rect(rect[0] + 23, rect[1] + 3, rect[2] - 3, rect[3]),
            title, fontsize=18, fontname="helv", color=(0, 0, 0))
        assert spare >= 0, (letter, spare)
        leader = page.new_shape()
        leader.draw_polyline([pymupdf.Point(*point) for point in points])
        leader.finish(color=blue, width=1.6, closePath=False)
        leader.commit()
        page.draw_circle(points[-1], 3, color=blue, fill=white, width=1.4)
    output = Path(output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    result.save(output, garbage=4, deflate=True)
    result.close()
    source.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=HERE.parent.parent / "sadb_bubble.pdf")
    build(parser.parse_args().output)
