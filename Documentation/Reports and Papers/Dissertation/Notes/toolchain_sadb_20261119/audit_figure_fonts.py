"""Reproduce font, object, text, geometry, and two-renderer checks for all repaired figures.

Run with Python, PyMuPDF, fontTools, Pillow, NumPy, pdffonts, and pdftoppm.
The two Poppler utilities can be selected with command-line options.
"""
import argparse
import hashlib
import json
import math
import re
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory

import fontTools
import numpy as np
import pymupdf as pdf
from PIL import Image, ImageDraw

from embed_figure_fonts import embed_figure


HERE = Path(__file__).resolve().parent
FIGURES = HERE.parents[1] / "Figures"
# PDF graphics transformation matrices read from the integrated November 19 proof.
PRINT_SCALE = {
    "10-introduction/roadmap.pdf": 1.0,
    "15-background/steeleknee.pdf": 0.85768,
    "30-torque/FlxPin_preopt_publication.pdf": 0.65814,
    "30-torque/FlxPin_optimized_publication.pdf": 0.65814,
    "30-torque/ExtPin_complete_publication.pdf": 0.58502,
    "30-torque/cam_current_legend.pdf": 0.84142,
    "30-torque/Ext20mm_physical_publication.pdf": 1.07578,
    "50-database/sadb_bubble.pdf": 0.4968,
    "60-sns/controller_retest_and_optimization.pdf": 1.00435,
    "60-sns/sns_knee_reflex_authoritative.pdf": 1.00006,
    "94-AppendixC/extensor_cam_AF_bones.pdf": 0.95514,
    "94-AppendixC/extensor_cam_GH_bones.pdf": 0.95514,
}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def run(arguments):
    result = subprocess.run(arguments, check=True, capture_output=True)
    return result.stdout.decode(errors="replace"), result.stderr.decode(errors="replace")


def font_resources(document):
    by_xref = {}
    for page in document:
        for xref, _, kind, name, _, encoding, _ in page.get_fonts(full=True):
            data = document.extract_font(xref)[3]
            by_xref[xref] = {"name": name, "kind": kind, "encoding": encoding,
                             "program_sha256": sha(data) if data else None}
    return by_xref


def poppler_render(path, prefix, dpi, executable):
    run([executable, "-singlefile", "-r", str(dpi), "-png", str(path), str(prefix)])
    with Image.open(prefix.with_suffix(".png")) as result:
        return np.array(result.convert("RGB"))


def audit(relative, scale, scratch, pdftoppm, pdffonts):
    current = FIGURES / relative
    original_name = "Ext20mm_physical_20261026_publication.pdf" if current.name == "Ext20mm_physical_publication.pdf" else current.name
    original = current.parent / "Deprecated/20261119_font_embedding" / original_name
    repeated = scratch / current.name
    embed_figure(original, repeated)
    old_data, new_data = original.read_bytes(), current.read_bytes()
    assert new_data.startswith(old_data) and new_data == repeated.read_bytes()
    with pdf.open(original) as before, pdf.open(current) as after:
        assert len(before) == len(after) == 1
        fonts_before, fonts_after = font_resources(before), font_resources(after)
        missing = sorted(xref for xref, font in fonts_before.items() if font["program_sha256"] is None)
        assert len(fonts_after) == len(fonts_before) and all(f["program_sha256"] for f in fonts_after.values())
        assert all(fonts_after[xref] == font for xref, font in fonts_before.items() if xref not in missing)
        changed_objects = [n for n in range(1, before.xref_length()) if before.xref_object(n) != after.xref_object(n)]
        streams = [n for n in range(1, before.xref_length()) if before.xref_is_stream(n)]
        assert changed_objects == missing
        assert all(before.xref_stream_raw(n) == after.xref_stream_raw(n) for n in streams)
        first, last = before[0], after[0]
        assert (first.rect, first.cropbox, first.mediabox, first.rotation) == (
            last.rect, last.cropbox, last.mediabox, last.rotation
        )
        assert first.get_text() == last.get_text() and first.get_text("words") == last.get_text("words")
        raster = {}
        for dpi in (72, 144, 216):
            matrix = pdf.Matrix(dpi / 72, dpi / 72)
            old = first.get_pixmap(matrix=matrix)
            new = last.get_pixmap(matrix=matrix)
            assert (old.width, old.height, old.samples) == (new.width, new.height, new.samples)
            raster[str(dpi)] = {"size_px": [old.width, old.height], "rgb_sha256": sha(old.samples),
                                "different_pixels": 0}
        old_print = poppler_render(original, scratch / "previous", 216 * scale, pdftoppm)
        new_print = poppler_render(current, scratch / "embedded", 216 * scale, pdftoppm)
        diff = np.any(old_print != new_print, axis=2)
        height, width = diff.shape
        glyph_area = Image.new("1", (width, height))
        draw = ImageDraw.Draw(glyph_area)
        font_names = {fonts_before[xref]["name"] for xref in missing}
        spans = [span for block in first.get_text("dict")["blocks"] for line in block.get("lines", [])
                 for span in line["spans"] if span["font"] in font_names]
        for span in spans:
            x0, y0, x1, y1 = span["bbox"]
            draw.rectangle((math.floor(x0 * 3 * scale) - 3, math.floor(y0 * 3 * scale) - 3,
                            math.ceil(x1 * 3 * scale) + 3, math.ceil(y1 * 3 * scale) + 3), fill=1)
        outside = int(np.count_nonzero(diff & ~np.asarray(glyph_area, dtype=bool)))
        assert outside == 0
        font_table = run([pdffonts, str(current)])[0].splitlines()
        # Poppler's table ends in embedding, subset, Unicode, object number, generation.
        assert all(line.split()[-5] == "yes" for line in font_table[2:] if line.strip())
        return {
            "figure": relative, "original_sha256": sha(old_data), "embedded_sha256": sha(new_data),
            "original_bytes": len(old_data), "embedded_bytes": len(new_data),
            "original_file_is_exact_prefix": True, "second_build_sha256_matches": True,
            "original_streams_byte_identical": len(streams), "old_objects_changed": changed_objects,
            "fonts_before_by_xref": fonts_before, "fonts_after_by_xref": fonts_after,
            "boxes_rotation_and_word_positions_match": True, "extracted_words": len(first.get_text("words")),
            "text_utf8_sha256": sha(first.get_text().encode()), "page_width_height_pt": list(first.rect)[2:],
            "mupdf_render_by_native_dpi": raster,
            "poppler_render_at_216_printed_dpi": {
                "manuscript_scale": scale, "size_px": [width, height], "different_pixels": int(diff.sum()),
                "different_pixels_outside_previous_helvetica_text": outside,
                "previous_helvetica_text_spans": len(spans),
                "old_rgb_sha256": sha(old_print.tobytes()), "new_rgb_sha256": sha(new_print.tobytes()),
            },
            "pdffonts": font_table,
        }


def audit_current_figure(relative, expected, pdffonts):
    path = FIGURES / relative
    assert sha(path.read_bytes()) == expected, relative
    with pdf.open(path) as document:
        assert len(document) == 1, relative
        fonts = font_resources(document)
        assert fonts and all(f["program_sha256"] for f in fonts.values()), relative
        bounds = tuple(document[0].rect)
    table = run([pdffonts, str(path)])[0].splitlines()
    assert table[2:] and all(row.split()[-5] == "yes" for row in table[2:] if row.strip()), relative
    return {"figure": relative, "sha256": expected, "page_box": bounds,
            "fonts_by_xref": fonts, "pdffonts": table}


def audit_dissertation(proof, pdffonts, expected_inputs):
    inputs = set()
    for line in proof.with_suffix(".fls").read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("INPUT "):
            match = re.search(r"(?:^|/)(?:Figures|figs)/(.+\.pdf)$", line[6:].replace("\\", "/"), re.I)
            if match:
                inputs.add(match.group(1))
    missing = set(expected_inputs) - inputs
    if missing:
        raise ValueError(f"Current figures absent from the full-dissertation recorder: {sorted(missing)}")
    for name in inputs:
        with pdf.open(FIGURES / name) as source:
            assert all(font["program_sha256"] for font in font_resources(source).values()), name
    with pdf.open(proof) as document:
        resources = font_resources(document)
        missing = [v["name"] for v in resources.values() if not v["program_sha256"]]
        assert not missing, f"Unembedded fonts in combined proof: {missing}"
        pages = len(document)
    text = run([pdffonts, str(proof)])[0].splitlines()
    assert all(row.split()[-5] == "yes" for row in text[2:] if row.strip())
    return {
        "build": proof.parent.name, "pdf_sha256": sha(proof.read_bytes()), "pages": pages,
        "input_pdf_count": len(inputs), "input_pdfs": sorted(inputs),
        "distinct_compiled_font_resources": len(resources), "unembedded_resources": 0,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdftoppm", default="pdftoppm")
    parser.add_argument("--pdffonts", default="pdffonts")
    parser.add_argument("--proof", type=Path, help="Optionally verify a combined proof and its sibling .fls")
    parser.add_argument("--output", type=Path, default=HERE / "ACTIVE_FIGURE_FONT_AUDIT.json")
    args = parser.parse_args()
    active = json.loads((HERE / "ACTIVE_FIGURE_FONTS.json").read_text(encoding="utf-8"))
    repairs = active["retained_incremental_repairs"]
    current = active["current_chapter4_artwork"]
    with TemporaryDirectory() as tmp:
        result = {
            "pymupdf": pdf.__version__, "fontTools": fontTools.__version__,
            "poppler_render": "".join(run([args.pdftoppm, "-v"])).splitlines()[0],
            "poppler_resource_check": "".join(run([args.pdffonts, "-v"])).splitlines()[0],
            "historical_audit": "FIGURE_FONTS.json",
            "retained_incremental_repairs": [audit(name, PRINT_SCALE[name], Path(tmp), args.pdftoppm, args.pdffonts)
                                             for name in repairs],
            "current_chapter4_artwork": [audit_current_figure(name, expected, args.pdffonts)
                                         for name, expected in current.items()],
        }
    if args.proof:
        result["combined_dissertation"] = audit_dissertation(args.proof, args.pdffonts, set(repairs) | set(current))
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
