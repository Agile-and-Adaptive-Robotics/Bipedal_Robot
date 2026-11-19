"""Measure rendered PDF ink, including visible hyperlink borders, on every page.

Requires Python 3.10+, NumPy, PyMuPDF, and Poppler's pdftoppm. Pass a PDF or a
build directory containing main.pdf. --pdftoppm selects a particular renderer.
"""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

import numpy as np
import pymupdf


def poppler_details(requested):
    program = requested or os.environ.get("PDFTOPPM") or shutil.which("pdftoppm")
    if not program:
        raise SystemExit("Supply --pdftoppm or put Poppler's pdftoppm on PATH.")
    output = subprocess.run([str(program), "-v"], capture_output=True, check=True)
    version = (output.stdout + output.stderr).decode(errors="replace")
    if "Poppler" not in version:
        raise SystemExit(f"The specified renderer is not Poppler: {program}")
    return str(program), version.splitlines()[0]


def read_pgm(stream):
    def field():
        line = stream.readline()
        while line.startswith(b"#"):
            line = stream.readline()
        return line.strip()

    if field() != b"P5":
        raise ValueError("Expected a binary grayscale PGM page from Poppler.")
    width, height = map(int, field().split())
    if field() != b"255":
        raise ValueError("Expected 8-bit grayscale samples.")
    data = stream.read(width * height)
    if len(data) != width * height:
        raise ValueError("Poppler output ended inside a page.")
    return np.frombuffer(data, dtype=np.uint8).reshape(height, width)


def inspect_render(pdf, labels, program, dpi):
    sides = ("left", "right", "top", "bottom")
    minima = [None] * 4
    limiting = [[] for _ in sides]
    violations = []
    empty_pages = []
    with tempfile.TemporaryFile() as stderr:
        # Default annotation rendering is intentional. Do not add -hide-annotations.
        args = [program, "-gray", "-r", str(dpi), str(pdf)]
        with subprocess.Popen(args, stdout=subprocess.PIPE, stderr=stderr) as proc:
            try:
                for page_number, label in enumerate(labels, 1):
                    pixels = read_pgm(proc.stdout)
                    height, width = pixels.shape
                    cols = np.flatnonzero(pixels.min(axis=0) < 255)
                    if not len(cols):
                        empty_pages.append(page_number)
                        continue
                    rows = np.flatnonzero(pixels.min(axis=1) < 255)
                    distances = [int(cols[0]), int(width-1-cols[-1]),
                                 int(rows[0]), int(height-1-rows[-1])]
                    page = {"pdf_page": page_number, "printed_page": label}
                    for i, distance in enumerate(distances):
                        if minima[i] is None or distance < minima[i]:
                            minima[i] = distance
                            limiting[i].clear()
                        if distance == minima[i]:
                            limiting[i].append(page)
                    if min(distances) < dpi:
                        violations.append(dict(page, margins_pt=dict(zip(
                            sides, [round(n * 72 / dpi, 4) for n in distances]))))
                if proc.stdout.read(1):
                    raise ValueError("Poppler rendered more pages than the PDF inventory.")
            except BaseException:
                proc.kill()
                raise
        stderr.seek(0)
        messages = stderr.read().decode(errors="replace").splitlines()
        if proc.returncode:
            raise RuntimeError("\n".join(messages))
    return {
        "dpi": dpi,
        "nonwhite_pixel_threshold": "<255",
        "minimum_pt": dict(zip(sides, [None if n is None else round(n*72/dpi, 4)
                                     for n in minima])),
        "limiting_pages": dict(zip(sides, limiting)),
        "pages_outside_one_inch": violations,
        "all_visible_ink_inside_one_inch": not violations,
        "blank_pdf_pages": empty_pages,
        "renderer_messages": messages,
    }


def inspect_links(doc):
    get = doc.xref_get_key
    named = doc.resolve_names()
    actions = Counter()
    border_styles = Counter()
    referenced = set()
    broken = []
    low_margin = []
    visible_count = 0
    citation_count = 0
    visible_citations = 0
    for page in doc:
        kind, value = get(page.xref, "Annots")
        if kind == "xref":
            value = doc.xref_object(int(value.split()[0]))
        elif kind != "array":
            continue
        for item in re.finditer(r"(\d+)\s+\d+\s+R", value):
            obj = int(item[1])
            if get(obj, "Subtype")[1] != "/Link":
                continue
            action = get(obj, "A/S")[1]
            key = "A/D" if action == "/GoTo" else "Dest"
            target_type, target = get(obj, key)
            if action == "null" and target_type != "null":
                action = "/GoTo"
            actions[action] += 1
            if action == "/GoTo" and target_type in ("string", "name"):
                target = target.lstrip("/") if target_type == "name" else target
                referenced.add(target)
                entry = named.get(target)
                if entry is None or not 0 <= entry.get("page", -1) < len(doc):
                    broken.append({"pdf_page": page.number+1, "target": target})
            if action == "/URI":
                referenced.add(get(obj, "A/URI")[1])
            border = get(obj, "Border")
            width = get(obj, "BS/W")
            numbers = re.findall(r"[-+]?(?:\d*\.?\d+)", border[1])
            stroke = float(width[1]) if width[0] in ("int", "real", "float") else (
                float(numbers[2]) if len(numbers) >= 3 else 1.0)
            flags = get(obj, "F")
            hidden = flags[0] == "int" and int(flags[1]) & (1 | 2 | 32)
            visible = stroke > 0 and not hidden and get(obj, "C") != ("array", "[]")
            visible_count += int(visible)
            if action == "/GoTo" and target.startswith("cite."):
                citation_count += 1
                visible_citations += int(visible)
            border_styles[f"{stroke:g} pt"] += 1
            # Rectangles are defined in unrotated, uncropped PDF coordinates.
            window = page.mediabox
            box = [float(t) for t in re.findall(r"[-+]?(?:\d*\.?\d+)", get(obj, "Rect")[1])]
            if len(box) == 4 and min(box[0]-window.x0, window.x1-box[2],
                                     box[1]-window.y0, window.y1-box[3]) < 72-1e-3:
                low_margin.append({"pdf_page": page.number+1, "printed_page": page.get_label(),
                                   "rect": box, "target": (get(obj, "A/URI")[1]
                                                         if action == "/URI" else target)})
    return {
        "annotations": sum(actions.values()), "actions": dict(actions),
        "visible_border_annotations": visible_count, "border_widths": dict(border_styles),
        "citation_annotations": citation_count, "visible_citation_borders": visible_citations,
        "all_citation_borders_visible": visible_citations == citation_count,
        "named_destinations": len(named), "unique_named_or_uri_targets": len(referenced),
        "unresolved_named_internal_links": broken,
        "annotation_rectangles_inside_one_inch": not low_margin,
        "annotation_rectangles_extending_outside_one_inch": low_margin,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdf", type=Path, help="PDF or build directory containing main.pdf")
    parser.add_argument("--dpi", nargs="+", type=int, default=[144, 288])
    parser.add_argument("--pdftoppm", help="Path to a specific Poppler pdftoppm executable")
    parser.add_argument("--output", type=Path, help="Write JSON here; default is standard output")
    args = parser.parse_args()
    pdf = args.pdf / "main.pdf" if args.pdf.is_dir() else args.pdf
    if not pdf.is_file() or any(d <= 0 for d in args.dpi):
        parser.error("A readable PDF and positive render resolutions are required.")
    program, version = poppler_details(args.pdftoppm)
    with pymupdf.open(pdf) as doc:
        labels = [p.get_label() for p in doc]
        geometry = Counter((tuple(p.mediabox), tuple(p.cropbox), p.rotation) for p in doc)
        result = {
            "pdf": str(pdf.resolve()), "sha256": hashlib.sha256(pdf.read_bytes()).hexdigest(),
            "pages": len(doc), "geometry": [{"media_box": media, "crop_box": crop,
                                             "rotation": angle, "pages": count}
                                            for (media, crop, angle), count in geometry.items()],
            "renderer": {"executable": program, "version": version, "annotations": "visible"},
            "links": inspect_links(doc),
        }
    result["renders"] = [inspect_render(pdf, labels, program, dpi) for dpi in args.dpi]
    result["pass"] = (all(r["all_visible_ink_inside_one_inch"] for r in result["renders"])
                      and not result["links"]["unresolved_named_internal_links"]
                      and result["links"]["all_citation_borders_visible"])
    serial = json.dumps(result, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        args.output.write_text(serial, encoding="utf-8")
    else:
        print(serial, end="")
    return 0 if result["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
