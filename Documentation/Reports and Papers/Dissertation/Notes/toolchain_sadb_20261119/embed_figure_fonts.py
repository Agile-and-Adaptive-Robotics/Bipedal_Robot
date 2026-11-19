"""Embed the Helvetica programs used by MuPDF in the dissertation figures.

Usage: python embed_figure_fonts.py input.pdf output.pdf

Use distinct paths. An incremental PDF update preserves the entire input file
and updates only the original unembedded font dictionaries. Requires PyMuPDF
and fontTools. The verified versions and output hashes are in FIGURE_FONTS.json.
"""
import argparse
import hashlib
import os
import shutil
from io import BytesIO
from pathlib import Path
from tempfile import TemporaryDirectory

import pymupdf as pdf
from fontTools.cffLib import CFFFontSet
from fontTools.encodings.StandardEncoding import StandardEncoding
from fontTools.pens.basePen import NullPen


FONT_PROGRAM_SHA256 = {
    "Helvetica": "366720e054067fa49d030f58005c739fc787dcc9727a9327a28767777791d706",
    "Helvetica-Bold": "d4ea80d1d0e29a1e950808009ca8f96fe07a0670b60b556f6b2cf0705337fa9f",
}


def winansi_widths(face):
    """PDF WinAnsi assigns bullets to unused upper slots and aliases space/hyphen."""
    widths = []
    for code in range(32, 256):
        if code in {127, 129, 141, 143, 144, 157}:
            character = "\u2022"
        elif code == 160:
            character = " "
        elif code == 173:
            character = "-"
        else:
            character = bytes([code]).decode("cp1252")
        widths.append(str(round(1000 * face.glyph_advance(ord(character)))))
    return "[" + " ".join(widths) + "]"


def standard_widths(top):
    """Use the original glyph encoding of Helvetica resources without /Encoding."""
    if top.Encoding != "StandardEncoding":
        raise ValueError("The embedded program uses a different glyph encoding")
    widths = []
    for glyph_name in StandardEncoding[32:]:
        charstring = top.CharStrings[glyph_name]
        charstring.draw(NullPen())
        widths.append(str(charstring.width))
    return "[" + " ".join(widths) + "]"


def embed_resource(document, xref, base_name, encoding):
    face = pdf.Font(base_name)
    program = face.buffer
    if hashlib.sha256(program).hexdigest() != FONT_PROGRAM_SHA256[base_name]:
        raise ValueError(f"Unverified {base_name} font program in this PyMuPDF installation")
    cff = CFFFontSet()
    cff.decompile(BytesIO(program), None)
    internal_name = cff.fontNames[0]
    top = cff[internal_name]
    program_xref = document.get_new_xref()
    document.update_object(program_xref, "<< /Subtype /Type1C >>")
    document.update_stream(program_xref, program)
    descriptor_xref = document.get_new_xref()
    bound = " ".join(str(int(value)) for value in top.FontBBox)
    cap_height = round(1000 * face.glyph_bbox(ord("H")).y1)
    x_height = round(1000 * face.glyph_bbox(ord("x")).y1)
    document.update_object(descriptor_xref, (
        f"<< /Type /FontDescriptor /FontName /{internal_name} /Flags 32 "
        f"/FontBBox [{bound}] /ItalicAngle 0 /Ascent {round(1000 * face.ascender)} "
        f"/Descent {round(1000 * face.descender)} /CapHeight {cap_height} "
        f"/XHeight {x_height} /StemV {top.Private.StdVW} /FontFile3 {program_xref} 0 R >>"
    ))
    for key, value in {
        "BaseFont": "/" + internal_name,
        "FirstChar": "32",
        "LastChar": "255",
        "Widths": winansi_widths(face) if encoding == "WinAnsiEncoding" else standard_widths(top),
        "FontDescriptor": f"{descriptor_xref} 0 R",
    }.items():
        document.xref_set_key(xref, key, value)


def embed_figure(source, destination):
    source, destination = Path(source).resolve(), Path(destination).resolve()
    if source == destination:
        raise ValueError("Keep the source and output at distinct paths")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=destination.parent) as work:
        candidate = Path(work) / "figure.pdf"
        shutil.copyfile(source, candidate)
        with pdf.open(candidate) as document:
            missing = {}
            for page in document:
                for xref, ext, kind, name, _, encoding, _ in page.get_fonts(full=True):
                    if ext == "n/a":
                        if kind != "Type1" or encoding not in {"", "WinAnsiEncoding"} or name not in FONT_PROGRAM_SHA256:
                            raise ValueError(f"Unsupported missing font {name}, {kind}, {encoding}")
                        missing[xref] = (name, encoding)
            for xref, (name, encoding) in sorted(missing.items()):
                embed_resource(document, xref, name, encoding)
            if missing:
                document.save(candidate, incremental=True, no_new_id=True, encryption=pdf.PDF_ENCRYPT_KEEP)
        with pdf.open(candidate) as checked:
            for page in checked:
                for font in page.get_fonts(full=True):
                    if not checked.extract_font(font[0])[3]:
                        raise ValueError(f"Unembedded font remains in {source.name}: {font}")
        if not candidate.read_bytes().startswith(source.read_bytes()):
            raise ValueError("Incremental update did not preserve the input")
        os.replace(candidate, destination)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    embed_figure(args.source, args.destination)
