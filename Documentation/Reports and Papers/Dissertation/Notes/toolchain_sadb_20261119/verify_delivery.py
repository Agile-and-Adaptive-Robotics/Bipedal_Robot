"""Compare independent builds, verify the source package, and publish the full dissertation."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
from zipfile import ZipFile

import pymupdf as pdf

from check_print_margins import inspect_links

here=Path(__file__).resolve().parent
diss=here.parents[1]
repo=diss.parents[2]
output=diss/"output/pdf"
full=output/"Bolen_Dissertation_Toolchain_SADb_20261119.pdf"


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


def verify_front_and_figures(book, local):
    toc = book.get_toc()
    intro = next(num - 1 for level, title, num in toc if level == 1 and title == "Introduction")
    assert book[2].get_text().splitlines()[:2] == ["i", "Abstract"] and book[intro].get_label() == "1"
    assert book[1].get_text().strip() == "© 2026 Ben Bolen"
    heading = next(num for level, title, num in toc if level == 1 and title == "Isometric Force Characterization of Braided Pneumatic Actuators")
    first, second = book[heading - 1], book[heading]
    heading_text = " ".join(first.get_text().split())
    assert "Author contributions." in heading_text and "Funding Acquisition, A.H." in heading_text
    assert "This chapter develops and validates" not in heading_text
    assert "This chapter develops and validates" in " ".join(second.get_text().split())
    assert {link["uri"] for link in first.get_links() if link.get("uri")} == {
        "https://doi.org/10.3390/act15050230", "https://pdxscholar.library.pdx.edu/mengin_fac/553/"}
    first_sns = next(num - 1 for level, title, num in toc if level == 1 and "Tools for Constructing" in title)
    after_sns = next(num - 1 for level, title, num in toc if level == 1 and num - 1 > first_sns)
    aux = "".join(path.read_text(encoding="utf-8", errors="replace") for path in (local / "chapters").glob("*.aux"))
    figures = {}
    for label in ["fig:pinned_extensor", "fig:cam_route", "fig:extensor_test", "fig:flexor_test",
                  "fig:sadb_clusters", "fig:sns_reflex_current", "fig:sns_syn6_overlay", "fig:bpa_exploded_wall",
                  "fig:sensorium_hierarchy", "fig:discussion_neuromechanical_workflow"]:
        figures[label] = list(re.search(r"\\newlabel\{" + label + r"\}\{\{([^}]+)\}\{(\d+)\}", aux).groups())
    expected = {"fig:pinned_extensor":"4.9", "fig:cam_route":"4.12", "fig:extensor_test":"4.13",
                "fig:flexor_test":"4.14", "fig:discussion_neuromechanical_workflow":"7.5"}
    assert all(figures[key][0] == number for key, number in expected.items())
    toc_pattern = re.compile(r"\\contentsline \{(chapter|section|subsection)\}\{\\numberline \{([^}]+)\}(.+)\}\{([^{}]+)\}\{[^{}]+\}%")
    headings = [dict(level=m[1], number=m[2], title=m[3], printed_page=m[4])
                for line in (local / "main.toc").read_text(encoding="utf-8").splitlines()
                if (m := toc_pattern.fullmatch(line))]
    numbers = {row["title"]:row["number"] for row in headings}
    assert numbers["Tracing a Timing Error in the Six-Synergy Controller"] == "6.5.3"
    assert numbers["SNS Tools and Neuromechanical Experiments"] == "7.2.4"
    named = book.resolve_names()
    front_lists = {}
    for category, extension in (("figure", "lof"), ("table", "lot")):
        pattern = re.compile(r"^\\contentsline \{" + category + r"\}\{\\numberline \{([^}]+)\}.*\}\{([^{}]+)\}\{([^{}]+)\}%")
        entries = []
        for line in (local / ("main." + extension)).read_text(encoding="utf-8").splitlines():
            if not line.startswith(r"\contentsline {" + category + "}"):
                continue
            number, folio, target = pattern.fullmatch(line).groups()
            assert book[named[target]["page"]].get_label() == folio, (number, folio, target)
            entries.append(number)
        front_lists[category] = len(entries)
    return {
        "front_list_entries_with_matching_folios": front_lists,
        "front_openings": {title:book[num - 1].get_label() for level, title, num in toc if level == 1 and num <= intro + 1},
        "chapter_3": {"publication_heading_printed_page":first.get_label(), "prose_start_printed_page":second.get_label()},
        "chapter_6": {"printed_pages":[book[first_sns].get_label(), book[after_sns - 1].get_label()],
                      "page_count":after_sns - first_sns},
        "figure_number_and_folio": figures, "numbered_chapters_and_sections": headings,
    }


def verify_package():
    update = json.loads((here / "TRANSFER.json").read_text(encoding="utf-8"))
    assert Path(update["simulated_update"]["stage"]).is_dir()
    package = here / "Dissertation_Overleaf_updates.zip"
    with ZipFile(package) as archive:
        assert set(archive.namelist()) == set(update["members"])
        assert all(hashlib.sha256(archive.read(name)).hexdigest() == digest for name, digest in update["members"].items())
    same_chapters = []
    from_export = diss / "Bolen_Dissertation.zip"
    assert sha(from_export) == update["native_snapshots"]["root"]
    with ZipFile(from_export) as archive:
        for file in sorted((diss / "Final/chapters").glob("*.tex")):
            key = "chapters/" + file.name
            if key in update["members"]:
                continue
            assert file.read_bytes().replace(b"\r\n", b"\n") == archive.read(key).replace(b"\r\n", b"\n"), key
            same_chapters.append(key)
    return {"package_simulation":update["simulated_update"], "native_exports":update["native_snapshots"],
            "chapter_sources_matching_author_export":same_chapters,
            "source_files_in_update":sorted(n for n in update["members"] if n.endswith((".tex", ".cls")))}


def verify_axis_revision(old_build, local, old_manifest):
    import ast
    import numpy as np
    from scipy.io import loadmat

    former = json.loads(old_manifest.read_text(encoding="utf-8"))["members"]
    current = json.loads((here / "TRANSFER.json").read_text(encoding="utf-8"))["members"]
    assert set(former) == set(current)
    changed_files = sorted(name for name in current if current[name] != former[name])
    pair = "figs/30-torque/pinned_extensor_torque"
    assert changed_files == ["ASSET_OPERATIONS.json", pair + ".alt.txt", pair + ".pdf"]
    folder = diss / "Figures/30-torque/Components/extensor_figures"
    builder = (folder / "draw_extensor_figures.py").read_text(encoding="utf-8")
    settings = re.search(r"ax\.set\(xlim=(\([^()]+\)), ylim=(\([^()]+\)),", builder)
    xlim, ylim = (ast.literal_eval(x) for x in settings.groups())
    assert xlim == (-135, 50) and ylim == (-2, 25)
    data = loadmat(folder / "extensor_data.mat", simplify_cells=True)["tests"]
    assert len(data) == 6
    series = {}
    for key in ("measured", "hybrid"):
        x = np.concatenate([np.asarray(row[key + "_angle"]).ravel() for row in data])
        y = np.concatenate([np.asarray(row[key]).ravel() for row in data])
        assert len(x) == len(y) == 113 and np.isfinite(x).all() and np.isfinite(y).all()
        assert xlim[0] < x.min() <= x.max() < xlim[1] and ylim[0] < y.min() <= y.max() < ylim[1]
        series[key] = {"points":len(x), "angle_range_deg":[float(x.min()), float(x.max())],
                       "torque_range_Nm":[float(y.min()), float(y.max())],
                       "angles_below_former_minimum":int((x < -120).sum()),
                       "torques_below_former_minimum":int((y < 0).sum())}
    drawing = json.loads((folder / "verification.json").read_text(encoding="utf-8"))
    assert drawing["graphics"]["pinned_extensor_torque"]["pdf_sha256"] == current[pair + ".pdf"]
    vectors = drawing["pinned_final_artwork"]
    assert tuple(vectors["common_angle_limits_deg"]) == xlim and tuple(vectors["common_torque_limits_Nm"]) == ylim
    panels = vectors["panels"]
    assert len(panels) == len(data)
    for key in ("measured", "hybrid"):
        assert sum(p["series"][key]["printed_markers"] for p in panels) == series[key]["points"]
    finite = sum(p["series"][key]["printed_finite_model_samples"] for p in panels for key in ("optimized", "unoptimized"))
    supplied = sum(np.isfinite(row[key]).sum() for row in data for key in ("optimized", "unoptimized"))
    assert finite == supplied
    outer_clearance = min(p["series"][key]["minimum_ink_clearance_pt"] for p in panels for key in ("hybrid", "measured"))
    assert outer_clearance > 0.0
    with pdf.open(old_build / "main.pdf") as before, pdf.open(local / "main.pdf") as final_doc:
        page_count = len(final_doc)
        assert len(before) == page_count and before.get_toc() == final_doc.get_toc()
        text_changes, raster_changes = [], []
        for old_page, page in zip(before, final_doc):
            if old_page.get_text() != page.get_text():
                text_changes.append(page.number + 1)
            if old_page.get_pixmap().digest != page.get_pixmap().digest:
                raster_changes.append(page.number + 1)
        found = [p.number + 1 for p in final_doc if p.get_label() == "105"]
        assert len(found) == 1 and raster_changes == found and set(text_changes) <= set(found)
        pix = final_doc[found[0] - 1].get_pixmap(matrix=pdf.Matrix(1.5, 1.5))
        pix.save(local / "figure_4_9_printed_page.png")
    return {"updated_package_members":changed_files, "former_final_page_count":page_count,
            "printed_page":105, "physical_page":found[0], "text_changed_physical_pages":text_changes,
            "raster_changed_physical_pages":raster_changes, "source_x_limits_deg":xlim,
            "source_y_limits_Nm":ylim, "plotted_series":series,
            "finite_model_vertices_verified":finite, "minimum_native_marker_ink_clearance_pt":outer_clearance,
            "previous_figure_sha256":former[pair + ".pdf"], "current_figure_sha256":current[pair + ".pdf"]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("local", type=Path)
    parser.add_argument("transfer", type=Path)
    parser.add_argument("--margins", required=True, type=Path,
                        help="Annotated Poppler margins report for the exact local PDF")
    parser.add_argument("--before-axis-proof", type=Path)
    parser.add_argument("--before-axis-package", type=Path)
    args = parser.parse_args()
    if bool(args.before_axis_proof) != bool(args.before_axis_package):
        parser.error("Supply both the preceding proof and its transfer manifest.")
    local, other = args.local.resolve(), args.transfer.resolve()
    evidence = {"local_build":check_build(local), "transfer_build":check_build(other),
                "printed_workflow_sizes":printed_workflow_sizes(local)}
    local_pdf = local / "main.pdf"
    margins = json.loads(args.margins.read_text(encoding="utf-8"))
    assert margins["pass"] and margins["sha256"] == sha(local_pdf)
    assert margins["links"]["annotation_rectangles_inside_one_inch"]
    assert all(not row["pages_outside_one_inch"] for row in margins["renders"])
    assert {144, 288, 300} <= {row["dpi"] for row in margins["renders"]}
    evidence["physical_margins"] = {row["dpi"]:row["minimum_pt"] for row in margins["renders"]}
    fonts = json.loads((here / "ACTIVE_FIGURE_FONT_AUDIT.json").read_text(encoding="utf-8"))
    assert fonts["combined_dissertation"]["pdf_sha256"] == sha(local_pdf)
    assert fonts["combined_dissertation"]["unembedded_resources"] == 0
    evidence["active_figure_fonts"] = {"retained_embedding_repairs":len(fonts["retained_incremental_repairs"]),
        "current_chapter4_artwork":len(fonts["current_chapter4_artwork"]),
        "compiled_input_pdfs":fonts["combined_dissertation"]["input_pdf_count"]}
    with pdf.open(local_pdf) as first, pdf.open(other / "main.pdf") as replay:
        assert len(first) == len(replay) and first.get_toc() == replay.get_toc()
        assert (local / "main.toc").read_bytes() == (other / "main.toc").read_bytes()
        for left, right in zip(first, replay):
            assert left.get_text() == right.get_text(), left.number + 1
            assert left.get_pixmap().digest == right.get_pixmap().digest, left.number + 1
        evidence["full_page_text_and_72dpi_RGB_renders_matched"] = len(first)
        evidence.update(verify_front_and_figures(first, local))
        links = margins["links"]
        assert links == inspect_links(replay) and links["all_citation_borders_visible"]
        assert not links["unresolved_named_internal_links"]
        evidence["visible_link_boxes_in_full_pdf"] = links["visible_border_annotations"]
        evidence["citation_boxes"] = {"all":links["citation_annotations"], "visible":links["visible_citation_borders"]}
    evidence.update(verify_package())
    if args.before_axis_proof:
        evidence["final_figure_4_9_revision"] = verify_axis_revision(args.before_axis_proof, local, args.before_axis_package)
    output.mkdir(parents=True, exist_ok=True)
    publish(local_pdf, full)
    if args.margins.resolve() != (here / "FINAL_PRINT_MARGINS.json").resolve():
        (here / "FINAL_PRINT_MARGINS.json").write_bytes(args.margins.read_bytes())
    evidence["sha256"] = {path.relative_to(diss).as_posix():sha(path)
                          for path in (full, here / "Dissertation_Overleaf_updates.zip")}
    evidence["build_directories"] = [path.relative_to(repo).as_posix() for path in (local, other)]
    (here / "FINAL_VERIFICATION.json").write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    print("Verified full-dissertation native-export replay:", evidence["full_page_text_and_72dpi_RGB_renders_matched"],
          "pages, every caption, all numbered sections, and", evidence["citation_boxes"]["visible"], "visible citation boxes.")


if __name__ == "__main__":
    main()
