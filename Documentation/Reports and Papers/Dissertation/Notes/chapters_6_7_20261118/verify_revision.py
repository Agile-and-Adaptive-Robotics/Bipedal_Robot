"""Check a full local build, score the retained W2L signals, and package the seven edited chapter files."""
import argparse
import hashlib
import json
import re
import shutil
import zipfile
from collections import Counter
from pathlib import Path

import numpy as np
import pymupdf
from PIL import Image, ImageDraw


HERE = Path(__file__).resolve().parent
DIS = HERE.parents[1]
ROOT = DIS.parents[2]
SOURCE = DIS / "Final"
FIGURES = DIS / "Figures"
CHAPTERS = ["60-sns", "70-overall", "92-AppendixA", "93-AppendixB", "94-AppendixC", "95-AppendixD", "96-AppendixE"]


def sha(data):
    return hashlib.sha256(data).hexdigest()


def measures(path):
    with np.load(path) as data:
        t = data["t"]
        keep = (t >= 4.5) & (t < 30)
        t, heels = t[keep], data["heelF"][keep]
        strikes = []
        for foot in heels.T:
            on = foot > 20
            primary = []
            for i in np.flatnonzero(on[1:] & ~on[:-1]) + 1:
                if not primary or t[i] - primary[-1] > .45:
                    primary.append(t[i])
            strikes.append(np.asarray(primary))
        left, right = strikes
        phases = [(event-a)/(b-a) for a,b in zip(right, right[1:]) for event in left if a <= event <= b]
        x = data["com"][keep, 0]
        return {"heel_strikes": [len(left), len(right)], "phase": float(np.median(phases)),
                "progress_m": float(x[-1]-x[0]), "source_sha256": sha(path.read_bytes())}


def log_summary(build):
    log = (build/"main.log").read_text(encoding="utf-8", errors="replace")
    blg = (build/"main.blg").read_text(encoding="utf-8", errors="replace")
    prohibited = [s for s in ("undefined on input line", "There were undefined", "multiply defined",
                                "Fatal error", "Overfull ", "Missing character:", "Rerun to get",
                                "Label(s) may have changed") if s in log]
    assert not prohibited, prohibited
    assert not re.search(r"^(?:! |LaTeX Error:)",log,re.M)
    assert not re.search(r"Warning--|error message|I couldn't open",blg,re.I)
    return {token: log.count(token) for token in ("pdfTeX warning", "Underfull \\hbox", "Underfull \\vbox",
                                                  "LaTeX Warning", "Overfull ")}


def check_figures(build):
    inventory = json.loads((DIS/"Notes/figure_cleanup_20261118/COMPILED_FIGURES.json").read_text(encoding="utf-8"))
    expected = {row["chapter"]+"/"+g for row in inventory["rows"] for g in row["graphic"]}
    used = set()
    for line in (build/"main.fls").read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("INPUT "):
            path = (SOURCE/line[6:].strip().strip('"')).resolve()
            if path.is_relative_to(FIGURES):
                relative = path.relative_to(FIGURES).as_posix()
                if relative.split("/")[0] in CHAPTERS:
                    used.add(relative)
    assert used == expected, sorted(expected ^ used)
    assets = {}
    with zipfile.ZipFile(DIS/"Notes/figure_cleanup_20261118/Revised_chapter_figures.zip") as z:
        for rel in sorted(used):
            data = (FIGURES/rel).read_bytes()
            assert z.read("figs/"+rel) == data
            assets[rel] = sha(data)
    return assets


def proof(build, previous):
    """Compare text and raster layout, and leave publication-size page renders for inspection."""
    render = build/"reviewed_pages"
    render.mkdir(exist_ok=True)
    with pymupdf.open(build/"main.pdf") as pdf, pymupdf.open(previous/"main.pdf") as old:
        assert len(pdf) == len(old)
        changed, pixel_changed, images = [], [], []
        mapping = {}
        for i, page in enumerate(pdf):
            prior = old[i]
            if page.get_text() != prior.get_text():
                changed.append(i+1)
            # Graphics, text placement, figure floats, and image content.
            new_pix = page.get_pixmap(dpi=96, alpha=False, colorspace=pymupdf.csGRAY)
            old_pix = prior.get_pixmap(dpi=96, alpha=False, colorspace=pymupdf.csGRAY)
            if new_pix.digest != old_pix.digest:
                pixel_changed.append(i+1)
                pix = page.get_pixmap(dpi=144, alpha=False)
                f = render/f"pdf_{i+1:03}.png"
                pix.save(f)
                images.append(f)
            mapping[i+1] = page.get_label()
        for j in range(0, len(images), 6):
            sheet = Image.new("RGB", (1200, 1070), "#dddddd")
            draw = ImageDraw.Draw(sheet)
            for n, f in enumerate(images[j:j+6]):
                with Image.open(f) as img:
                    img.thumbnail((394, 500))
                    x, y = 400*(n%3), 535*(n//3)
                    sheet.paste(img, (x+(400-img.width)//2,y))
                page = int(f.stem.split("_")[1])
                draw.text((x+16,y+503),f"PDF {page} / printed {mapping[page]}", fill="black")
            sheet.save(render/f"overview_{1+j//6}.jpg",quality=90)
        return {"pages":len(pdf), "pages_with_text_changes":changed, "pages_with_raster_changes_96dpi":pixel_changed,
                "revised_layout":{p:mapping[p] for p in pixel_changed}, "pdf_sha256":sha((build/"main.pdf").read_bytes())}


def main():
    cli = argparse.ArgumentParser()
    cli.add_argument("build", type=Path)
    cli.add_argument("previous_build", type=Path)
    args = cli.parse_args()
    build, previous = args.build.resolve(), args.previous_build.resolve()
    report = {"final_log":log_summary(build), "baseline_log":log_summary(previous), "graphics":check_figures(build)}
    results = ROOT/"Code/MuJoCo_SNS/spinal/reports_20261027_codex"
    recordings={}
    for name, expected in [("control", (.498, 2.239)), ("del_Ib", (.494, 2.302)), ("del_all_feedback", (.223, 1.996))]:
        data = measures(results/f"panel_{name}.npz")
        assert tuple(round(data[k],3) for k in ("phase","progress_m")) == expected
        recordings[name]=data
    report["w2l_regenerated"] = recordings
    report["w2l_protocol"] = {"window_seconds":"4.5 <= t < 30","heel_rising_crossing_N":20,
                              "successive_same_foot_interval_s":">0.45","common_pulse_N":40,
                              "common_pulse_time_s":"8 <= t < 8.15"}
    report.update(proof(build, previous))
    record = json.loads((HERE/"PASSAGE_REVIEW.json").read_text(encoding="utf-8"))
    for p in record["passages"]:
        actual = (SOURCE/"chapters"/p["file"]).read_text(encoding="utf-8").splitlines()[p["line"]-1]
        assert sha(actual.encode()) == p["after_sha256"]
    nums = re.compile(r"(?<![A-Za-z_0-9.])\d+(?:\.\d+)?(?![A-Za-z0-9])")
    change = {}
    for p in record["passages"]:
        formatting = r"\\interlinepenalty=\d+"
        old = Counter(nums.findall(re.sub(formatting,"",p["original"])))
        new = Counter(nums.findall(re.sub(formatting,"",p["proposed_and_applied"])))
        if old != new:
            change[p["id"]] = {"removed_numeric_tokens":dict(old-new),"added_numeric_tokens":dict(new-old)}
    report["numeral_changes_in_edited_paragraphs"] = change
    report["model_identifier_added"] = {"T08":"Gait2392"}
    destination = HERE/"Revised_chapters_and_appendices.zip"
    with zipfile.ZipFile(destination, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for chapter in CHAPTERS:
            data = (SOURCE/"chapters"/(chapter+".tex")).read_bytes()
            info = zipfile.ZipInfo("chapters/"+chapter+".tex", (2026,11,18,0,0,0))
            info.compress_type = zipfile.ZIP_DEFLATED
            z.writestr(info, data)
    with zipfile.ZipFile(destination) as z:
        assert z.testzip() is None and len(z.namelist()) == len(CHAPTERS)
        for name in z.namelist():
            assert z.read(name) == (SOURCE/name).read_bytes()
        old_bundle = DIS/"Notes/revision_20261117/Bolen_Dissertation_revised_sources.zip"
        with zipfile.ZipFile(old_bundle) as old:
            changed = [n for n in z.namelist() if old.read("Final/"+n) != z.read(n)]
            assert changed == ["chapters/60-sns.tex","chapters/70-overall.tex"]
        report["chapter_zip_changes_since_november_17"] = changed
    placements = json.loads((DIS/"Notes/figure_cleanup_20261118/COMPILED_FIGURES.json").read_text(encoding="utf-8"))["rows"]
    graphics_pages = {row["pdf_page"] for row in placements}
    report["previous_figure_pages_with_raster_changes"] = sorted(graphics_pages & set(report["pages_with_raster_changes_96dpi"]))
    assert not report["previous_figure_pages_with_raster_changes"]
    report["chapter_zip"] = {"members":len(CHAPTERS),"sha256":sha(destination.read_bytes())}
    shutil.copyfile(build/"main.pdf", HERE/"Bolen_Dissertation_revised.pdf")
    (HERE/"VERIFICATION.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in report.items() if k not in {"graphics","revised_layout"}},indent=2))


if __name__ == "__main__":
    main()
