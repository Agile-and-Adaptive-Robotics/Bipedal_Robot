"""Build an additive Overleaf update and reproduce it against the untouched native export."""
import argparse
from collections import defaultdict
from difflib import unified_diff
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
from tempfile import NamedTemporaryFile
from zipfile import ZipFile, ZIP_DEFLATED

HERE = Path(__file__).resolve().parent
DISS = HERE.parents[1]
REPO = DISS.parents[2]
FINAL = DISS / "Final"
FIGS = DISS / "Figures"
SCRATCH = REPO / "tmp/figure_names_overview"
HISTORY = HERE.parent / "figure_names_and_overview_20261119"
DATE = re.compile(r"[_-]?(?:19|20)\d{2}[-_]?[01]\d[-_]?[0-3]\d")
INCLUDE = re.compile(r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}")
CURATED = {"50-database", "60-sns", "70-overall", "92-AppendixA", "93-AppendixB",
           "94-AppendixC", "95-AppendixD", "96-AppendixE"}
PANEL_GROUPS = ("pinned_extensor_and_cam", "physical_validation_20261119")
VALIDATION = ("pinned_extensor_torque", "cam_extensor_route", "extensor_validation", "flexor_validation")
HISTORICAL_SUPERSEDED = {
    "30-torque/ExtPin_complete_publication.pdf": "30-torque/pinned_extensor_torque.pdf",
    "30-torque/cam_current_legend.pdf": "30-torque/cam_extensor_route.pdf",
    "30-torque/Ext20mm_physical_publication.pdf": "30-torque/extensor_validation.pdf",
}
EDITABLE = "30-torque/Components/flexor_comparison.fig"


def sha(data):
    return hashlib.sha256(data).hexdigest()


def normal(data):
    return data.replace(b"\r\n", b"\n")


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path, content):
    path.write_text(json.dumps(content, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def read_native():
    archive = DISS / "Bolen_Dissertation.zip"
    captured = SCRATCH / "live_before.zip"
    with ZipFile(archive) as root, ZipFile(captured) as live:
        members = {m.filename: root.read(m) for m in root.infolist() if not m.is_dir()}
        live_names = [m.filename for m in live.infolist() if not m.is_dir()]
        assert len(members) == len(set(live_names)) == len(live_names)
        assert members.keys() == set(live_names)
        assert all(value == live.read(name) for name, value in members.items())
    return members, {"root": sha(archive.read_bytes()), "live": sha(captured.read_bytes()),
                     "identical_member_contents": len(members)}


def safe_path(base, relative):
    name = PurePosixPath(relative)
    if name.is_absolute() or ".." in name.parts or ":" in relative or "\\" in relative:
        raise ValueError(relative)
    path = (base / relative).resolve()
    if not path.is_relative_to(base.resolve()):
        raise ValueError(relative)
    return path


def artwork(relative):
    return safe_path(FIGS, relative.removeprefix("figs/"))


def publication_sources(base):
    output = {}
    candidates = [FINAL / n for n in ("main.tex", "psutex.cls")]
    candidates += sorted((FINAL / "chapters").glob("*.tex"))
    for path in candidates:
        name = path.relative_to(FINAL).as_posix()
        data = path.read_bytes()
        if name not in base or normal(data) != normal(base[name]):
            output[name] = data
    entry = (HERE / "Chomsky_2023.bib").read_bytes()
    native_bib = base["thesis.bib"]
    local_bib = (FINAL / "thesis.bib").read_bytes()
    native = normal(native_bib).rstrip(b"\n")
    assert normal(local_bib).rstrip(b"\n").startswith(native + b"\n")
    assert normal(local_bib)[len(native):].strip() == normal(entry).strip()
    assert b"@misc{chomsky_false_promise_2023," not in native_bib
    assert normal(local_bib).count(b"@misc{chomsky_false_promise_2023,") == 1
    output["Chomsky_2023.bib"] = entry
    return output


def graphic_names(source):
    body = re.sub(r"(?<!\\)%[^\n]*", "", source)
    return set(INCLUDE.findall(body))


def project_art():
    active = {}
    files = {}
    for source in sorted((FINAL / "chapters").glob("*.tex")):
        folder = FIGS / source.stem
        names = graphic_names(source.read_text(encoding="utf-8"))
        if names:
            active[source.stem] = sorted(names)
        elif not folder.is_dir():
            continue
        if source.stem in CURATED:
            roots = {p.name for p in folder.iterdir() if p.is_file()
                     and p.suffix.lower() in {".pdf", ".png", ".jpg", ".jpeg", ".eps", ".svg"}}
            assert roots == names, (source.stem, roots ^ names)
        for name in names:
            assert Path(name).name == name and not DATE.search(name), (source, name)
            candidate = folder / name
            assert candidate.is_file(), candidate
            files["figs/" + source.stem + "/" + name] = ("active artwork", candidate)
            for leaf in (candidate.stem + ".alt.txt", candidate.stem + "_alt.txt"):
                alt = folder / leaf
                if alt.is_file():
                    files["figs/" + source.stem + "/" + leaf] = ("active alternative text", alt)
            converted = candidate.with_name(candidate.stem + "-eps-converted-to.pdf")
            if candidate.suffix.lower() == ".eps" and converted.is_file():
                files["figs/" + source.stem + "/" + converted.name] = ("compiled EPS companion", converted)
    assert set(stem + ".pdf" for stem in VALIDATION) <= set(active["30-torque"])
    for stem in VALIDATION:
        assert "figs/30-torque/" + stem + ".alt.txt" in files
    files["figs/" + EDITABLE] = ("editable flexor comparison", FIGS / EDITABLE)
    for folder in PANEL_GROUPS:
        source = FIGS / "30-torque/Deprecated" / folder
        for file in sorted(source.iterdir()):
            if file.is_file() and file.suffix in {".pdf", ".jpg", ".png", ".eps", ".txt"}:
                files["figs/" + file.relative_to(FIGS).as_posix()] = ("former Chapter 4 panel source", file)
    return active, files


def active_font_record(active, assets):
    historic = read_json(HERE / "FONT_EMBEDDING_FILES.json")
    retained = {name: row for name, row in historic.items() if name not in HISTORICAL_SUPERSEDED}
    assert len(historic) == 12 and len(retained) == 9
    for name, row in retained.items():
        assert name.split("/", 1)[1] in active[name.split("/", 1)[0]]
        assert sha((FIGS / name).read_bytes()) == row["after_sha256"], name
    replacements = {}
    for name in VALIDATION:
        path = "figs/30-torque/" + name + ".pdf"
        replacements[path[5:]] = sha(assets[path][1].read_bytes())
    return {"retained_incremental_repairs": retained, "superseded_historical_names": HISTORICAL_SUPERSEDED,
            "current_chapter4_artwork": replacements}


def rename_plan(base, files):
    mapped = {}
    for row in read_json(HISTORY / "FIGURE_RENAMES.json"):
        old = "figs/" + row["chapter"] + "/" + row["old"]
        new = "figs/" + row["chapter"] + "/" + row["new"]
        if row["new"] == "FlxGrp_2x2_publication.fig":
            new = "figs/" + EDITABLE
        elif row["new"] == "Ext20mm_physical_publication.pdf":
            new = "figs/30-torque/extensor_validation.pdf"
        assert old in base and new in files, (old, new)
        mapped[old] = new
    archived = []
    for row in read_json(HISTORY / "ARCHIVED_ART_NAMES.json"):
        old, new = ("figs/" + row[n] for n in ("old", "new"))
        assert sha(artwork(new).read_bytes()) == row["sha256"], new
        archived.append({"from": old, "to": new, "sha256": row["sha256"],
                         "old_path_in_native_export": old in base})
        if old in base:
            assert sha(base[old]) == row["sha256"]
            mapped[old] = new
    return mapped, archived


def original_destination(name, data):
    relative = PurePosixPath(name.removeprefix("figs/"))
    chapter, leaf = relative.parts[0], DATE.sub("", relative.name)
    folder = FIGS / chapter / "Deprecated"
    # The actual native bytes may predate an intermediate local font repair.
    choices = []
    if folder.is_dir():
        for candidate in folder.rglob(leaf):
            if candidate.is_file() and candidate.stat().st_size == len(data) and candidate.read_bytes() == data:
                choices.append(candidate)
    if choices:
        path = min(choices, key=lambda x: (len(x.parts), x.as_posix()))
        return "figs/" + path.relative_to(FIGS).as_posix()
    return f"figs/{chapter}/Deprecated/preceding_live_export/{leaf}"


def movements(base, files, mapped, active):
    old_active = {}
    for source in sorted((FINAL / "chapters").glob("*.tex")):
        source_key = "chapters/" + source.name
        if source_key in base:
            old_active[source.stem] = graphic_names(base[source_key].decode("utf-8"))
    retire = set()
    for chapter, previous in old_active.items():
        for name in previous - set(active.get(chapter, [])):
            stem = Path(name).stem
            for leaf in (name, stem + ".alt.txt", stem + "_alt.txt"):
                if (path := f"figs/{chapter}/{leaf}") in base:
                    retire.add(path)
    for name in base:
        if not name.startswith("figs/") or len(PurePosixPath(name).parts) != 3:
            continue
        if not artwork(name).exists():
            retire.add(name)
    retire.update(mapped)
    changes = []
    targets = {}
    for old, old_data in sorted(base.items()):
        if not old.startswith("figs/"):
            continue
        new = mapped.get(old)
        local = files.get(old)
        replaced = local is not None and old_data != local[1].read_bytes()
        if not replaced and old not in retire:
            continue
        if replaced and old.endswith(".txt") and normal(old_data) == normal(local[1].read_bytes()):
            continue
        destination = None
        if new:
            future = files[new][1].read_bytes() if new in files else artwork(new).read_bytes()
            if old_data == future:
                destination = new
        if destination is None:
            destination = original_destination(old, old_data)
        prior = targets.get(destination, base.get(destination))
        assert prior is None or prior == old_data, (old, destination)
        targets[destination] = old_data
        changes.append({"from": old, "to": destination, "source_sha256": sha(old_data),
                        "operation": "move before uploading", "new_active_path": new or (old if replaced else None)})
    return changes


def write_package(members, target):
    with NamedTemporaryFile(dir=target.parent, suffix=".pending", delete=False) as stream:
        candidate = Path(stream.name)
    try:
        with ZipFile(candidate, "w", compression=ZIP_DEFLATED) as archive:
            for name, data in sorted(members.items()):
                archive.writestr(name, data)
        with ZipFile(candidate) as archive:
            assert archive.testzip() is None and set(archive.namelist()) == set(members)
            assert all(archive.read(name) == data for name, data in members.items())
        candidate.replace(target)
    finally:
        candidate.unlink(missing_ok=True)


def stage_merge(stage, base, members, moves, active, files):
    stage = stage.resolve()
    allowed = (REPO / "tmp").resolve()
    if not stage.is_relative_to(allowed) or stage == allowed or stage.exists():
        raise ValueError(f"Stage in a new directory under {allowed}")
    stage.mkdir(parents=True)
    for name, data in base.items():
        target = safe_path(stage, name)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
    for item in moves:
        source, dest = (safe_path(stage, item[key]) for key in ("from", "to"))
        assert sha(source.read_bytes()) == item["source_sha256"]
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists():
            assert source.read_bytes() == dest.read_bytes()
            source.unlink()
        else:
            source.rename(dest)
    for name, data in members.items():
        if name in {"READ_ME_FIRST.txt", "Chomsky_2023.bib", "ASSET_OPERATIONS.json"}:
            continue
        target = safe_path(stage, name)
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.is_file() and target.read_bytes() != data and name.startswith("figs/"):
            assert name.endswith(".txt") and normal(target.read_bytes()) == normal(data), name
        target.write_bytes(data)
    raw = base["thesis.bib"]
    ending = b"\r\n" if raw.endswith(b"\r\n") else b"\n"
    entry = normal(members["Chomsky_2023.bib"]).replace(b"\n", ending)
    bib = raw + (ending if raw.endswith(ending) else ending * 2) + entry
    (stage / "thesis.bib").write_bytes(bib)
    assert bib.startswith(raw) and normal(bib).strip() == normal((FINAL / "thesis.bib").read_bytes()).strip()
    for name, (_, path) in files.items():
        expected, actual = path.read_bytes(), safe_path(stage, name).read_bytes()
        assert expected == actual, name
    for item in moves:
        assert sha(safe_path(stage, item["to"]).read_bytes()) == item["source_sha256"], item
    tree = stage / "figs"
    tree_files = [item for item in tree.rglob("*") if item.is_file()]
    dated = [item.relative_to(tree).as_posix() for item in tree_files if DATE.search(item.name)]
    assert not dated, dated
    archived_files = [item for item in tree_files if "Deprecated" in item.relative_to(tree).parts]
    for source in sorted((FINAL / "chapters").glob("*.tex")):
        assert normal(source.read_bytes()) == normal((stage / "chapters" / source.name).read_bytes()), source
    return {"stage": str(stage), "source_files": len(list((stage / "chapters").glob("*.tex"))),
            "verified_graphic_and_source_files": len(files), "verified_preserved_moves": len(moves),
            "native_bibliography_bytes_retained": len(raw),
            "complete_figure_tree": {"path": str(tree), "files": len(tree_files),
                                     "total_bytes": sum(item.stat().st_size for item in tree_files),
                                     "archived_files": len(archived_files), "dated_basenames": dated}}


def main():
    cli = argparse.ArgumentParser()
    cli.add_argument("--stage", type=Path)
    cli.add_argument("--manifest-only", action="store_true")
    args = cli.parse_args()
    base, snapshot = read_native()
    sources = publication_sources(base)
    active, files = project_art()
    fonts = active_font_record(active, files)
    mapped, archives = rename_plan(base, files)
    moves = movements(base, files, mapped, active)
    members = dict(sources)
    directory_paths = defaultdict(list)
    assets = {}
    for name, (kind, path) in sorted(files.items()):
        data = path.read_bytes()
        if data == base.get(name):
            continue
        members[name] = data
        directory = str(PurePosixPath(name).parent)
        directory_paths[directory].append(str(path.resolve()))
        assets[name] = {"local_path": str(path.resolve()), "sha256": sha(data), "bytes": len(data), "kind": kind,
                        "only_line_endings_changed": name.endswith(".txt") and name in base and normal(data) == normal(base[name])}
    manifest = {"baseline": snapshot, "upload_paths_by_live_directory": dict(directory_paths),
                "files": assets, "move_existing_before_upload": moves,
                "archived_art_filename_history": archives,
                "sources": sorted(n for n in sources if n.endswith((".cls", ".tex"))),
                "bibliography_addendum": "Chomsky_2023.bib"}
    write_json(SCRATCH / "live_assets.json", manifest)
    print("Live asset manifest:", len(assets), "upload files,", len(moves), "preservation moves,", len(archives), "local archive records.")
    if args.manifest_only:
        return
    write_json(HERE / "ACTIVE_FIGURE_FONTS.json", fonts)
    package_manifest = {"source_export": snapshot, "move_existing_before_upload": moves,
                        "uploads": {k: {a: v for a, v in row.items() if a != "local_path"} for k, row in assets.items()},
                        "conditional_archive_renames": archives}
    members["ASSET_OPERATIONS.json"] = (json.dumps(package_manifest, indent=2) + "\n").encode()
    instructions = [
        "Updating the existing Bolen-Dissertation project", "",
        "1. Preserve native artwork by applying the FROM -> TO moves in ASSET_OPERATIONS.json before uploading files. Its conditional history includes local archive paths only if they exist in your project.",
        *["   " + item["from"] + " -> " + item["to"] for item in moves], "",
        "2. Upload the supplied files under figs/ to the corresponding project folders. The Deprecated/ copies include earlier panel artwork; Components/ holds the editable MATLAB flexor comparison.", "",
        "3. Apply the supplied source files at the following project paths. An exact source diff is in the review folder.",
        *["   " + name for name in manifest["sources"]], "",
        "4. Append the supplied Chomsky_2023.bib entry to the existing thesis.bib if chomsky_false_promise_2023 is absent. Keep the native bibliography.", "",
        "5. Compile main.tex. The coauthored force paper has its own heading page. Chapter 4 uses the pinned extensor and cam composites and the paired physical-validation figures.", ""]
    members["READ_ME_FIRST.txt"] = "\n".join(instructions).encode("utf-8")
    changes = []
    for name in manifest["sources"]:
        before = normal(base.get(name, b"")).decode("utf-8").splitlines(keepends=True)
        after = normal(members[name]).decode("utf-8").splitlines(keepends=True)
        changes.extend(unified_diff(before, after, fromfile="author_export/" + name, tofile="local_revision/" + name))
    (HERE / "OVERLEAF_SOURCE_CHANGES.diff").write_text("".join(changes), encoding="utf-8")
    target = HERE / "Dissertation_Overleaf_updates.zip"
    write_package(members, target)
    record = {"source_export_sha256": snapshot["root"], "native_snapshots": snapshot,
              "total_graphics_audited": sum(len(v) for v in active.values()), "active_images_by_chapter": active,
              "members": {name: sha(data) for name, data in sorted(members.items())},
              "graphic_changes": {name: item["kind"] for name, item in assets.items()},
              "retire_or_rehome": moves, "local_archive_history": archives, "font_repairs": fonts,
              "native_bibliography_is_additive": True}
    if args.stage:
        record["simulated_update"] = stage_merge(args.stage, base, members, moves, active, files)
    write_json(HERE / "TRANSFER.json", record)
    print("Verified package:", len(members), "members;", record["total_graphics_audited"], "included image files.")


if __name__ == "__main__":
    main()
