import json
import os

RC = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022\Revision_C"

# a) Grillner 1985 citation counts from the SADb citation graph
cites = json.load(open(r"D:\Github\Bipedal_Robot\SADb_audit\export\sadb_cites.json", encoding="utf-8"))
print("cites file type:", type(cites).__name__, "keys:" if isinstance(cites, dict) else "len:", list(cites.keys())[:6] if isinstance(cites, dict) else len(cites))
export = json.load(open(r"D:\Github\Bipedal_Robot\SADb_audit\export\sadb_export.json", encoding="utf-8"))
grill = [r for r in export if "grillner" in (r.get("primary") or "").lower() and (r.get("year") in (1985, "1985"))]
print("Grillner 1985 record:", [(r["id"], r["title"][:40]) for r in grill])
gid = grill[0]["id"] if grill else None
if gid and isinstance(cites, dict):
    for k in list(cites.keys())[:8]:
        v = cites[k]
        print("key", k, type(v).__name__, (len(v) if hasattr(v, "__len__") else v))
    # try common structures
    for k, v in cites.items():
        if isinstance(v, dict) and gid in v:
            print(k, "->", v[gid] if not isinstance(v[gid], list) else f"list[{len(v[gid])}]")
        elif isinstance(v, list) and v and isinstance(v[0], dict) and "cites" in v[0]:
            pass

# b) figure files exist
need = [
    r"figs\50-database\sadb_pivot_app.png",
    r"figs\50-database\sadb_clusters_app.png",
    r"figs\50-database\sadb_explorer.png",
]
for root, dirs, files in os.walk(os.path.join(RC, "figs")):
    for f in files:
        p = os.path.join(root, f)
        rel = os.path.relpath(p, RC)
        if any(n.split("\\")[-1] in f for n in need) or f in (
            "extensor_cam_route.pdf",
            "PulleySurfaces.pdf",
            "PulleyMargins_resolved.pdf",
            "design_workflow.pdf",
        ):
            print("fig:", rel, os.path.getsize(p))

# g) 'the author' occurrences across Revision_C chapters
for f in sorted(os.listdir(os.path.join(RC, "chapters"))):
    if not f.endswith(".tex"):
        continue
    s = open(os.path.join(RC, "chapters", f), encoding="utf-8", errors="replace").read()
    n = s.count("the author")
    if n:
        i = s.find("the author")
        print("PERSON:", f, n, "x ::", s[max(0, i - 50) : i + 50].replace("\n", " "))

# c) build log
rev = os.path.join(RC, "review")
print("review/ contents:", sorted(os.listdir(rev)))
for cand in ("build.log", "build_log.txt", "main.log"):
    p = os.path.join(rev, cand)
    if os.path.exists(p):
        log = open(p, encoding="utf-8", errors="replace").read()
        print(cand, "chars", len(log), "| undefined:", log.lower().count("undefined"), "| errors '!':", log.count("\n!"))
