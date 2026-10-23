"""READ-ONLY audit parser (2026-10-23, SADb audit): verify payload counts and
per-record recordType/archived presence in BOTH the embedded HTML payload and
export/sadb_export.json. Prints only; writes nothing."""
import json, os, re

APP = os.path.join("SADb_audit", "app", "sadb_app.html")
EXP = os.path.join("SADb_audit", "export", "sadb_export.json")

def stats(rows, rt_key, ar_key):
    total = len(rows)
    miss_rt = [r.get("id", "?") for r in rows if rt_key not in r]
    miss_ar = [r.get("id", "?") for r in rows if ar_key not in r]
    papers = [r for r in rows if r.get(rt_key) == "paper"]
    nonp = [r for r in rows if r.get(rt_key) != "paper"]
    ap = sum(1 for r in papers if r.get(ar_key))
    an = sum(1 for r in nonp if r.get(ar_key))
    byt = {}
    for r in nonp:
        byt[r[rt_key]] = byt.get(r[rt_key], 0) + 1
    ids = [r["id"] for r in rows]
    print("  total=%d | papers=%d (arch %d) | non-paper=%d (arch %d) | arch_total=%d | live=%d"
          % (total, len(papers), ap, len(nonp), an, ap + an, total - ap - an))
    print("  non-paper by type: %s" % byt)
    print("  records MISSING %s: %d | records MISSING %s: %d | duplicate ids: %d"
          % (rt_key, len(miss_rt), ar_key, len(miss_ar), total - len(set(ids))))
    return {r["id"] for r in rows}, (total, len(papers), ap, len(nonp), an)

html = open(APP, encoding="utf-8").read()
m = re.search(r'<script type="application/json" id="sadb-data">(.*?)</script>', html, re.S)
assert m, "sadb-data block not found"
rows = json.loads(m.group(1))
print("HTML payload (sadb-data, compact keys rt/ar):")
ids_html, s_html = stats(rows, "rt", "ar")

exp = json.load(open(EXP, encoding="utf-8"))
print("export/sadb_export.json (recordType/archived):")
ids_exp, s_exp = stats(exp, "recordType", "archived")

print("id sets identical HTML vs export:", ids_html == ids_exp)
print("counts identical HTML vs export:", s_exp == s_html)
print("EXPECTATION (1393, 943/382, 450/54) HTML:", s_html == (1393, 943, 382, 450, 54),
      "| export:", s_exp == (1393, 943, 382, 450, 54))
