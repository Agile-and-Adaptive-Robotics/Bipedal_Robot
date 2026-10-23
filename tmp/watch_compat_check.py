# Compatibility check: the watch shortlist must keep zotero_add_digest.py and
# the batch3 pattern working unchanged. We exec the REAL creators() source from
# zotero_add_digest.py and run it over the watch shortlist's "authors" fields,
# and compare row key sets against batch3\rr_digest_shortlist.json.
import json, io

SRC = r"D:\Github\Bipedal_Robot\SADb_audit\zotero_add_digest.py"
WATCH = r"D:\Github\Bipedal_Robot\SADb_audit\watch\reports\2026-10_shortlist.json"
BATCH3 = r"D:\Github\Bipedal_Robot\SADb_audit\batch3\rr_digest_shortlist.json"

# extract creators() verbatim from zotero_add_digest.py
src = open(SRC, encoding="utf-8").read()
a = src.index("def creators(")
b = src.index("items = []")
ns = {}
exec(src[a:b], ns)
creators = ns["creators"]

rows = json.load(open(WATCH, encoding="utf-8"))
b3 = json.load(open(BATCH3, encoding="utf-8"))

need = {"title", "authors", "container", "year", "doi"}  # what zotero_add_digest reads
print("batch3 row keys:", sorted(b3[0].keys()))
print("watch  row keys:", sorted(rows[0].keys()))
print("watch rows missing zotero-needed keys:", [r["doi"] for r in rows if not need <= set(r)])
print("watch rows missing fetch field:", [r["doi"] for r in rows if "fetch" not in r])

ok = True
for r in rows[:4]:
    c = creators(r.get("authors", ""))[:6]
    print("creators(%s...) ->" % r["authors"][:30], json.dumps(c[:2], ensure_ascii=False))
    for x in c:
        if "lastName" not in x:
            ok = False
print("all creators parsed with lastName:", ok)
