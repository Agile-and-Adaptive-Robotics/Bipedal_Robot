# One-off: verify the group-stamping regen changed ONLY grp/groups
# (nodes byte-identical otherwise, edges/notes untouched).
# BASELINE REFRESHED 2026-09-25 to the pre-w2lvar/syn6 sidecar (the
# original 2026-09-24 snapshot predates Ben's corrected rules/
# shevtsova/shinohara exports of that same night, so it could no longer
# reproduce; semantics unchanged: pre-existing entries must show ZERO
# content diffs, grp fields excepted).
import json, io, sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
import os
a = json.load(open(os.path.join(os.environ["TEMP"], "tpl_pre_groups.json"),
                   encoding="utf-8"))
b = json.load(open("connectome_templates.json", encoding="utf-8"))
# 2026-09-25: the library gained the w2lvar/syn6 variant entries, so the
# regenerated file may hold MORE keys — the guard's semantics stay
# "pre-existing templates unchanged" (subset check + per-key content).
assert set(a) <= set(b), \
    "template key sets differ: pre-groups keys missing: %s" % \
    (set(a) - set(b),)
strip = lambda n: {k: v for k, v in n.items() if k not in ("grp",)}
diffs = 0
for k in a:
    na, nb = a[k]["nodes"], b[k]["nodes"]
    assert len(na) == len(nb), (k, "node count changed")
    assert a[k].get("edges", []) == b[k].get("edges", []), (k, "edges")
    assert a[k].get("_note") == b[k].get("_note"), (k, "note")
    for x, y in zip(na, nb):
        if strip(x) != strip(y):
            diffs += 1
            print("CONTENT DIFF", k, x.get("label"), x, "->", y)
    stamped = sum(1 for n in nb if n.get("grp"))
    print("%-12s grp-stamped %3d/%3d" % (k, stamped, len(nb)))
print("content diffs (excl grp):", diffs)
