"""Smoke-check the generated app (2026-09-28): payloads parse, counts match,
and the main JS block is extracted to _app_main.js for `node --check`.
2026-10-23: asserts the FULL pre-prune base counts — 1,393 records = 943 papers
(382 archived by the prune) + 450 non-paper records (54 archived deleted
orphans); archived total 436, live 957."""
import json, re, os

HERE = os.path.dirname(os.path.abspath(__file__))
html = open(os.path.join(HERE, "sadb_app.html"), encoding="utf-8").read()

m = re.search(r'<script type="application/json" id="sadb-data">(.*?)</script>', html, re.S)
assert m, "payload block not found"
data = json.loads(m.group(1))
print("records embedded:", len(data))

papers = [r for r in data if r.get("rt") == "paper"]
nonp = [r for r in data if r.get("rt") != "paper"]
arch_p = sum(1 for r in papers if r["ar"])
arch_n = sum(1 for r in nonp if r["ar"])
by_type, by_type_arch = {}, {}
for r in nonp:
    by_type[r["rt"]] = by_type.get(r["rt"], 0) + 1
    if r["ar"]:
        by_type_arch[r["rt"]] = by_type_arch.get(r["rt"], 0) + 1
print("papers:", len(papers), f"({arch_p} archived)")
print("non-paper:", len(nonp), f"({arch_n} archived) by type:", by_type,
      "| archived by type:", by_type_arch)
print("archived total:", arch_p + arch_n, "| live:", len(data) - arch_p - arch_n)
print("with notes:", sum(1 for r in data if r["n"]))
print("with afferents:", sum(1 for r in data if r["af"]))
print("clusters:", len({r["cl"] for r in data if r["cl"] >= 0}))

# Full pre-prune base (verified 2026-10-02 prune accounting) must hold exactly.
ids = [r["id"] for r in data]
assert len(ids) == len(set(ids)), "duplicate record ids in payload"
assert len(data) == 1393, f"total {len(data)} != 1393"
assert len(papers) == 943, f"papers {len(papers)} != 943"
assert arch_p == 382, f"archived papers {arch_p} != 382"
assert len(nonp) == 450, f"non-paper {len(nonp)} != 450"
assert arch_n == 54, f"archived non-paper {arch_n} != 54"
assert arch_p + arch_n == 436, f"archived total {arch_p + arch_n} != 436"
assert len(data) - (arch_p + arch_n) == 957, "live != 957"
assert by_type == {"model": 247, "review": 157, "feedback": 42, "recorder": 4}, \
    f"non-paper by type {by_type} != model 247 / review 157 / feedback 42 / recorder 4"
assert by_type_arch == {"model": 3, "review": 48, "feedback": 3}, \
    f"archived non-paper by type {by_type_arch} != model 3 / review 48 / feedback 3 / recorder 0"
assert all(r["rt"] in ("paper", "review", "model", "feedback", "recorder") for r in data), \
    "record type outside paper/review/model/feedback/recorder"
assert "</script" not in m.group(1).replace("<\\/", "<"), "unescaped </script in payload"

c = re.search(r'<script type="application/json" id="sadb-cites">(.*?)</script>', html, re.S)
cites = json.loads(c.group(1))
print("citing papers:", len(cites), "| pairs:", sum(len(v) for v in cites.values()))

l = re.search(r'<script type="application/json" id="cluster-labels">(.*?)</script>', html, re.S)
labels = json.loads(l.group(1))
print("cluster labels:", len(labels))
assert all("</" not in v for v in labels.values())

blocks = re.findall(r"<script>(.*?)</script>", html, re.S)
assert blocks, "no plain <script> block found"
js = blocks[-1]
open(os.path.join(HERE, "_app_main.js"), "w", encoding="utf-8").write(js)
print(f"extracted {len(js)/1024:.0f} KB of JS; placeholders left:",
      len(re.findall(r"__(DATA|CITES|CLUSTER)__", js)))
print("OK")
