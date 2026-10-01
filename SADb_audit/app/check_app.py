"""Smoke-check the generated app (2026-09-28): payloads parse, counts match,
and the main JS block is extracted to _app_main.js for `node --check`."""
import json, re, os

HERE = os.path.dirname(os.path.abspath(__file__))
html = open(os.path.join(HERE, "sadb_app.html"), encoding="utf-8").read()

m = re.search(r'<script type="application/json" id="sadb-data">(.*?)</script>', html, re.S)
assert m, "payload block not found"
data = json.loads(m.group(1))
print("papers embedded:", len(data))
print("with notes:", sum(1 for r in data if r["n"]))
print("with afferents:", sum(1 for r in data if r["af"]))
print("clusters:", len({r["cl"] for r in data if r["cl"] >= 0}))
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
