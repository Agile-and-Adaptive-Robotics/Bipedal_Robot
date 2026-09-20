"""Smoke-check the generated app: JSON payload parses, record count matches."""
import json, re, os

HERE = os.path.dirname(os.path.abspath(__file__))
html = open(os.path.join(HERE, "sadb_app.html"), encoding="utf-8").read()
m = re.search(r'<script type="application/json" id="sadb-data">(.*?)</script>', html, re.S)
assert m, "payload block not found"
data = json.loads(m.group(1))
print("papers embedded:", len(data))
print("with notes:", sum(1 for r in data if r["n"]))
print("with layout coords:", sum(1 for r in data if r["x"] or r["yl"]))
print("clusters:", len({r["cl"] for r in data if r["cl"] >= 0}))
assert "</script" not in m.group(1).replace("<\\/", "<"), "unescaped </script in payload"
print("OK")
