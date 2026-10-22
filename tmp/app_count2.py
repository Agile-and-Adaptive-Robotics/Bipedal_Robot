import collections
import json

d = json.load(open(r"D:\Github\Bipedal_Robot\SADb_audit\export\sadb_export.json", encoding="utf-8"))
prune = collections.Counter(str(r.get("prune")) for r in d)
print("prune flag counts over 943:", prune.most_common())
arch = collections.Counter(str(r.get("archived")) for r in d)
print("(json has 'archived'? )", arch.most_common()[:2])

# built app payload: extract the const DATA array length
html = open(r"D:\Github\Bipedal_Robot\SADb_audit\app\sadb_app.html", encoding="utf-8", errors="replace").read()
i = html.find("const DATA")
seg = html[i : i + 200]
print("payload head:", seg[:120].replace("\n", " "))
# count title occurrences inside payload region as a rough record count
payload = html[i:]
print("records in built payload (by '\"title\":'):", payload.count('"title":'))
