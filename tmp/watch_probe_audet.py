import json, urllib.request

TARGET = "10.1152/jn.00240.2026"
url = ("https://api.crossref.org/works?query.author=Audet&filter=from-pub-date:2026-09-01,until-pub-date:2026-10-22"
       "&rows=25&select=DOI,title,author")
req = urllib.request.Request(url, headers={"User-Agent": "SADbMonthlyWatch/1.0"})
d = json.load(urllib.request.urlopen(req, timeout=30))
items = d["message"]["items"]
print("items returned:", len(items), "| total-results:", d["message"]["total-results"])
found = [i for i in items if i.get("DOI", "").lower() == TARGET]
print("target in top-25:", bool(found))
if found:
    print("target authors:", [(a.get("family"), a.get("given")) for a in found[0].get("author", [])][:8])
    print("target title:", (found[0].get("title") or [""])[0][:90])
# also show first few returned DOIs
for i in items[:6]:
    fams = [a.get("family") for a in (i.get("author") or [])][:3]
    print("  item:", i.get("DOI"), fams)
