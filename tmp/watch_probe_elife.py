import json, urllib.request

for doi in ("10.7554/elife.108656", "10.7554/elife.108656.2"):
    try:
        d = json.load(urllib.request.urlopen(urllib.request.Request(
            "https://api.crossref.org/works/" + doi,
            headers={"User-Agent": "SADbMonthlyWatch/1.0"}), timeout=30))["message"]
        print(doi, "| type", d.get("type"))
        print("   title:", (d.get("title") or [""])[0][:120])
        rel = d.get("relation") or {}
        print("   relations:", {k: [x.get("id") for x in v] for k, v in rel.items()})
    except Exception as e:
        print(doi, "ERROR", e)
