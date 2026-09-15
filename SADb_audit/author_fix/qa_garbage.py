import json, glob
for f in sorted(glob.glob(r"D:\Github\Bipedal_Robot\SADb_audit\author_fix\payloads_*.json")):
    data = json.load(open(f, encoding="utf-8"))
    for i, rec in enumerate(data["records"]):
        fields = rec.get("fields", {})
        sec = fields.get("Secondary Authors") or []
        bad = [s for s in sec if len(s) <= 2 and not any(v in s.lower() for v in "aeiou")]
        prim = fields.get("Primary Author", "")
        if bad or (len(prim) <= 2 and not any(v in prim.lower() for v in "aeiou")):
            print(f, i, rec.get("id") or fields.get("DOI"), prim, bad)
