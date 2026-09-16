import json, urllib.request

keys = ['WJ93DY7W','JMJIQ9ZZ','6CEVM6YE','6NIJ7VIE','8YCE28D9','4HP4LJDQ','6Z7RBVRC','PXUSZ7NR','ZYA3GGJF','SUEWLVPW','QG6A6BH4','E3U29CXK','R29RM7JI','WGZSJTQH','TQYBQQUG','TBHKUUAK','GM9GET82','FNQCAAU7','4S96Y3FG','YE5VCDMJ','JEJUI2NS','F55LBFZG','7NTJY5GQ','BLL4S52K','X43H35G9','KUYUIBA7','IZ2IA2ML','JEUZYHH8','DSTCLWMI','PV8Q9QMY','6NNY7PSJ','F6HULPJH','7ZJV4GTK']

def get_surname(creators):
    for c in creators or []:
        if c.get("lastName"):
            return c["lastName"]
        if c.get("name"):
            return c["name"].strip().split(" ")[-1]
    return ""

out = []
for k in keys:
    try:
        with urllib.request.urlopen(f"http://localhost:23119/api/groups/735051/items/{k}?format=json", timeout=30) as r:
            it = json.loads(r.read().decode())
    except Exception as e:
        print("FAIL", k, str(e)[:80])
        continue
    d = it["data"]
    out.append({
        "key": k, "itemType": d.get("itemType"),
        "doi": (d.get("DOI") or "").strip(),
        "title": " ".join((d.get("title") or "").split()),
        "surname": get_surname(d.get("creators")),
        "year": str(d.get("date") or ""),
        "pub": " ".join((d.get("publicationTitle") or "").split()),
        "url": d.get("url") or "",
        "abstract": d.get("abstractNote") or "",
        "numChildren": (it.get("meta") or {}).get("numChildren", 0),
    })

with open(r"D:\Github\Bipedal_Robot\SADb_audit\new_refs_full.json", "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("pulled:", len(out))
for o in out:
    print("  %s %-14s abs=%-5d kids=%s :: %s" % (o["key"], o["itemType"], len(o["abstract"]), o["numChildren"], o["title"][:60]))
