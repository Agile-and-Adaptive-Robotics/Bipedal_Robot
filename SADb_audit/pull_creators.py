import json, urllib.request

keys = ['JMJIQ9ZZ','6CEVM6YE','6NIJ7VIE','8YCE28D9','4HP4LJDQ','6Z7RBVRC','PXUSZ7NR','ZYA3GGJF','SUEWLVPW','E3U29CXK','R29RM7JI','WGZSJTQH','TQYBQQUG','TBHKUUAK','FNQCAAU7','4S96Y3FG','BLL4S52K','IZ2IA2ML','7ZJV4GTK']

fresh = {o["key"]: o for o in json.load(open(r"D:\Github\Bipedal_Robot\SADb_audit\new_refs_fresh.json", encoding="utf-8"))}

for k in keys:
    with urllib.request.urlopen(f"http://localhost:23119/api/groups/735051/items/{k}?format=json", timeout=30) as r:
        it = json.loads(r.read().decode())
    creators = []
    for c in it["data"].get("creators") or []:
        if c.get("lastName"):
            creators.append(c["lastName"])
        elif c.get("name"):
            creators.append(c["name"].strip().split(" ")[-1])
    fresh[k]["creators"] = creators

json.dump([fresh[k] for k in keys], open(r"D:\Github\Bipedal_Robot\SADb_audit\new_refs_fresh.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
for k in keys:
    print(k, fresh[k]["surname"], "| creators:", fresh[k]["creators"])
