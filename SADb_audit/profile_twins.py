"""Profile the Models / Review Papers twin records (2026-10-01, minimal API).

Answers Ben's question with numbers: how many twins are stubs (Name + Paper
link only — safe to delete and re-express as Papers fields) vs. rich (carry
Notes / Rules Used / Papers Cited / Attachments — would lose real data).
Also checks the Feedback table's unused "Models copy" field, and where
"Rules Used" links. ~6 API calls total. Writes twins_profile.json locally.
"""
import json, os, re, urllib.parse, urllib.request

BASE = "appMQTnobUNRytIp7"


def pat():
    return re.search(r"PAT\s*:\s*(pat[^\s]+)",
                     open(r"D:\Github\api_credentials_local.txt", encoding="utf-8").read()).group(1)


def air(path, params=None):
    q = "?" + urllib.parse.urlencode(params, doseq=True) if params else ""
    req = urllib.request.Request("https://api.airtable.com/v0/" + path + q,
                                 headers={"Authorization": "Bearer " + pat()})
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read().decode())


def fetch_all(tid, fields):
    out, offset, calls = [], None, 0
    while True:
        p = [("pageSize", 100), ("fields[]", fields)]
        if offset:
            p.append(("offset", offset))
        d = air(f"{BASE}/{tid}", p); calls += 1
        out.extend(d.get("records", []))
        offset = d.get("offset")
        if not offset:
            return out, calls


profile = {}
total_calls = 0

for label, tid in (("models", "tblsBq9IEv7dZe6fn"), ("reviews", "tblSEubKcRId4wYMK")):
    recs, c = fetch_all(tid, ["Name", "Paper", "Notes", "Rules Used", "Papers Cited", "Attachments"])
    total_calls += c
    stub = rich = orphan = 0
    rich_list, orphan_list = [], []
    for r in recs:
        f = r.get("fields", {})
        has_paper = bool(f.get("Paper"))
        has_data = bool((f.get("Notes") or "").strip() or f.get("Rules Used")
                        or f.get("Papers Cited") or f.get("Attachments"))
        if not has_paper:
            orphan += 1
            orphan_list.append(f.get("Name", r["id"]))
        elif has_data:
            rich += 1
            rich_list.append(f.get("Name", r["id"]))
        else:
            stub += 1
    profile[label] = {"total": len(recs), "stub_name_link_only": stub, "rich": rich,
                      "orphan_no_paper": orphan, "rich_names": rich_list,
                      "orphan_names": orphan_list}
    print(f"{label}: {len(recs)} total | stubs {stub} | rich {rich} | orphans {orphan}")

fb, c = fetch_all("tblot5mo4s5KgN5le", ["Name", "Papers", "Models", "Models copy"])
total_calls += c
fb_mc = sum(1 for r in fb if r["fields"].get("Models copy"))
profile["feedback_models_copy_used"] = fb_mc
print(f"Feedback: {len(fb)} records; 'Models copy' used on {fb_mc}")

json.dump(profile, open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                     "twins_profile.json"), "w"), ensure_ascii=False, indent=1)
print(f"API calls used: {total_calls}; wrote twins_profile.json")
