"""Create Ben's approved Feedback tags + link their papers (2026-10-02 'yes').

1. New Feedback records (Name + Paper link only — no Recorders, per standing rule):
   - 'Ia disynaptic excitation'                -> Angel 1996  (selective Ia activation)
   - 'I disynaptic excitation (mixed group I)' -> Angel 2005  (group I, Ia/Ib undissociated)
   - 'Lateral postural reflex'                 -> Karayannidou 2009 (Ben's exact name)
   - 'Vestibular feedback to postural control' -> Kooij 2000
   - 'Visual feedback to postural control'     -> Kooij 2000
   - 'Somatosensory feedback to postural control' -> Kooij 2000
2. Re-point Angel 1996/2005 off 'Ib disynaptic excitation' to the new tags.
"""
import json, re, urllib.request

BASE, PAPERS, FEEDBACK = "appMQTnobUNRytIp7", "tblnnMrZszhboU4uD", "tblot5mo4s5KgN5le"
PAT = re.search(r"PAT\s*:\s*(pat[^\s]+)",
                open(r"D:\Github\api_credentials_local.txt", encoding="utf-8").read()).group(1)


def req(path, body=None, method="GET"):
    r = urllib.request.Request("https://api.airtable.com/v0/" + path,
                               data=json.dumps(body).encode() if body else None,
                               headers={"Authorization": "Bearer " + PAT,
                                        "Content-Type": "application/json"}, method=method)
    with urllib.request.urlopen(r, timeout=90) as resp:
        return json.loads(resp.read().decode())


ANGEL96, ANGEL05, KARAY, KOOIJ = ("recWccr9LDvwi2WzR", "reca10selgQlFrpHA",
                                  "recFd4pICpFKMQZ4h", "rec4qLdYPmyjgG02Z")
new_tags = []
ANGEL96, ANGEL05, KARAY, KOOIJ = ("recWccr9LDvwi2WzR", "reca10selgQlFrpHA",
                                  "recFd4pICpFKMQZ4h", "rec4qLdYPmyjgG02Z")
new_tags = [
    ("Ia disynaptic excitation", [ANGEL96]),
    ("I disynaptic excitation (mixed group I)", [ANGEL05]),
    # Karayannidou 2009 was cut in the 290 -> its paper link is impossible;
    # an unlinked tag would be a junk record. Skipped (restore-on-request).
    ("Vestibular feedback to postural control", [KOOIJ]),
    ("Visual feedback to postural control", [KOOIJ]),
    ("Somatosensory feedback to postural control", [KOOIJ]),
]

try:
    created0 = req(f"{BASE}/{FEEDBACK}",
                  {"records": [{"fields": {"Name": n, "Papers": p}} for n, p in new_tags],
                   "typecast": False}, "POST")
except urllib.error.HTTPError as e:
    print("CREATE FAILED", e.code, e.read().decode()[:600])
    raise SystemExit(1)
fb_ids = {"Ia disynaptic excitation": "recbqTL7wmqdombxL", "I disynaptic excitation (mixed group I)": "recT2o3boyYUgmCjN"}
# re-point the two Angel papers (replace their Ib tag; keep any other links)
for rid, tag in ((ANGEL96, "Ia disynaptic excitation"),
                 (ANGEL05, "I disynaptic excitation (mixed group I)")):
    cur = req(f"{BASE}/{PAPERS}/{rid}")["fields"].get("Feedback", [])
    keep = [x["id"] if isinstance(x, dict) else x for x in cur
            if not (isinstance(x, dict) and x.get("name") == "Ib disynaptic excitation")]
    keep.append(fb_ids[tag])
    try:
        req(f"{BASE}/{PAPERS}",
            {"records": [{"id": rid, "fields": {"Feedback": keep}}], "typecast": False}, "PATCH")
        print(f"re-pointed {rid} -> {tag} (kept {len(keep)-1} other links)")
    except urllib.error.HTTPError as e:
        print(f"RE-POINT FAIL {rid}:", e.read().decode()[:250])

# verify echoes
for rid in (ANGEL96, ANGEL05, KARAY, KOOIJ):
    f = req(f"{BASE}/{PAPERS}/{rid}")["fields"].get("Feedback", [])
    print(rid, "feedback now:", [x.get("name") if isinstance(x, dict) else x for x in f])
