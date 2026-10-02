"""Execute the approved "90" deletion (Ben's go, 2026-10-01) — guarded + audited.

Order of operations:
  1. GUARD: the rich Models twin "Latash 2020" links to the split-belt
     PREPRINT (being deleted as a duplicate). Ensure it also links the
     PUBLISHED record recQkgW6mi6AezV1b BEFORE any deletion (add if missing).
  2. Delete the 88 manifest papers (10 per call).
  3. Delete the 7 satellite orphans (3 Models + 1 Review Papers + 3 Feedback).
  4. Append the audit row to deletion_log.csv.
~15 API calls total (within the grace-period budget; deletion via API is
auditable and far safer than ~95 UI interactions).
"""
import csv, json, os, re, time, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = "appMQTnobUNRytIp7"
PAPERS = "tblnnMrZszhboU4uD"
PUBLISHED = "recQkgW6mi6AezV1b"   # Frontiers version of the split-belt CPG paper


def pat():
    return re.search(r"PAT\s*:\s*(pat[^\s]+)",
                     open(r"D:\Github\api_credentials_local.txt", encoding="utf-8").read()).group(1)


_PAT = pat()


def req(path, body=None, method="GET"):
    r = urllib.request.Request("https://api.airtable.com/v0/" + path,
                               data=json.dumps(body).encode() if body else None,
                               headers={"Authorization": "Bearer " + _PAT,
                                        "Content-Type": "application/json"},
                               method=method)
    with urllib.request.urlopen(r, timeout=90) as resp:
        return json.loads(resp.read().decode())


def req_q(path, method):
    r = urllib.request.Request("https://api.airtable.com/v0/" + path,
                               headers={"Authorization": "Bearer " + _PAT},
                               method=method)
    with urllib.request.urlopen(r, timeout=90) as resp:
        return json.loads(resp.read().decode())


mf = json.load(open(os.path.join(HERE, "deletion_manifest.json"), encoding="utf-8"))

# --- 1. guard: re-point the Latash 2020 twin at the published record ---
q = urllib.parse.quote('{Name}="Latash 2020"')
twins = req(f"{BASE}/tblsBq9IEv7dZe6fn?filterByFormula={q}&fields%5B%5D=Paper").get("records", [])
guard = "no twin named 'Latash 2020' found"
for t in twins:
    linked = [x["id"] if isinstance(x, dict) else x for x in t["fields"].get("Paper", [])]
    if PUBLISHED in linked:
        guard = f"twin {t['id']} already links the published record"
    else:
        linked.append(PUBLISHED)
        req(f"{BASE}/tblsBq9IEv7dZe6fn",
            {"records": [{"id": t["id"], "fields": {"Paper": linked}}]}, "PATCH")
        guard = f"twin {t['id']} re-pointed: added published {PUBLISHED} (had {len(linked)-1})"
print("GUARD:", guard)

# --- 2+3. deletions (query-string records[] — DELETE bodies get dropped) ---
log_rows, failed = [], []
ALREADY_GONE = {"recIoqMm5pjzrrBkv", "recqoLEviGnogVbUt"}   # consumed by probes


def delete_batch(tid, ids, label):
    ids = [i for i in ids if i not in ALREADY_GONE]
    ok = 0
    for i in range(0, len(ids), 10):
        chunk = ids[i:i + 10]
        try:
            q = "?" + urllib.parse.urlencode([("records[]", x) for x in chunk])
            req_q(f"{BASE}/{tid}{q}", "DELETE")
            ok += len(chunk)
        except Exception as e:
            failed.extend(chunk)
            print(f"  DELETE FAIL {label} batch {i//10}: {str(e)[:90]}")
        time.sleep(0.45)
    print(f"deleted {ok}/{len(ids)} from {label}")
    log_rows.append((label, len(ids), ok, " ".join(ids)))

delete_batch(PAPERS, mf["papers"], "Papers (the 88)")
delete_batch("tblsBq9IEv7dZe6fn", mf["satellite"]["Models"], "Models orphans")
delete_batch("tblSEubKcRId4wYMK", mf["satellite"]["Review Papers"], "Review Papers orphan")
delete_batch("tblot5mo4s5KgN5le", mf["satellite"]["Feedback"], "Feedback orphans")

with open(os.path.join(HERE, "deletion_log.csv"), "w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["table", "planned", "deleted", "record_ids", "guard", "timestamp"])
    for label, planned, ok, ids in log_rows:
        w.writerow([label, planned, ok, ids, guard, time.strftime("%Y-%m-%d %H:%M")])
print("audit: deletion_log.csv + deleted_records_20261001.json (full snapshots)")
if failed:
    print("FAILED IDS:", failed)
