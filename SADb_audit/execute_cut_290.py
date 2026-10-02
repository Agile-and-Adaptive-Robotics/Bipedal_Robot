"""Execute Ben's approved cut list (2026-10-02): delete the 290 from Airtable.

- Snapshots all 290 full export rows first -> archive/deleted_records_20261002.json
  (the app's keep-them mechanism reads the archive/ folder).
- Deletes via the working query-string batch DELETE (bodies get dropped).
- Appends to deletion_log.csv. Papers only; twins stay (Ben holds stubs).
"""
import csv, json, os, re, time, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
BASE, PAPERS = "appMQTnobUNRytIp7", "tblnnMrZszhboU4uD"
ARCHIVE = os.path.join(HERE, "archive")
os.makedirs(ARCHIVE, exist_ok=True)

PAT = re.search(r"PAT\s*:\s*(pat[^\s]+)",
                open(r"D:\Github\api_credentials_local.txt", encoding="utf-8").read()).group(1)


def req_q(path):
    r = urllib.request.Request("https://api.airtable.com/v0/" + path,
                               headers={"Authorization": "Bearer " + PAT}, method="DELETE")
    with urllib.request.urlopen(r, timeout=90) as resp:
        return json.loads(resp.read().decode())


# snapshot from the cached export (full rows)
records = {r["id"]: r for r in json.load(open(os.path.join(HERE, "export", "sadb_export.json"),
                                              encoding="utf-8"))}
ids = [row["record_id"] for row in csv.DictReader(
    open(os.path.join(HERE, "cut_list_290.csv"), encoding="utf-8-sig"))]
snap = [records[i] for i in ids if i in records]
print(f"cut list: {len(ids)} rows | snapshotting {len(snap)} from export")
# move the 10-01 snapshot into archive/ alongside
old = os.path.join(HERE, "deleted_records_20261001.json")
if os.path.exists(old):
    os.replace(old, os.path.join(ARCHIVE, "deleted_records_20261001.json"))
json.dump(snap, open(os.path.join(ARCHIVE, "deleted_records_20261002.json"), "w"),
          ensure_ascii=False, indent=1)

ok, failed = 0, []
t0 = time.time()
for i in range(0, len(ids), 10):
    chunk = ids[i:i + 10]
    try:
        q = "?" + urllib.parse.urlencode([("records[]", x) for x in chunk])
        req_q(f"{BASE}/{PAPERS}{q}")
        ok += len(chunk)
    except Exception as e:
        failed.extend(chunk)
        print(f"  DELETE FAIL batch {i//10}: {str(e)[:90]}")
    time.sleep(0.45)
print(f"deleted {ok}/{len(ids)} in {time.time()-t0:.0f}s; failed {len(failed)}")

with open(os.path.join(HERE, "deletion_log.csv"), "a", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["Papers (the 290 cut list)", len(ids), ok, " ".join(ids),
                "Ben: ax them from Airtable, keep in the app (archive/)",
                time.strftime("%Y-%m-%d %H:%M")])
if failed:
    print("FAILED IDS:", failed)
