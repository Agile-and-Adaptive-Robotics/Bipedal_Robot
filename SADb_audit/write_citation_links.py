"""Write the in-corpus citation graph into Airtable (2026-09-28).

For every paper that cites other corpus papers, write its cited record ids
into the Papers self-link field "Cites (in corpus)" (fldqQUpWF6Lbi6hdp).
Airtable auto-reciprocates into "Cited by (in corpus)" (fldeiA2v7KuPeR0b6),
so one write populates both directions.

Source of truth for direction = export/sadb_cites.json (OpenAlex-derived,
{citerRecordId: [citedRecordIds]}). Re-runnable: each run overwrites the
field with the freshly computed list. Stdlib only, PAT from
D:\Github\api_credentials_local.txt. Run: myo python write_citation_links.py
"""
import json, os, re, time, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = "appMQTnobUNRytIp7"
PAPERS = "tblnnMrZszhboU4uD"
FLD_CITES = "fldqQUpWF6Lbi6hdp"   # "Cites (in corpus)"


def pat():
    txt = open(r"D:\Github\api_credentials_local.txt", encoding="utf-8").read()
    m = re.search(r"PAT\s*:\s*(pat[^\s]+)", txt)
    if not m:
        raise SystemExit("PAT not found in credentials file")
    return m.group(1)


_PAT = pat()


def air_patch(body):
    req = urllib.request.Request(
        f"https://api.airtable.com/v0/{BASE}/{PAPERS}",
        data=json.dumps(body).encode(),
        headers={"Authorization": "Bearer " + _PAT, "Content-Type": "application/json"},
        method="PATCH")
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read().decode())


cites = json.load(open(os.path.join(HERE, "export", "sadb_cites.json"), encoding="utf-8"))
work = {ra: rbs for ra, rbs in cites.items() if rbs}
print(f"papers with in-corpus references: {len(work)}; "
      f"total directed pairs: {sum(len(v) for v in work.values())}")

items = list(work.items())
done, failed = 0, 0
for i in range(0, len(items), 10):
    chunk = items[i:i + 10]
    body = {"records": [{"id": ra, "fields": {FLD_CITES: rbs}} for ra, rbs in chunk],
            "typecast": False}
    for attempt in range(3):
        try:
            air_patch(body)
            done += len(chunk)
            break
        except Exception as e:
            print(f"  batch {i//10} retry {attempt}: {e}")
            time.sleep(3 + 3 * attempt)
    else:
        failed += len(chunk)
        print(f"  FAILED batch starting at {i}")
    time.sleep(0.35)

print(f"wrote citation links on {done} papers ({failed} failed)")
