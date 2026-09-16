"""Create Airtable records for all classified gap papers.
Includes Notes ONLY where extraction quality passes whitespace check.
Uses Airtable REST API directly (PAT). Creates Review/Models records too.
Env: AT_PAT. Batches of 40.
"""
import csv, json, os, re, time, urllib.request

HERE = os.path.join(r"D:\Github\Bipedal_Robot\SADb_audit", "task5_progress")
SAD = r"D:\Github\Bipedal_Robot\SADb_audit"
PAT = os.environ["AT_PAT"]
BASE = "appMQTnobUNRytIp7"

def air_post(path, body):
    req = urllib.request.Request(
        f"https://api.airtable.com/v0/{BASE}/{path}",
        data=json.dumps(body).encode(),
        headers={"Authorization": "Bearer " + PAT, "Content-Type": "application/json"},
        method="POST")
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode())

def text_quality_ok(text):
    if not text or len(text) < 300:
        return False
    words = text[:2000].split()
    if len(words) < 50:
        return False
    # check that most "words" are real words (have vowels, reasonable length)
    real = sum(1 for w in words[:100] if 2 <= len(w) <= 20 and re.search(r"[aeiouAEIOU]", w))
    return real / max(len(words[:100]), 1) > 0.6

def note_quality_ok(note):
    if not note or len(note) < 60:
        return False
    words = note.split()
    # check for garbled text (words with no spaces between them)
    long_words = sum(1 for w in words if len(w) > 25)
    return long_words / max(len(words), 1) < 0.1

classified = json.load(open(os.path.join(HERE, "task5_classified.json"), encoding="utf-8"))
print("classified:", len(classified))

# build record payloads
records = []
review_paper_ids = {}  # att_key -> (surname, year) for Review Papers creation
model_paper_ids = {}

for t in classified:
    fields = {"Name": t["title"][:200]}
    doi = t["fields"].get("DOI", "")
    if doi:
        fields["DOI"] = doi
    yr = t["fields"].get("Year")
    if yr:
        fields["Year"] = yr
    pa = t["fields"].get("Primary Author", "")
    if pa:
        fields["Primary Author"] = pa
    animals = t["fields"].get("Animal", [])
    if animals:
        fields["Animal"] = animals

    note = t.get("note", "")
    if note and note_quality_ok(note):
        fields["Notes"] = note

    records.append({"att_key": t["att_key"], "fields": fields,
                    "is_review": t["is_review"], "is_model": t["is_model"],
                    "surname": t["fields"].get("Primary Author", ""),
                    "year": str(t["fields"].get("Year", "")),
                    "has_note": bool(fields.get("Notes"))})

created_keys = set()
import json as _j
if os.path.exists(os.path.join(HERE, "task5_created.json")):
    for c in _j.load(open(os.path.join(HERE, "task5_created.json"), encoding="utf-8")):
        if "fields" in c:
            k = (c["fields"].get("DOI",""), c["fields"].get("Name",""))
        else:
            k = c.get("att_key","")
        created_keys.add(k)
# deduplicate by Name+DOI within the batch
seen = set(created_keys)
deduped = []
for r in records:
    key_doi = r["fields"].get("DOI", "")
    key_name = r["fields"].get("Name", "")
    key_att = r["att_key"]
    if key_doi and key_doi in {k for k in created_keys if k and "/" in str(k)}:
        continue
    if key_att and key_att in created_keys:
        continue
    if (key_doi, key_name) in created_keys:
        continue
    deduped.append(r)
records = deduped
print("deduped:", len(records), "from", len(records))

# create in batches of 10
created = []
fail_count = 0
for i in range(0, len(records), 10):
    chunk = records[i:i + 10]
    body = {"records": [{"fields": r["fields"]} for r in chunk]}
    try:
        resp = air_post("Papers", body)
        for j, rec in enumerate(resp.get("records", [])):
            created.append({"att_key": chunk[j]["att_key"], "record_id": rec["id"],
                            "is_review": chunk[j]["is_review"],
                            "is_model": chunk[j]["is_model"],
                            "surname": chunk[j]["surname"],
                            "year": chunk[j]["year"],
                            "has_note": chunk[j]["has_note"]})
        print(f"batch {i//10+1}: created {len(resp.get('records', []))}")
    except urllib.error.HTTPError as e:
        print(f"batch {i//10+1} FAIL: {e.code} {e.read().decode()[:200]}")
        fail_count += len(chunk)
    except Exception as e:
        print(f"batch {i//10+1} FAIL: {e}")
        fail_count += len(chunk)
    time.sleep(1)

print(f"created: {len(created)}, failed: {fail_count}")

# create Review Papers and Models records
reviews_to_make = [c for c in created if c["is_review"]]
models_to_make = [c for c in created if c["is_model"]]

for r in reviews_to_make:
    name = f"{r['surname']} {r['year']}" if r["year"] else r["surname"]
    try:
        air_post("Review%20Papers", {"fields": {"Name": name, "Paper": [r["record_id"]]}})
        print("Review:", name)
    except Exception as e:
        print("Review FAIL:", name, str(e)[:80])
    time.sleep(0.3)

for r in models_to_make:
    name = f"{r['surname']} {r['year']}" if r["year"] else r["surname"]
    try:
        air_post("Models", {"fields": {"Name": name, "Paper": [r["record_id"]]}})
        print("Model:", name)
    except Exception as e:
        print("Model FAIL:", name, str(e)[:80])
    time.sleep(0.3)

# save completion state
with open(os.path.join(HERE, "task5_created.json"), "w", encoding="utf-8") as f:
    json.dump(created, f, ensure_ascii=False, indent=1)
with_note = sum(1 for c in created if c["has_note"])
print(f"\nFINAL: created {len(created)} | with notes {with_note} | without notes {len(created)-with_note}")
print(f"Review records: {len(reviews_to_make)} | Model records: {len(models_to_make)}")



