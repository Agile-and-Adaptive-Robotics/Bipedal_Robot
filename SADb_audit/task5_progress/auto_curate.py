"""Task 5 auto-curation: classify, extract quotes, write notes, and create
Airtable records for all 374 gap papers. Uses PAT for bulk writes.
Env: AT_PAT. Run with anaconda base python.
"""
import csv, json, os, re, sys, time, urllib.parse, urllib.request

HERE = os.path.join(r"D:\Github\Bipedal_Robot\SADb_audit", "task5_progress")
SAD = r"D:\Github\Bipedal_Robot\SADb_audit"
PAT = os.environ["AT_PAT"]
BASE = "appMQTnobUNRytIp7"

def norm_doi(d):
    if not d:
        return ""
    d = str(d).strip().lower()
    for p in ("https://doi.org/", "http://doi.org/", "doi:"):
        if d.startswith(p):
            d = d[len(p):]
    return d

def air(method, path, body=None):
    req = urllib.request.Request(
        ("https://api.airtable.com/v0/" + path) if not path.startswith("http") else path,
        data=json.dumps(body).encode() if body is not None else None,
        headers={"Authorization": "Bearer " + PAT, "Content-Type": "application/json"},
        method=method)
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode())

# load all extracted papers
all_papers = []
for fn in sorted(os.listdir(HERE)):
    if fn.startswith("all_extract_") and fn.endswith(".json"):
        with open(os.path.join(HERE, fn), encoding="utf-8") as f:
            all_papers.extend(json.load(f))
print("loaded:", len(all_papers))

# existing Airtable DOIs (for dedup)
existing = set()
offset = ""
hdr = {"Authorization": "Bearer " + PAT}
while True:
    q = f"{BASE}/Papers?pageSize=100&fields%5B%5D=DOI"
    if offset:
        q += "&offset=" + offset
    req = urllib.request.Request("https://api.airtable.com/v0/" + q, headers=hdr)
    with urllib.request.urlopen(req, timeout=40) as r:
        page = json.loads(r.read().decode())
    for r in page.get("records", []):
        existing.add(norm_doi((r.get("fields") or {}).get("DOI")))
    offset = page.get("offset")
    if not offset:
        break
print("existing Airtable DOIs:", len(existing))

ANIMALS = {"cat": "Cat", "cats": "Cat", "feline": "Cat",
           "dog": "Dog", "canine": "Dog",
           "frog": "Frog", "turtle": "Turtle",
           "human": "Human", "humans": "Human", "people": "Human", "patient": "Human", "subjects": "Human",
           "mouse": "Mice", "mice": "Mice", "murine": "Mice",
           "rat": "Rat", "rats": "Rat",
           "lamprey": "Lamprey", "lampreys": "Lamprey",
           "salamander": "Salamander",
           "zebrafish": "Zebrafish", "larval zebrafish": "Zebrafish",
           "cockroach": "Cockroach", "cockroaches": "Cockroach",
           "stick insect": "Stick Insect", "carausius": "Stick Insect",
           "drosophila": "Insects", "insect": "Insects", "insects": "Insects",
           "locust": "Insects", "cricket": "Insects",
           "crayfish": "Arthropods", "crab": "Arthropods",
           "monkey": "Mammals", "macaque": "Mammals",
           "guinea fowl": "Vertebrates", "bird": "Vertebrates"}

def classify(text, title):
    tl = (title + " " + text[:2000]).lower()
    is_review = bool(re.search(r"\breview\b|\bwe review\b|\bthis review\b|\bsynthesi[sz]e\b|\boverview\b|\bsurvey\b", tl))
    is_model = bool(re.search(r"\bmodel\b|\bsimulation\b|\bsimulat(ed|ing|ions?)\b|\bcomputational\b|\bneuromechanical\b|\bneuromusculoskeletal\b|\bneural network\b|\bcpg\b.*\bmodel\b", tl))
    animals = []
    for pat, label in ANIMALS.items():
        if re.search(r"\b" + pat + r"\b", tl) and label not in animals:
            animals.append(label)
    if "human" in tl or "patient" in tl or "subject" in tl:
        if "Human" not in animals:
            animals.append("Human")
    return is_review, is_model, animals

def extract_quote(text, max_words=35):
    """Extract a meaningful verbatim quote from the text."""
    # find abstract or key findings section
    abstract_match = re.search(r"(?:Abstract|ABSTRACT|Summary)\s*[:\.]?\s*(.*?)(?:Introduction|INTRODUCTION|Keywords|1\s+Introduction|\n\n\n)", text, re.DOTALL | re.IGNORECASE)
    source = abstract_match.group(1).strip() if abstract_match else text[:2500]
    # split into sentences
    sentences = re.split(r"(?<=[.!?])\s+", source)
    # prefer sentences with finding-verbs
    best = ""
    best_score = 0
    for s in sentences:
        s = s.strip()
        words = s.split()
        if len(words) < 8 or len(words) > max_words:
            continue
        score = 0
        for kw in ["show", "found", "demonstrate", "reveal", "suggest", "indicate",
                   "result", "conclude", "establish", "identify", "demonstrate that",
                   "we show", "we found", "these results", "our results"]:
            if kw in s.lower():
                score += 2
        if any(kw in s.lower() for kw in ["proprioceptive", "sensory", "feedback", "afferent",
                                            "locomot", "gait", "walk", "spinal", "central pattern"]):
            score += 1
        if score > best_score:
            best = s
            best_score = score
    if not best and sentences:
        # fallback: longest mid-length sentence
        for s in sentences:
            w = s.split()
            if 10 <= len(w) <= max_words:
                best = s.strip()
                break
    return best[:300] if best else ""

def write_note(title, text, is_review, is_model, animals):
    """Write a 1-2 sentence note grounded in the text."""
    abstract_match = re.search(r"(?:Abstract|ABSTRACT)\s*[:\.]?\s*(.*?)(?:Introduction|INTRODUCTION|1\s+Introduction|\n\n)", text, re.DOTALL | re.IGNORECASE)
    source = abstract_match.group(1).strip() if abstract_match else text[:2500]
    sentences = re.split(r"(?<=[.!?])\s+", source)
    finding = ""
    for s in sentences:
        sl = s.lower().strip()
        if any(kw in sl for kw in ["we show", "we found", "we demonstrate", "we conclude",
                                     "these results", "our results", "here we", "we present",
                                     "this study", "we investigated", "our findings"]):
            finding = s.strip()
            break
    if not finding:
        for s in sentences:
            sl = s.lower()
            if any(kw in sl for kw in ["show", "found", "demonstrate", "reveal", "suggest",
                                        "establish", "propose", "develop"]):
                finding = s.strip()
                break
    if not finding and sentences:
        finding = sentences[min(2, len(sentences)-1)].strip()
    parts = []
    if is_review:
        parts.append("Review of " + _topic(title).lower() + ".")
    if is_model:
        parts.append("Computational model study.")
    if finding:
        # trim to reasonable length
        words = finding.split()
        if len(words) > 45:
            finding = " ".join(words[:45]) + "..."
        parts.append(finding)
    if not parts:
        return ""
    note = " ".join(parts)
    # ensure it ends with period
    if not note.endswith("."):
        note += "."
    return note

def _topic(title):
    words = title.split()
    skip = {"a", "an", "the", "of", "in", "on", "for", "and", "to", "with", "by", "from"}
    return " ".join(w for w in words[:8] if w.lower() not in skip)

def get_verbatim(text):
    return extract_quote(text)

# process
to_create = []
insufficient = []
duplicates = []
review_records = []
model_records = []

for p in all_papers:
    doi = norm_doi(p["doi"])
    if doi in existing:
        duplicates.append(p["att_key"])
        continue
    text = p["text"]
    has_text = len(text) > 300 and not text.startswith("EXTRACT-FAIL") and not text.startswith("NO-STAGED")
    is_review, is_model, animals = classify(text, p["title"]) if has_text else (False, False, [])
    note = write_note(p["title"], text, is_review, is_model, animals) if has_text else ""
    quote = get_verbatim(text) if has_text else ""

    fields = {"Name": p["title"][:200]}
    if doi:
        fields["DOI"] = p["doi"]
    if p["year"] and re.match(r"\d{4}", p["year"]):
        fields["Year"] = int(p["year"][:4])
    if p["surname"]:
        fields["Primary Author"] = p["surname"]
    if animals:
        fields["Animal"] = animals
    if note:
        fields["Notes"] = note
    if p["folders"]:
        fields["Secondary Authors"] = []  # populated later from Zotero if needed

    to_create.append({"att_key": p["att_key"], "doi": doi, "fields": fields,
                      "is_review": is_review, "is_model": is_model,
                      "has_text": has_text, "quote": quote, "animals": animals,
                      "note": note, "title": p["title"], "path": p.get("path", "")})

    if not has_text or not note:
        insufficient.append(p["att_key"])
    if is_review:
        review_records.append(p["att_key"])
    if is_model:
        model_records.append(p["att_key"])

print("to create:", len(to_create))
print("duplicates (skip):", len(duplicates))
print("insufficient (no note):", len(insufficient))
print("reviews:", len(review_records))
print("models:", len(model_records))
print("with notes:", sum(1 for t in to_create if t["note"]))

# write batch payload files for Airtable creation
CH = 50
os.makedirs(os.path.join(HERE, "airtable_batches"), exist_ok=True)
for i in range(0, len(to_create), CH):
    chunk = to_create[i:i + CH]
    records = []
    for t in chunk:
        records.append({"fields": t["fields"]})
    with open(os.path.join(HERE, "airtable_batches", f"batch_{i:03d}.json"), "w", encoding="utf-8") as f:
        json.dump({"records": records}, f, ensure_ascii=False)
print("wrote", (len(to_create) + CH - 1) // CH, "Airtable batch files")

# write full classification for review + Resume support
with open(os.path.join(HERE, "task5_classified.json"), "w", encoding="utf-8") as f:
    json.dump(to_create, f, ensure_ascii=False, indent=1)
with open(os.path.join(HERE, "task5_insufficient.json"), "w", encoding="utf-8") as f:
    json.dump(insufficient, f, ensure_ascii=False, indent=1)

# print first 5 for QC
print("\n=== QC SAMPLE (first 5 with notes) ===")
count = 0
for t in to_create:
    if t["note"] and count < 5:
        print(f"\n[{t['fields'].get('Primary Author','?')} {t['fields'].get('Year','')}] {t['title'][:60]}")
        print(f"  TYPE: {'review' if t['is_review'] else ''}{' model' if t['is_model'] else ''}{': animal study' if t['animals'] and not t['is_model'] and not t['is_review'] else ''}")
        print(f"  ANIMALS: {t['animals']}")
        print(f"  NOTE: {t['note'][:250]}")
        print(f"  QUOTE: {t['quote'][:200]}")
        count += 1
