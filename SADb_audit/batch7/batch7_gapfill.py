"""Batch 7 fallback grounding: Crossref + OpenAlex for Yakovenko 2001 chapter,
OpenAlex title-search for the DOI-less Lewinger 2006 paper."""
import json, re, urllib.parse, urllib.request

def get_json(url):
    req = urllib.request.Request(url, headers={"User-Agent": "SADb-curation/1.0 (mailto:benjamin.bolen@pdx.edu)"})
    with urllib.request.urlopen(req, timeout=45) as r:
        return json.loads(r.read().decode("utf-8", "replace"))

def invert(inv):
    if not inv:
        return ""
    pos = {}
    for w, idxs in inv.items():
        for i in idxs:
            pos[i] = w
    return " ".join(pos[i] for i in sorted(pos))

# 1) Yakovenko chapter
for label, url in [
    ("CROSSREF", "https://api.crossref.org/works/10.1201/9781420038705.ch6"),
    ("OPENALEX", "https://api.openalex.org/works/https://doi.org/10.1201/9781420038705.ch6"),
]:
    try:
        w = get_json(url)
        m = w.get("message", w)
        ab = re.sub(r"<[^>]+>", " ", m.get("abstract", "") or "")
        print(f"=== {label} Yakovenko | container={m.get('container-title','')} ===")
        print(ab[:1500] if ab else "(no abstract)")
    except Exception as e:
        print(f"=== {label} Yakovenko FAILED: {e} ===")
    print()

# 2) Lewinger by title search
q = urllib.parse.quote('title.search:"Sensory Coupled Action Switching Modules"')
try:
    d = get_json(f"https://api.openalex.org/works?filter={q}&per-page=5")
    print("=== OPENALEX title hits:", d.get("meta", {}).get("count"), "===")
    for w in d.get("results", [])[:5]:
        print("-", w.get("doi"), "|", w.get("display_name"), "|", w.get("publication_year"),
              "| cited_by:", w.get("cited_by_count"))
        ab = invert(w.get("abstract_inverted_index"))
        if ab:
            print("  ABSTRACT:", ab[:1600])
except Exception as e:
    print("Lewinger search FAILED:", e)
