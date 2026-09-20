"""OpenAlex abstract lookup for the two Elsevier chapters (last grounding fallback)."""
import json, urllib.request

DOIS = {
    "Thiry 2020": "10.1016/b978-0-12-816477-8.00011-9",
    "Prochazka 2017": "10.1016/b978-0-12-803766-9.00008-7",
}

for label, doi in DOIS.items():
    url = "https://api.openalex.org/works/https://doi.org/" + doi
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "SADb-curation/1.0"})
        w = json.loads(urllib.request.urlopen(req, timeout=45).read().decode())
        inv = w.get("abstract_inverted_index")
        print(f"=== OPENALEX {label} | cited_by={w.get('cited_by_count')} ===")
        if inv:
            pos = {}
            for word, idxs in inv.items():
                for i in idxs:
                    pos[i] = word
            print(" ".join(pos[i] for i in sorted(pos)))
        else:
            print("(no abstract_inverted_index)")
    except Exception as e:
        print(f"=== OPENALEX {label} FAILED: {e} ===")
    print()
