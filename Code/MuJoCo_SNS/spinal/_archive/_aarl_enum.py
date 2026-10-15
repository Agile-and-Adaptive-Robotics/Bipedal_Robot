"""Enumerate ALL items in the locomotor-relevant AARL collections.
Writes lit_aarl_enum.txt with one line per item, tagged by collection."""
import io
import json
import sys
import urllib.request

KEY = "34tQExoeiKS1yRe3Z1fbTF5z"
GID = "735051"
HDR = {"Zotero-API-Key": KEY, "Zotero-API-Version": "3"}
COLLECTIONS = [
    ("Y9RHVYHT", "Afferent"),
    ("P7B3YT3N", "Sensory Feedback"),
    ("GVN8ELFX", "Ia Afferent Cite"),
    ("MA3TJQ29", "CPGs"),
    ("FJVWT8R9", "Dual-Layered CPGs"),
    ("8QZLFJ8X", "Network Architecture"),
    ("PUS8CILG", "Spikes to Muscle Activation"),
    ("68IMS87T", "Spiking/NonSpiking Networks"),
    ("RK95XFHP", "FSA_Extensions"),
    ("7UEJ9QHI", "Aim 2"),
    ("8ZJVZTWZ", "CHAPTER 2. BACKGROUND"),
    ("5UEQLNWE", "CHAPTER 4. CALCULATING MN ACTIVATIONS"),
    ("DU88N5P6", "Damping in Locomotion"),
    ("68MAL8TP", "npjRobotics"),
    ("T8DEMIYL", "Literature Reviews"),
]

out = io.open(sys.argv[1] if len(sys.argv) > 1 else "lit_aarl_enum.txt",
              "w", encoding="utf-8")
total = 0
for ck, cname in COLLECTIONS:
    start = 0
    n_coll = 0
    while True:
        url = (f"https://api.zotero.org/groups/{GID}/collections/{ck}/items"
               f"?format=json&limit=100&start={start}&itemType=-attachment")
        req = urllib.request.Request(url, headers=HDR)
        with urllib.request.urlopen(req, timeout=60) as r:
            items = json.loads(r.read().decode("utf-8"))
        if not items:
            break
        for it in items:
            d = it["data"]
            if d["itemType"] in ("attachment", "note", "annotation"):
                continue
            yr = (d.get("date", "") or "")[:4]
            cr = it["meta"].get("creatorSummary", "") or ""
            title = (d.get("title", "") or "").replace("\n", " ")[:100]
            out.write(f"[{cname}] {it['key']} {yr} {cr} | {title}\n")
            n_coll += 1
        total += len(items)
        if len(items) < 100:
            break
        start += 100
    out.write(f"--- end {cname} ---\n")
    print(f"{cname}: {n_coll} items", flush=True)
out.write(f"TOTAL items enumerated: {total}\n")
out.close()
print(f"TOTAL: {total}")
