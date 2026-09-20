"""Batch 6 grounding: fetch abstract/full-text evidence for the 10 target papers.

Ladder per spec: Zotero LOCAL API (abstract + .zotero-ft-cache) -> Europe PMC
-> Zotero WEB API (read-only). Saves one txt per paper in batch6\\.
"""
import json, os, re, urllib.parse, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
SAD = os.path.dirname(HERE)
ZOT_LOCAL = "http://localhost:23119/api/users/0"
ZOT_WEB = "https://api.zotero.org/users/631450"
ZOT_KEY = None  # read from credentials file, used only for GET (read-only)
STORAGE = os.path.join(os.path.expanduser("~"), "Zotero", "storage")

txt = open(r"D:\Github\api_credentials_local.txt", encoding="utf-8").read()
m = re.search(r"API key\s*:\s*(\S+)", txt)
ZOT_KEY = m.group(1)

targets = json.load(open(os.path.join(HERE, "batch6_targets.json"), encoding="utf-8"))


def get(url, headers=None):
    try:
        req = urllib.request.Request(url, headers=headers or {})
        with urllib.request.urlopen(req, timeout=45) as r:
            return json.loads(r.read().decode("utf-8", "replace"))
    except Exception as e:
        return {"_err": str(e)}


def epmc(doi):
    q = urllib.parse.quote(f'DOI:"{doi}"')
    d = get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={q}"
            f"&format=json&resultType=core")
    try:
        hit = d["resultList"]["result"][0]
        return hit.get("title", ""), hit.get("authorString", ""), hit.get("journalInfo", {}).get("journal", {}).get("title", ""), hit.get("abstractText", ""), hit.get("id", "")
    except (KeyError, IndexError, TypeError):
        return "", "", "", "", ""


for t in targets:
    key, doi, title = t["zotero_key"], t["doi"], t["title"]
    out, source = [], []

    it = get(f"{ZOT_LOCAL}/items/{key}?format=json")
    abs_local, creators = "", ""
    if "data" in it:
        abs_local = it["data"].get("abstractNote", "")
        creators = "; ".join(
            (c.get("lastName", "") or c.get("name", "")) for c in it["data"].get("creators", []))
        out.append(f"TITLE: {it['data'].get('title', title)}")
        out.append(f"AUTHORS: {creators}")
        out.append(f"ITEMTYPE: {it['data'].get('itemType','')}  YEAR: {it['data'].get('date','')}")
        out.append(f"JOURNAL: {it['data'].get('publicationTitle','') or it['data'].get('proceedingsTitle','') or it['data'].get('bookTitle','')}")
        if abs_local:
            out.append(f"\n[ABSTRACT - Zotero local]\n{abs_local}")
            source.append("zot-local-abs")
        kids = get(f"{ZOT_LOCAL}/items/{key}/children?format=json")
        for ch in kids if isinstance(kids, list) else []:
            ak = ch.get("key", "")
            ft = os.path.join(STORAGE, ak, ".zotero-ft-cache")
            if os.path.exists(ft):
                try:
                    body = open(ft, encoding="utf-8", errors="replace").read()
                    body = re.sub(r"\s+", " ", body)
                    out.append(f"\n[FULLTEXT cache {ak}: {len(body)} chars total; first 3500]\n{body[:3500]}")
                    source.append(f"ft-cache:{ak}")
                except OSError:
                    pass
    else:
        out.append(f"(zotero local miss: {it.get('_err', it)})")

    if not abs_local and doi:
        ttl, auth, jr, ab, pmid = epmc(doi)
        if ab:
            out.append(f"\n[ABSTRACT - Europe PMC pmid={pmid}]\n{jr} | {auth}\n{ab}")
            source.append("europepmc")
        else:
            out.append(f"(europepmc: no abstract found)")

    if not abs_local and not any(s.startswith("ft-cache") for s in source):
        w = get(f"{ZOT_WEB}/items/{key}?format=json", {"Zotero-API-Key": ZOT_KEY})
        if "data" in w:
            wa = w["data"].get("abstractNote", "")
            if wa:
                out.append(f"\n[ABSTRACT - Zotero web]\n{wa}")
                source.append("zot-web-abs")

    out.insert(0, f"=== {key} | {title} | doi={doi} ===")
    out.append(f"\n[SOURCES: {', '.join(source) if source else 'NONE — NO GROUNDING'}]")
    fn = os.path.join(HERE, f"ground_{key}.txt")
    open(fn, "w", encoding="utf-8").write("\n".join(out))
    print(f"{key}: {' + '.join(source) if source else '*** NO GROUNDING ***'}")
