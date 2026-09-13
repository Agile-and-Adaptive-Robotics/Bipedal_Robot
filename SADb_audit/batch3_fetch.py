#!/usr/bin/env python
# Batch 3 grounding fetch: Zotero local API -> PubMed -> Europe PMC, one text file per paper.
import json, subprocess, os, re

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "batch3")
os.makedirs(OUT, exist_ok=True)

ROWS = [
    ("LWWSHR63", "10.1002/cne.23904"),
    ("ACANM79H", "10.1016/0301-0082(96)00028-7"),
    ("8AUBC3V7", "10.1177/0278364905055381"),
    ("YSB9TKJX", "10.1016/j.conb.2009.09.002"),
    ("JU4WBGTJ", "10.1152/jn.90338.2008"),
    ("QU973HAA", "10.1016/s0079-6123(06)57016-5"),
    ("YJVT4HJU", ""),
    ("DQLYCBH2", "10.1113/jphysiol.2013.261115"),
    ("LNCVWZ8E", "10.1146/annurev.physiol.62.1.723"),
    ("WA46PM2Y", "10.3389/fncom.2013.00014"),
]

def curl_json(url):
    r = subprocess.run(["curl.exe", "-s", "-m", "25", url], capture_output=True)
    try:
        return json.loads(r.stdout.decode("utf-8", "replace"))
    except Exception:
        return None

def curl_text(url):
    r = subprocess.run(["curl.exe", "-s", "-m", "25", url], capture_output=True)
    return r.stdout.decode("utf-8", "replace")

def strip_html(s):
    return re.sub(r"<[^>]+>", " ", s or "")

def zotero_search(lib, query):
    from urllib.parse import quote
    url = "http://localhost:23119/api/%s/items?format=json&limit=8&qmode=everything&q=%s" % (lib, quote(query, safe=""))
    return curl_json(url) or []

def pick(items, doi):
    for it in items:
        d = it.get("data", {})
        if doi and (d.get("DOI", "").strip().lower() == doi.lower()):
            return it
    for it in items:
        d = it.get("data", {})
        if d.get("itemType") not in ("attachment", "note"):
            return it
    return None

def doi_url(doi):
    return doi.replace("/", "%2F").replace("(", "%28").replace(")", "%29")

summary = []
for key, doi in ROWS:
    f = os.path.join(OUT, "ground_%s.txt" % key)
    parts = ["==== %s  DOI=%s ====" % (key, doi or "NONE")]
    if doi:
        zpers = zotero_search("users/0", doi)
        zgroup = zotero_search("groups/735051", doi)
    else:
        t = "Highly mobile robots that run and jump"
        zpers = zotero_search("users/0", t)
        zgroup = zotero_search("groups/735051", t)
    zp = pick(zpers, doi)
    zg = pick(zgroup, doi)
    for label, it in (("PERSONAL", zp), ("GROUP", zg)):
        if it:
            d = it["data"]
            parts.append("\n--- ZOTERO %s key=%s itemType=%s" % (label, it.get("key"), d.get("itemType")))
            parts.append("TITLE: %s" % d.get("title", ""))
            cr = d.get("creators", [])
            parts.append("AUTHORS: %s" % "; ".join(((c.get("lastName", "") or c.get("name", "")) + " " + c.get("firstName", "")).strip() for c in cr[:8]))
            parts.append("DATE: %s  PUB: %s" % (d.get("date", ""), d.get("publicationTitle", "")))
            a = d.get("abstractNote", "")
            parts.append("ZOTERO_ABSTRACT(%d chars): %s" % (len(a), a[:4500]))
        else:
            parts.append("\n--- ZOTERO %s: NOT FOUND" % label)
    if zp:
        ch = curl_json("http://localhost:23119/api/users/0/items/%s/children?format=json" % zp.get("key")) or []
        for c in ch:
            cd = c.get("data", {})
            if cd.get("itemType") == "attachment":
                parts.append("PERSONAL_ATTACHMENT: %s (%s)" % (c.get("key"), cd.get("title", "")))
    if doi:
        es = curl_json("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=" + doi_url(doi) + "%5BAID%5D+OR+" + doi_url(doi) + "%5BDOI%5D&retmode=json")
        ids = (es or {}).get("esearchresult", {}).get("idlist", [])
        if ids:
            pmid = ids[0]
            parts.append("\n--- PUBMED PMID=%s" % pmid)
            parts.append(curl_text("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=%s&rettype=abstract&retmode=text" % pmid)[:6500])
        else:
            parts.append("\n--- PUBMED: no hit")
        ep = curl_json("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:%22" + doi_url(doi) + "%22&format=json&pageSize=2")
        res = ((ep or {}).get("resultList") or {}).get("result") or []
        if res:
            r0 = res[0]
            parts.append("\n--- EUROPEPMC (%s) cited:%s" % (r0.get("pubYear"), r0.get("citedByCount")))
            at = strip_html(r0.get("abstractText", ""))
            parts.append("EPMC_ABSTRACT(%d chars): %s" % (len(at), at[:5500]))
        else:
            parts.append("\n--- EUROPEPMC: no hit")
    txt = "\n".join(parts)
    with open(f, "w", encoding="utf-8") as fh:
        fh.write(txt)
    summary.append("%s: %d chars" % (key, len(txt)))
print("\n".join(summary))
