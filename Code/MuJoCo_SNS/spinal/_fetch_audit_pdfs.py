"""Search + download the audit's key papers from the AARL group library
into spinal/lit_pdfs/ (untracked working folder)."""
import io
import os
import sys
import urllib.parse
import urllib.request

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
KEY = "34tQExoeiKS1yRe3Z1fbTF5z"
GID = "735051"
HDR = {"Zotero-API-Key": KEY, "Zotero-API-Version": "3"}
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "lit_pdfs")
os.makedirs(OUT, exist_ok=True)

import json
import urllib.request as rq


def get(url):
    req = rq.Request(url, headers=HDR)
    with rq.urlopen(req, timeout=90) as r:
        return r.read()


def jget(url):
    return json.loads(get(url).decode("utf-8"))


QUERIES = {
    "pratt1987_iain_renshaw": "Ia inhibitory interneurons Renshaw cells contributors spinal mechanisms fictive locomotion",
    "hultborn1971_recurrent_ia": "Recurrent inhibition motor axon collaterals transmission Ia inhibitory pathway",
    "gossard1994_ib_pathway": "Transmission locomotor-related group Ib pathway hindlimb extensor",
    "conway1987_grI_reset": "Proprioceptive input resets central locomotor rhythm spinal cat",
    "perreault1999_grI_extensor": "Proprioceptive Control of Extensor Activity during Fictive Scratching and Weight Support",
    "perreault2011_grii_reset": "Effects of stimulation of hindlimb flexor group II afferents during fictive locomotion",
    "noga1987_mecamylamine": "role of Renshaw cells in locomotion antagonism of their excitation",
    "jankowska2010_ib_ii": "Functional subdivision of feline spinal interneurons reflex pathways group Ib and II",
    "pearson1998_muscle_afferents": "Enhancement and Resetting of Locomotor Activity by Muscle Afferents",
    "dominguez2020_reset_ins": "Candidate Interneurons Mediating the Resetting of the Locomotor Rhythm by Extensor Group I",
    "mccrea1980_renshaw": "Renshaw cell activity and recurrent effects on motoneurons during fictive locomotion",
}

for tag, q in QUERIES.items():
    out_pdf = os.path.join(OUT, tag + ".pdf")
    if os.path.exists(out_pdf) and os.path.getsize(out_pdf) > 50000:
        print(f"{tag}: already downloaded")
        continue
    try:
        url = (f"https://api.zotero.org/groups/{GID}/items?format=json"
               f"&q={urllib.parse.quote(q)}&limit=5&itemType=-attachment")
        items = jget(url)
        got = False
        for it in items:
            kids = jget(f"https://api.zotero.org/groups/{GID}/items/"
                        f"{it['key']}/children?format=json")
            pdfs = [c for c in kids
                    if c["data"].get("contentType") == "application/pdf"]
            if pdfs:
                raw = get(f"https://api.zotero.org/groups/{GID}/items/"
                          f"{pdfs[0]['key']}/file")
                with open(out_pdf, "wb") as f:
                    f.write(raw)
                print(f"{tag}: PDF {len(raw)//1024} KB  <- {it['key']}")
                got = True
                break
        if not got:
            print(f"{tag}: no PDF in group library")
    except Exception as e:
        print(f"{tag}: ERROR {e}")
