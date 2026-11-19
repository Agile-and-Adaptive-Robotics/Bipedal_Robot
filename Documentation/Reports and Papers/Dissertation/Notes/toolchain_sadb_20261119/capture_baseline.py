"""Freeze the November 19 author readback used for the Chapter 5 and toolchain review."""
from pathlib import Path
from datetime import datetime, timezone
from zipfile import ZipFile, ZIP_DEFLATED
import hashlib
import json

here = Path(__file__).resolve().parent
diss = here.parents[1]
repo = diss.parents[2]
scratch = repo / "tmp/toolchain_sadb_20261119"
captured = {"10-introduction": (17472,"a391e6e8"), "50-database": (17621,"fdb3209c"),
            "60-sns": (42546,"2336212d"), "70-overall": (18584,"f9abf0b3"),
            "93-AppendixB": (21978,"d1f0755a"), "95-AppendixD": (70707,"8008e2c2"),
            "96-AppendixE": (6646,"1e97f4d1")}
old = "The legs cycled while suspended, and the model tipped when released."
author = ("The legs performed stepping-like motion while it was suspended in air. However, but the "
          "model fell over when lowered to the ground and support was removed.")

def fingerprint(s):
    data=s.encode("utf-16-le")
    v=2166136261
    for i in range(0,len(data),2):
        v=((v ^ int.from_bytes(data[i:i+2],"little"))*16777619)&0xffffffff
    return len(data)//2,f"{v:08x}"

archive=here/"LIVE_CHAPTER_BASELINE.zip"
if archive.exists():
    with ZipFile(archive) as z:
        for name,expected in captured.items():
            assert fingerprint(z.read(name+".tex").decode())==expected
    print("Verified the frozen editor captures.")
else:
    sources = {}
    for name,expected in captured.items():
        source = (diss/"Final/chapters"/(name+".tex")).read_text(encoding="utf-8")
        if name=="70-overall":
            assert source.count(old)==1
            source=source.replace(old,author)
        assert fingerprint(source)==expected,(name,fingerprint(source),expected)
        sources[name+".tex"]=source.encode()
    with ZipFile(archive,"x",compression=ZIP_DEFLATED) as z:
        for name,source in sources.items():
            z.writestr(name,source)
    (scratch/"thesis_prior.bib").write_bytes((diss/"Final/thesis.bib").read_bytes())
    payload={"captured_utc":datetime.now(timezone.utc).isoformat(),"project":"Bolen-Dissertation",
             "method":"Read-only file selection and copy from the current Overleaf editor in Chrome. FNV-1a over UTF-16 code units after CRLF normalization.",
             "files":captured,"bibliography":{"fingerprint":[149968,"244191ac"],"ong_predicting_2019":True,"Chomsky":False},
             "author_change":{"file":"70-overall.tex","line":61,"local_before":old,"live":author},
             "previous_local_bibliography_sha256":hashlib.sha256((scratch/"thesis_prior.bib").read_bytes()).hexdigest()}
    (here/"LIVE_CAPTURE.json").write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print("Verified all seven chapter readbacks, including Ben's changes on line 61.")
