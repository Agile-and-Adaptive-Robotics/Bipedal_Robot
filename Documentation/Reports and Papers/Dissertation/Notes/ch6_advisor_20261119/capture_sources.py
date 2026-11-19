"""Reconstruct the editor readbacks, checking their Chrome extension fingerprints."""
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
import hashlib
import json

here = Path(__file__).resolve().parent
diss = here.parents[1]
capture = json.loads((here / "LIVE_REVIEW.json").read_text(encoding="utf-8"))
base_zip = diss / "Bolen_Dissertation.zip"
base_sha = hashlib.sha256(base_zip.read_bytes()).hexdigest()
assert base_sha == "bdf32431a4b823961bad1fa9d9c93976ab0957cc071915ef6302d41d4039b5fc"
def fingerprint(s):
    data = s.encode("utf-16-le")
    v = 2166136261
    for i in range(0, len(data), 2):
        v = ((v ^ int.from_bytes(data[i:i+2], "little")) * 16777619) & 0xffffffff
    return [len(data)//2, f"{v:08x}"]
sources = {}
with ZipFile(base_zip) as z:
    for name, expected in capture["files"].items():
        source = z.read("chapters/"+name).decode("utf-8-sig").replace("\r\n", "\n")
        if name == "60-sns.tex":
            for e in capture["edits"]:
                if e["file"] != name:
                    continue
                p, old, new = e["offset"], e["original"], e["proposed"]
                assert source[p:p+len(old)] == old, (p, source[p:p+len(old)], old)
                source = source[:p] + new + source[p+len(old):]
        if name == "10-introduction.tex":
            lines = source.split("\n")
            lines[15] = capture["intro16"]
            lines[55] = lines[55].replace("further route development.", "further BPA routing optimization.")
            source = "\n".join(lines)
        print(name, fingerprint(source), expected, flush=True)
        assert fingerprint(source) == expected, name
        sources[name] = source
with ZipFile(here / "LIVE_SOURCES.zip", "w", compression=ZIP_DEFLATED) as z:
    for name, source in sources.items():
        z.writestr(name,source.encode("utf-8"))
