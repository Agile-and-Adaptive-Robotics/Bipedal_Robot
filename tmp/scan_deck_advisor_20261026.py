# Scan defense pptx decks for advisor-name strings (Hunt / Alexand).
import re
import zipfile

DECKS = [
    r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation"
    r"\Defense\Bolen_Dissertation_Defense.pptx",
    r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation"
    r"\Defense\AE_Revision_20261023\output"
    r"\Bolen_Defense_Assertion_Evidence_20261023_v6.pptx",
]
for f in DECKS:
    z = zipfile.ZipFile(f)
    hits = []
    for n in z.namelist():
        if not (n.startswith("ppt/slides/slide") and n.endswith(".xml")):
            continue
        s = z.read(n).decode("utf8", "replace")
        for m in re.findall(r"[^<>]{0,60}(?:[Hh]unt|lexand)[^<>]{0,60}", s):
            if "untrigger" not in m.lower():
                hits.append((n.split("/")[-1], m.strip()))
    print("===", f.split("\\")[-1], "===")
    for h in hits[:12]:
        print("  ", h[0], "|", h[1])
    if not hits:
        print("   (no Hunt/Alexand strings)")
