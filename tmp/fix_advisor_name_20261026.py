# Fix advisor name in defense decks, tolerant of a locked file.
import os
import shutil
import zipfile

DECKS = [
    r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation"
    r"\Defense\Bolen_Dissertation_Defense.pptx",
    r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation"
    r"\Defense\AE_Revision_20261023\output"
    r"\Bolen_Defense_Assertion_Evidence_20261023_v6.pptx",
]
BAD = "M. Anthony Hunt"
GOOD = "Alexander J Hunt"
for f in DECKS:
    tmp = f + ".tmp"
    if os.path.exists(tmp):
        os.remove(tmp)
    bak = f + ".bak_wrongname_20261026"
    if not os.path.exists(bak):
        shutil.copy2(f, bak)
    zin = zipfile.ZipFile(f)
    n_fix = 0
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            data = zin.read(item.filename)
            if item.filename.startswith("ppt/slides/slide") and \
               item.filename.endswith(".xml"):
                s = data.decode("utf8", "replace")
                if BAD in s:
                    n_fix += s.count(BAD)
                    data = s.replace(BAD, GOOD).encode("utf8")
            zout.writestr(item, data)
    zin.close()
    try:
        os.replace(tmp, f)
        print(f.split("\\")[-1], ":", n_fix, "replacements OK")
    except PermissionError:
        os.remove(tmp)
        print(f.split("\\")[-1], ": LOCKED by another process -",
              n_fix, "replacements PREPARED but NOT applied;",
              "close PowerPoint/preview and rerun.")
