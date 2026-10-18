"""Page count of the published PDF + Methods voice probe."""
import re

try:
    from pypdf import PdfReader
except ImportError:
    try:
        from PyPDF2 import PdfReader
    except ImportError:
        PdfReader = None

for p in [
    r"D:/GitHub/Bipedal_Robot/Documentation/Reports and Papers/Isometric_Force_Test/IsometricForceCharacterizationofBPAs.pdf",
    r"D:/GitHub/Bipedal_Robot/Documentation/Reports and Papers/Isometric_Force_Test/Final/Final_proofread/manuscript.v5/actuators-4207821-edited.pdf",
]:
    if PdfReader:
        try:
            r = PdfReader(p)
            print(p.split("/")[-1], "->", len(r.pages), "pages; title:", (r.metadata or {}).get("/Title"))
            continue
        except Exception as e:
            print(p, "pypdf failed:", e)
    # fallback: count /Page objects
    data = open(p, "rb").read()
    n = data.count(b"/Type /Page") - data.count(b"/Type /Pages")
    print(p.split("/")[-1], "-> crude /Page count:", n)

# Methods voice: print the Methods section sentences
SRC = r"D:/GitHub/Bipedal_Robot/Documentation/Reports and Papers/Isometric_Force_Test/Final/Final_proofread/manuscript.v5/actuators-4207821-edited.tex"
text = open(SRC, encoding="utf-8", errors="replace").read()
def strip_comments(s):
    return "\n".join(re.sub(r"(?<!\\)%.*$", "", l) for l in s.splitlines())
body = strip_comments(text)
m = re.search(r"\\section\{Materials and Methods\}(.*?)\\section\{Results\}", body, re.S).group(1)
m = re.sub(r"\\begin\{figure\}.*?\\end\{figure\}", " [FIGURE] ", m, flags=re.S)
m = re.sub(r"\\[a-zA-Z]+(\[[^\]]*\])?", " ", m)
m = re.sub(r"[{}\\$]", "", m)
print("\n--- METHODS TEXT (cleaned) ---")
print(re.sub(r"\n+", "\n", m).strip()[:2400])
