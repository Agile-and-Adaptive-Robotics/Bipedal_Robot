"""Corrected counts (handles escaped \\%) for the Actuators style analysis."""
import re

SRC = r"D:/GitHub/Bipedal_Robot/Documentation/Reports and Papers/Isometric_Force_Test/Final/Final_proofread/manuscript.v5/actuators-4207821-edited.tex"
text = open(SRC, encoding="utf-8", errors="replace").read()

def strip_comments(s):
    out = []
    for l in s.splitlines():
        # remove % comments but keep escaped \%
        l = re.sub(r"(?<!\\)%.*$", "", l)
        if l.strip().startswith("%"):
            l = ""
        out.append(l)
    return "\n".join(out)

body_start = text.index(r"\begin{document}")
body_end = text.index(r"\begin{thebibliography}")
body = strip_comments(text[body_start:body_end])

figs = re.findall(r"\\label\{(fig:[^}]*)\}", body)
tabs = re.findall(r"\\label\{(tab:[^}]*)\}", body)
eqs = re.findall(r"\\label\{(eq:[^}]*)\}", body)
print("figures:", len(figs), figs)
print("tables:", len(tabs), tabs)
print("labeled equations:", len(eqs), eqs)

cites = re.findall(r"\\cite[tp]?\{([^}]*)\}", body)
keys = [k.strip() for c in cites for k in c.split(",")]
print("citation calls:", len(cites), "key instances:", len(keys), "unique:", len(set(keys)))
narrative = len(re.findall(r"\\citet\{", body))
parenth = len(re.findall(r"\\citep\{", body))
bare = len(re.findall(r"\\cite\{", body))
print(f"narrative citet={narrative}, parenthetical citep={parenth}, bare cite={bare}")

figrefs = len(re.findall(r"Figure~\\ref\{", body))
tabrefs = len(re.findall(r"Table~\\ref\{", body))
eqrefs = len(re.findall(r"\\ref\{eq:", body))
print(f"in-text refs: Figure={figrefs}, Table={tabrefs}, eq={eqrefs}")

# abstract word count (abstract is one line in the file)
for l in text.splitlines():
    if l.startswith("\\abstract{"):
        a = l[len("\\abstract{"):].rstrip().rstrip("}")
        a = re.sub(r"\\qtylist?\{[^}]*\}\{[^}]*\}", " X ", a)
        a = re.sub(r"\\[a-zA-Z]+", " ", a)
        a = re.sub(r"[{}\\]", " ", a)
        words = a.split()
        print("ABSTRACT WORDS:", len(words))
        sents = [s for s in re.split(r"\. ", a) if s.strip()]
        print("ABSTRACT SENTENCES:", len(sents))
        break

# sentence-initial patterns in Results (voice check)
sec = re.search(r"\\section\{Results\}(.*?)\\section\{Discussion\}", body, re.S).group(1)
sents = re.split(r"(?<=[.])\s+", re.sub(r"\\[a-zA-Z]+(\[[^\]]*\])?", " ", sec))
openers = []
for s in sents:
    s2 = re.sub(r"[{}\\$]", "", s).strip()
    if s2:
        openers.append(" ".join(s2.split()[:3]))
print("\nResults sentence openers (first 20):")
for o in openers[:20]:
    print("   ", o)
