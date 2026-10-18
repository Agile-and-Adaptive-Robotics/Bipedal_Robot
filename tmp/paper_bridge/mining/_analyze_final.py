"""Final per-section word counts with correct comment stripping."""
import re

SRC = r"D:/GitHub/Bipedal_Robot/Documentation/Reports and Papers/Isometric_Force_Test/Final/Final_proofread/manuscript.v5/actuators-4207821-edited.tex"
text = open(SRC, encoding="utf-8", errors="replace").read()

def strip_comments(s):
    return "\n".join(re.sub(r"(?<!\\)%.*$", "", l) for l in s.splitlines())

body = strip_comments(text[text.index(r"\begin{document}"):text.index(r"\begin{thebibliography}")])

def wordcount(seg, drop_floats=True):
    s = seg
    if drop_floats:
        s = re.sub(r"\\begin\{figure\}.*?\\end\{figure\}", " ", s, flags=re.S)
        s = re.sub(r"\\begin\{table\}.*?\\end\{table\}", " ", s, flags=re.S)
    s = re.sub(r"\\begin\{(equation|align|adjustwidth|tabularx)\*?\}.*?\\end\{\1\*?\}", " ", s, flags=re.S)
    s = re.sub(r"\\(?:sub)?section\*?\{([^}]*)\}", r" ", s)
    s = re.sub(r"\\[a-zA-Z]+(\[[^\]]*\])?(\{[^{}]*\})*", " ", s)
    s = re.sub(r"[{}\\$&~]", " ", s)
    return len(s.split())

sec_pat = re.compile(r"\\section\{([^}]*)\}")
ms = list(sec_pat.finditer(body))
total = 0
for i, m in enumerate(ms):
    seg = body[m.end(): ms[i+1].start() if i+1 < len(ms) else len(body)]
    w = wordcount(seg)
    total += w
    # paragraph count: blank-line separated blocks of >=40 chars of prose
    paras = [p for p in re.split(r"\n\s*\n", seg) if len(p.strip()) > 40]
    print(f"{m.group(1):<25} {w:>5} words   ~{len(paras)} paragraphs")
print(f"{'TOTAL body prose':<25} {total:>5} words")

# equation locations
for env in ["equation", "align"]:
    for mm in re.finditer(r"\\begin\{" + env + r"\}(.*?)\\end\{" + env + r"\}", body, re.S):
        lab = re.search(r"\\label\{([^}]*)\}", mm.group(1))
        # which section?
        pos = mm.start()
        sec = [x.group(1) for x in ms if x.start() < pos][-1]
        print(f"  {env:8} {lab.group(1) if lab else '(unlabeled)':<15} in {sec}")
