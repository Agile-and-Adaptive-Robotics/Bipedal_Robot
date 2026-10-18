"""Quantitative style analysis of the published Actuators paper (bolen_2026).

Source: Documentation/Reports and Papers/Isometric_Force_Test/Final/Final_proofread/
        manuscript.v5/actuators-4207821-edited.tex  (final proofread = published text)
Run with: D:/Anaconda/python.exe tmp/paper_bridge/mining/_analyze_actuators_style.py
"""
import re, sys

SRC = r"D:/GitHub/Bipedal_Robot/Documentation/Reports and Papers/Isometric_Force_Test/Final/Final_proofread/manuscript.v5/actuators-4207821-edited.tex"

text = open(SRC, encoding="utf-8", errors="replace").read()

# --- isolate body (between \begin{document} and \end{thebibliography}) ---
body_start = text.index(r"\begin{document}")
body_end = text.index(r"\begin{thebibliography}")
body = text[body_start:body_end]

def strip_comments(s):
    return "\n".join(l.split("%", 1)[0] if not l.lstrip().startswith("%") else "" for l in s.splitlines())

def wordcount(s):
    # remove tex commands / braces / math, then count words
    s = strip_comments(s)
    s = re.sub(r"\\begin\{(equation|align|adjustwidth|figure|table|tabularx|center|cases|aligned)\*?\}.*?\\end\{\1\*?\}", " MATHBLOCK ", s, flags=re.S)
    s = re.sub(r"\\(?:sub)?section\*?\{([^}]*)\}", r" \1 ", s)
    s = re.sub(r"\\(?:caption|abstract|Title|keyword)\{", " ", s)
    s = re.sub(r"\\[a-zA-Z]+(\[[^\]]*\])?(\{[^{}]*\})*", " ", s)
    s = re.sub(r"[{}\\$&~]", " ", s)
    return len(s.split())

# --- split by \section ---
sec_pat = re.compile(r"\\section\{([^}]*)\}")
matches = list(sec_pat.finditer(body))
print("=== SECTION SKELETON (order, word counts excl. floats/math) ===")
sections = []
for i, m in enumerate(matches):
    name = m.group(1)
    seg = body[m.end(): matches[i+1].start() if i+1 < len(matches) else len(body)]
    # front matter between begin{document} and first section
    sections.append((name, seg))
    print(f"  {name:<45} ~{wordcount(seg):>5} words  (tex lines {m.start()})")
front = body[:matches[0].start()]
print(f"  {'<front matter: title/abstract/highlights>':<45} ~{wordcount(front):>5} words")

# abstract specifically
am = re.search(r"\\abstract\{(.*)\}", text[body_start:], flags=re.S)
if am:
    aw = re.sub(r"\\[a-zA-Z]+", " ", am.group(1))
    aw = re.sub(r"[{}\\]", " ", aw)
    print(f"\nABSTRACT word count: {len(aw.split())}")

# --- citations ---
cites = re.findall(r"\\cite[tp]?\{([^}]*)\}", strip_comments(body))
keys = [k.strip() for c in cites for k in c.split(",")]
print(f"\n=== CITATIONS ===\n  \\citep/\\citet calls in body: {len(cites)}; total key instances: {len(keys)}; unique keys: {len(set(keys))}")
narrative = len(re.findall(r"\\citet\{", strip_comments(body)))
parenth = len(re.findall(r"\\citep\{", strip_comments(body)))
bare = len(re.findall(r"(?<!\\cite[tp])\\cite\{", strip_comments(body)))
print(f"  narrative \\citet: {narrative}; parenthetical \\citep: {parenth}; bare \\cite: {bare}")
bodyw = wordcount(body)
print(f"  body word count (excl floats/math): {bodyw}; citation calls per 100 words: {100*len(cites)/bodyw:.2f}")

# citation density per section
print("\n  per-section citation calls:")
for name, seg in sections:
    c = re.findall(r"\\cite[tp]?\{", strip_comments(seg))
    print(f"    {name:<45} {len(c)}")

# --- equations ---
eqs = re.findall(r"\\label\{(eq:[^}]*)\}", strip_comments(body))
print(f"\n=== EQUATIONS ===\n  labeled equations: {len(eqs)}")
print("  labels:", ", ".join(eqs))
eq_envs = len(re.findall(r"\\begin\{equation\}", strip_comments(body))) + len(re.findall(r"\\begin\{align\}", strip_comments(body)))
print(f"  equation/align environments: {eq_envs}")

# --- figures & tables ---
figs = re.findall(r"\\label\{(fig:[^}]*)\}", strip_comments(body))
tabs = re.findall(r"\\label\{(tab:[^}]*)\}", strip_comments(body))
print(f"\n=== FLOATS ===\n  figures: {len(figs)} -> {figs}")
print(f"  tables: {len(tabs)} -> {tabs}")
figrefs = re.findall(r"Figure~\\ref\{([^}]*)\}", strip_comments(body))
tabrefs = re.findall(r"Table~\\ref\{([^}]*)\}", strip_comments(body))
print(f"  in-text 'Figure~\\ref' mentions: {len(figrefs)}; 'Table~\\ref': {len(tabrefs)}")
eqrefs = re.findall(r"\(?\\ref\{(eq:[^}]*)\}\)?", strip_comments(body))
print(f"  in-text equation \\ref mentions: {len(eqrefs)}")

# --- paragraph counts per section ---
print("\n=== PARAGRAPH COUNTS ===")
for name, seg in sections:
    seg2 = strip_comments(seg)
    paras = [p for p in re.split(r"\n\s*\n", seg2) if p.strip() and not p.strip().startswith("\\begin") and "%" not in p[:2]]
    print(f"  {name:<45} ~{len(paras)} blocks")

# --- subsections ---
print("\n=== SUBSECTIONS ===")
for name, seg in sections:
    subs = re.findall(r"\\subsection\{([^}]*)\}", strip_comments(seg))
    if subs:
        print(f"  {name}:")
        for s in subs:
            print(f"     - {s}")
