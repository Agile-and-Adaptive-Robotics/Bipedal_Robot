# Checker for the abstract+intro fragment (run with D:/Anaconda/python.exe).
# 1) Every \\citet/\\citep key in the fragment must exist in ProofFinal/BolenFrontiers22.bib.
# 2) Banned-pattern scan per BRIEF editorial contract + style contract.
# 3) Word/sentence counts vs plan budgets (abstract ~150 w / 5 sentences; intro ~800 w).
import re, sys

FRAG = r"D:\GitHub\Bipedal_Robot\tmp\paper_bridge\sections\10_abstract_intro.tex"
BIB = r"D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Knee_Torque_Test\ProofFinal\BolenFrontiers22.bib"

frag = open(FRAG, encoding="utf-8").read()
bib = open(BIB, encoding="utf-8", errors="replace").read()

bib_keys = set(re.findall(r"@\w+\s*\{\s*([^,\s]+)\s*,", bib))

# strip comment lines (those starting with %) for citation scanning
active = "\n".join(l for l in frag.splitlines() if not l.lstrip().startswith("%"))

cite_calls = re.findall(r"\\cite[tp]\{([^}]*)\}", active)
keys_used = []
for c in cite_calls:
    for k in c.split(","):
        k = k.strip()
        if k:
            keys_used.append(k)

missing = [k for k in keys_used if k not in bib_keys]
print(f"citation calls: {len(cite_calls)}; keys used: {len(keys_used)}; unique: {len(set(keys_used))}")
print(f"MISSING KEYS: {missing if missing else 'NONE'}")
dupe = sorted({k for k in keys_used if keys_used.count(k) > 1})
print(f"repeated keys (fine, informational): {dupe if dupe else 'none'}")

# ---- banned patterns ----
banned = [
    r"\bdelve\b", r"it is important to note", r"plays? a crucial role",
    r"\blandscape\b", r"\bmoreover\b", r"\bfurthermore\b", r"\badditionally\b",
    r"\bplant\b", r"\bplants\b",
    r"\bbehaviour\b", r"\bfibre\b", r"\bmodelled\b", r"\bcolour\b", r"\bcentre\b",
]
hits = []
for pat in banned:
    for m in re.finditer(pat, active, re.I):
        hits.append((pat, active[max(0, m.start()-40):m.end()+40].replace("\n", " ")))
print(f"banned-pattern hits: {hits if hits else 'NONE'}")

emdash = active.count("---")
print(f"em-dashes (---): {emdash}")
dbl = len(re.findall(r"[a-z]  +[a-z]", active))
print(f"double spaces between words: {dbl}")

# bare units: number+unit without siunitx (e.g. "10mm", "10 mm", "620 kPa" outside \qty)
bare = re.findall(r"(?<!\\qty\{)\b\d+(?:\.\d+)?\s*(?:mm|cm|kPa|N\bm|N \* m)\b", active)
print(f"bare-unit candidates (check each is inside \\qty): {bare if bare else 'NONE'}")

# ---- counts ----
# abstract body = between \section{} (first, inside abstract) and \tiny
m_abs = re.search(r"\\begin\{abstract\}.*?\\section\{\}\s*(.*?)\s*\\tiny", frag, re.S)
abs_text = m_abs.group(1)
# strip latex commands lightly for word counting
def words(t):
    t = re.sub(r"\\qty\{([^}]*)\}\{[^}]*\}", r"\1", t)
    t = re.sub(r"\\qtylist\{([^}]*)\}\{[^}]*\}", lambda mm: mm.group(1).replace(";", " and "), t)
    t = re.sub(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?(\{[^{}]*\})?", " ", t)
    t = re.sub(r"[${}~\\]", " ", t)
    return t.split()

aw = words(abs_text)
sents = [s for s in re.split(r"(?<=[.!?])\s+", abs_text.strip()) if s]
print(f"ABSTRACT: {len(aw)} words, {len(sents)} sentences")
abs_cites = len(re.findall(r"\\cite[tp]", abs_text))
print(f"  abstract citations: {abs_cites}")

# intro = after \section{Introduction}
m_int = re.search(r"\\section\{Introduction\}\s*(.*)$", frag, re.S)
int_text = m_int.group(1)
iw = words(int_text)
paras = [p for p in re.split(r"\n\s*\n", int_text.strip()) if p.strip()]
print(f"INTRO: {len(iw)} words, {len(paras)} paragraphs")
cite_pat = re.compile(r"\\cite[tp]")
for i, p in enumerate(paras, 1):
    pw = words(p)
    print(f"  P{i}: {len(pw)} words, {len(cite_pat.findall(p))} citation calls")
