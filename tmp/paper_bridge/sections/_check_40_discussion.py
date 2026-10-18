"""QC checks for 40_discussion_conclusion.tex (run with D:/Anaconda/python.exe)."""
import re, sys, io

FRAG = r"D:\GitHub\Bipedal_Robot\tmp\paper_bridge\sections\40_discussion_conclusion.tex"
BIB = r"D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Knee_Torque_Test\ProofFinal\BolenFrontiers22.bib"
SIBLINGS = [
    r"D:\GitHub\Bipedal_Robot\tmp\paper_bridge\sections\20_methods.tex",
    r"D:\GitHub\Bipedal_Robot\tmp\paper_bridge\sections\30_results.tex",
    r"D:\GitHub\Bipedal_Robot\tmp\paper_bridge\sections\10_abstract_intro.tex",
    r"D:\GitHub\Bipedal_Robot\tmp\paper_bridge\sections\tables_main.tex",
]
SUPP = r"D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Knee_Torque_Test\ProofFinal\frontiers_SupplementaryMaterial.tex"

text = io.open(FRAG, encoding="utf-8").read()
lines = text.splitlines()

# strip comment lines (start with % after optional whitespace)
active = [l for l in lines if not l.lstrip().startswith("%")]

# ---- 1. per-section prose word counts (exclude floats: figure/table envs, captions) ----
def strip_env(s, env):
    return re.sub(r"\\begin\{%s\}.*?\\end\{%s\}" % (env, env), " ", s, flags=re.S)

prose = "\n".join(active)
for env in ("figure", "table", "equation"):
    prose = strip_env(prose, env)
prose = re.sub(r"\\caption\{.*?\}\}?", " ", prose, flags=re.S)
prose = re.sub(r"%.*", " ", prose)
prose = re.sub(r"\\(label|ref|eqref)\{[^}]*\}", " REF ", prose)
prose = re.sub(r"\\cite[tp]?\{[^}]*\}", " CITE ", prose)
prose = re.sub(r"\\(qty|num|qtylist|ang|si|unit)\{[^{}]*\}", " Q ", prose)
prose = re.sub(r"\$\$?.*?\$\$?", " MATH ", prose, flags=re.S)
prose = re.sub(r"\\[a-zA-Z]+", " ", prose)
prose = re.sub(r"[{}~]", " ", prose)

sections = re.split(r"(\\section\{[^}]*\})", "\n".join(active))
# simpler: walk active lines, bucket by \section/\subsection markers
buckets, cur = [], ("(front)", [])
for l in active:
    m = re.match(r"\\(sub)?section\{([^}]*)\}", l.strip())
    if m:
        if cur[1]:
            buckets.append(cur)
        cur = (m.group(2), [])
    else:
        cur[1].append(l)
buckets.append(cur)

def wc(ls):
    s = "\n".join(ls)
    for env in ("figure", "table", "equation"):
        s = strip_env(s, env)
    s = re.sub(r"%.*", " ", s)
    s = re.sub(r"\\(label|ref|eqref)\{[^}]*\}", " ", s)
    s = re.sub(r"\\cite[tp]?\{[^}]*\}", " ", s)
    s = re.sub(r"\\(qty|qtylist|ang|num|si|unit)\{[^{}]*\}(\{[^{}]*\})?", " N ", s)
    s = re.sub(r"\$\$?.*?\$\$?", " M ", s, flags=re.S)
    s = re.sub(r"\\[a-zA-Z]+", " ", s)
    s = re.sub(r"[{}~]", " ", s)
    return len(s.split())

print("== section prose word counts ==")
tot = 0
disc = 0
for name, ls in buckets:
    n = wc(ls)
    tot += n
    print(f"  {name}: {n}")
    if name != "Conclusion":
        disc += n
print(f"  TOTAL prose: {tot}  (Discussion-part {disc}, Conclusion {tot-disc})")

# ---- 2. banned phrases / style scan ----
print("\n== banned/style scan ==")
banned = [r"\bdelve\b", r"it is important to note", r"plays? a crucial role",
          r"\blandscape\b", r"\bmoreover\b", r"\bfurthermore\b", r"\badditionally\b",
          r"\bplant\b", r"\bcrucial\b"]
for pat in banned:
    hits = [i+1 for i, l in enumerate(lines) if re.search(pat, l, re.I)]
    print(f"  {'OK  ' if not hits else 'HIT '} {pat!r}: lines {hits}")
ndash = text.count("---")
print(f"  em-dashes (---): {ndash}")

# ---- 3. citation keys vs bib ----
print("\n== citations vs bib ==")
bib = io.open(BIB, encoding="utf-8", errors="replace").read()
bibkeys = set(re.findall(r"^@\w+\{([^,\s]+),", bib, re.M))
cited = re.findall(r"\\cite[tp]?\{([^}]*)\}", "\n".join(active))
keys = [k.strip() for c in cited for k in c.split(",")]
missing = [k for k in keys if k not in bibkeys]
print(f"  citation calls: {len(keys)}; unique: {len(set(keys))}; missing from bib: {missing or 'NONE'}")
print("  unique keys:", sorted(set(keys)))

# ---- 4. refs vs labels (mine + siblings + supp) ----
print("\n== \\ref resolution ==")
labels = set(re.findall(r"\\label\{([^}]+)\}", text))
for p in SIBLINGS + [SUPP]:
    s = io.open(p, encoding="utf-8", errors="replace").read()
    labels |= set(re.findall(r"\\label\{([^}]+)\}", s))
refs = set(re.findall(r"\\(?:ref|eqref)\{([^}]+)\}", "\n".join(active)))
unres = sorted(r for r in refs if r not in labels)
print(f"  refs used: {sorted(refs)}")
print(f"  UNRESOLVED: {unres or 'NONE'}")
mylabels = set(re.findall(r"\\label\{([^}]+)\}", text))
print(f"  my labels: {sorted(mylabels)}")
dup = sorted(mylabels & (labels - mylabels))
print(f"  collisions with siblings/supp: {dup or 'NONE'}")

# ---- 5. diameter usage (mathabx dependency) ----
print("\n\\diameter count:", len(re.findall(r"\\diameter", "\n".join(active))))

# ---- 6. bare-unit spot scan (cm/mm/kPa without siunitx, outside math) ----
bare = []
for i, l in enumerate(active):
    for m in re.finditer(r"(?<![\w\\])(-?\d+(?:\.\d+)?)\s*(cm|mm|kPa|N)\b(?!\\)", l):
        frag = l[max(0, m.start()-30):m.end()+10]
        if "$" not in frag.split(m.group(0))[0][-40:]:
            bare.append((i+1, m.group(0)))
print("bare-unit suspects (verify manually):", bare or "NONE")
