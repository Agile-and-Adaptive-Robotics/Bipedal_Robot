import re, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

frag = open(r'D:\GitHub\Bipedal_Robot\tmp\paper_bridge\sections\20_methods.tex', encoding='utf-8').read()
bib = open(r'D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Knee_Torque_Test\ProofFinal\BolenFrontiers22.bib', encoding='utf-8').read()

bibkeys = set(re.findall(r'^@[a-zA-Z]+\{([^,\s]+),', bib, re.M))

# 1. citation keys in fragment
cites = re.findall(r'\\cite[tp]?\{([^}]*)\}', frag)
keys = []
for c in cites:
    keys += [k.strip() for k in c.split(',')]
missing = sorted(set(keys) - bibkeys)
print("citation calls:", len(cites), "| unique keys:", len(set(keys)))
print("MISSING KEYS:", missing if missing else "NONE - all present")

# 2. banned AI-slop sweep (case-insensitive)
banned = [r'delve', r'it is important to note', r'plays? a crucial role',
          r'\blandscape\b', r'moreover', r'furthermore', r'additionally',
          r'\bplant\b', r'crucial', r'\brevolutioniz']
body = '\n'.join(l for l in frag.splitlines() if not l.lstrip().startswith('%'))
for pat in banned:
    hits = [(body[:m.start()].count('\n')+1, body[max(0,m.start()-30):m.end()+30].replace('\n',' '))
            for m in re.finditer(pat, body, re.I)]
    if hits: print(f"BANNED '{pat}':", hits)
emdash = frag.count('\u2014')
print("em-dashes:", emdash)

# 3. labels defined / refs used
labels = re.findall(r'\\label\{([^}]*)\}', frag)
refs = re.findall(r'\\(?:eq)?ref\{([^}]*)\}', frag)
print("labels:", sorted(set(labels)))
print("internal refs unresolved in-fragment (expected: figure/results labels):",
      sorted(set(r for r in refs if r not in labels)))

# 4. environment balance
for env in ['equation','cases','aligned','bmatrix','figure','center']:
    b = len(re.findall(r'\\begin\{'+env+r'\}', frag)); e = len(re.findall(r'\\end\{'+env+r'\}', frag))
    print(f"env {env}: begin={b} end={e} {'OK' if b==e else 'MISMATCH'}")
print("brace balance:", frag.count('{') - frag.count('}'))

# 5. includegraphics targets
print("graphics:", re.findall(r'\\includegraphics\[[^\]]*\]\{([^}]*)\}', frag))

# 6. prose word count (strip comments, commands, math, captions of floats kept? count body text only per style card: excl. floats/equations/captions)
lines = [l for l in frag.splitlines() if not l.lstrip().startswith('%')]
text = '\n'.join(lines)
# drop equations
text = re.sub(r'\\begin\{equation\}.*?\\end\{equation\}', ' ', text, flags=re.S)
# drop figure floats entirely
text = re.sub(r'\\begin\{figure\}.*?\\end\{figure\}', ' ', text, flags=re.S)
text = re.sub(r'\\(sub)*section\{[^}]*\}', ' ', text)
text = re.sub(r'\\[a-zA-Z]+\*?(\[[^\]]*\])?', ' ', text)
text = re.sub(r'[{}~$&\\]', ' ', text)
words = [w for w in re.split(r'\s+', text) if re.search(r'[A-Za-z]', w)]
print("prose words (excl. equations, figure floats, headers):", len(words))
