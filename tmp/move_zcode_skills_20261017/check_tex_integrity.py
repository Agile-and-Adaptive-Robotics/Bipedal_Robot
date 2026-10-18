import re, glob, io, os
base = r"D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal"
files = sorted(glob.glob(os.path.join(base, "chapters", "*.tex"))) + [os.path.join(base, "main.tex")]
labels, refs = set(), []
env_stack_problems = []
for f in files:
    txt = io.open(f, encoding="latin-1").read()
    # strip comments (unescaped % to EOL)
    txt_nc = re.sub(r'(?<!\\)%.*', '', txt)
    for m in re.finditer(r'\\label\{([^}]+)\}', txt_nc):
        labels.add(m.group(1))
    for m in re.finditer(r'\\(?:ref|eqref|pageref|autoref)\{([^}]+)\}', txt_nc):
        for lbl in m.group(1).split(','):
            refs.append((os.path.basename(f), lbl.strip()))
    # brace balance
    ob, cb = txt_nc.count('{'), txt_nc.count('}')
    if ob != cb:
        print(f"BRACE-IMBALANCE {os.path.basename(f)}: {{={ob} }}={cb}")
    # env matching
    begins = re.findall(r'\\begin\{([^}]+)\}', txt_nc)
    ends = re.findall(r'\\end\{([^}]+)\}', txt_nc)
    if sorted(begins) != sorted(ends):
        from collections import Counter
        cb_, ce_ = Counter(begins), Counter(ends)
        env_stack_problems.append((os.path.basename(f), cb_ - ce_, ce_ - cb_))
missing = {}
for src, lbl in refs:
    if lbl not in labels:
        missing.setdefault(src, []).append(lbl)
if missing:
    print("MISSING-REFS:")
    for src, ls in missing.items():
        print(" ", src, "->", sorted(set(ls)))
else:
    print("REFS-OK (every \\ref resolves to a \\label)")
for f, b, e in env_stack_problems:
    print("ENV-MISMATCH", f, "unmatched-begin:", dict(b), "unmatched-end:", dict(e))
if not env_stack_problems:
    print("ENVIRONMENTS-OK")
# duplicate labels
all_lab = []
for f in files:
    txt = io.open(f, encoding="latin-1").read()
    txt_nc = re.sub(r'(?<!\\)%.*', '', txt)
    all_lab += re.findall(r'\\label\{([^}]+)\}', txt_nc)
dups = {l for l in all_lab if all_lab.count(l) > 1}
print("DUP-LABELS:", sorted(dups) if dups else "none")
