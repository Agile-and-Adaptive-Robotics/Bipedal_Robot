import hashlib, os

ROOT = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation"
zc = os.path.join(ROOT, "zip_extract", "chapters")
pc = os.path.join(ROOT, "ProofFinal", "chapters")

def h(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()

print("--- chapter diffs (zip vs local ProofFinal):")
for f in sorted(os.listdir(zc)):
    zp = os.path.join(zc, f)
    pp = os.path.join(pc, f)
    if not os.path.exists(pp):
        print(f, " MISSING locally")
    elif h(zp) != h(pp):
        print(f, " DIFFERS")
zf = os.path.join(ROOT, "zip_extract", "main.tex")
pf = os.path.join(ROOT, "ProofFinal", "main.tex")
print("main.tex:", "same" if h(zf) == h(pf) else "DIFFERS")
figs = os.path.join(ROOT, "zip_extract", "figs")
n = sum(len(fs) for _, _, fs in os.walk(figs))
print("zip figs file count:", n)
