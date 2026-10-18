# Gate-1 helper: 10-word n-gram overlap between 60-sns intro paragraph (SNS community)
# and 15-background sec:bg_sns paragraphs (case/punct-normalized), per rubric 1b (own text).
import re

sns = open(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\60-sns.tex", encoding="utf-8").read()
bg = open(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\15-background.tex", encoding="utf-8").read()

def norm(t):
    t = re.sub(r"%.*", " ", t)
    t = re.sub(r"\\cite[tp]?\{[^}]*\}", " ", t)
    t = re.sub(r"\\[a-zA-Z]+", " ", t)
    t = re.sub(r"[^A-Za-z0-9\s-]", " ", t)
    return re.sub(r"\s+", " ", t.lower()).split()

def grams(w, n):
    return {" ".join(w[i:i+n]) for i in range(len(w)-n+1)}

a, b = norm(sns), norm(bg)
for n in (7, 10):
    ov = grams(a, n) & grams(b, n)
    print(f"{n}-gram overlaps between 60-sns and 15-background: {len(ov)}")
    for s in sorted(ov)[:20]:
        print("   ::", s)
