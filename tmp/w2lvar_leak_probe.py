"""Which edge classes changed in w2lvar_flat after the spiking build was
inserted before the variant builds?"""
import io, json, sys
from collections import Counter
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

OLD = json.load(open(r"D:\Github\Bipedal_Robot\tmp\connectome_templates_before_20261003.json",
                     encoding="utf-8"))
NEW = json.load(open(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\connectome_templates.json",
                     encoding="utf-8"))


def census(spec):
    c = Counter()
    for e in spec["edges"]:
        c[(e["from"].split("_")[0], e["to"].split("_")[0], e["sign"],
           e.get("tag", ""))] += 1
    return c


a, b = census(OLD["w2lvar_flat"]), census(NEW["w2lvar_flat"])
print("classes only/more in NEW (polluted):")
for k in sorted(set(a) | set(b)):
    d = b.get(k, 0) - a.get(k, 0)
    if d:
        print("  %+5d  %s" % (d, k))
print("old total", sum(a.values()), "new total", sum(b.values()))
