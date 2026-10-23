import re
import sys
import collections

log = open(sys.argv[1], encoding="latin-1", errors="replace").read()
counts = collections.Counter()
examples = collections.defaultdict(list)
depth_stack = []
# crude: track open parens with file-like names
cur = "?"
i = 0
pat = re.compile(r"\(([^\s()]+\.(?:tex|cls|sty|cfg|def|fd|clo|bbl))|Overfull \\[hv]box \(([\d.]+)pt[^\n]*")
for m in pat.finditer(log):
    if m.group(1):
        depth_stack.append(m.group(1))
        cur = m.group(1)
    elif m.group(2):
        key = cur.split("\\")[-1].split("/")[-1]
        counts[key] += 1
        if float(m.group(2)) >= 20 and len(examples[key]) < 6:
            examples[key].append(m.group(0)[:100])
for k, v in counts.most_common():
    print(f"{v:>4}  {k}")
print("--- examples >=20pt ---")
for k in examples:
    for e in examples[k]:
        print(k, "::", e)
