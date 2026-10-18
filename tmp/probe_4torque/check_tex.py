# Structural sanity check of the revised 30-torque.tex (read-only).
import io
import re
import collections

P = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\30-torque.tex"
t = io.open(P, encoding="utf-8").read()

begins = collections.Counter(re.findall(r"\\begin\{(\w+\*?)\}", t))
ends = collections.Counter(re.findall(r"\\end\{(\w+\*?)\}", t))
print("begins:", dict(begins))
print("ends:  ", dict(ends))
print("env MATCH" if begins == ends else "env MISMATCH")

labs = re.findall(r"\\label\{([^}]+)\}", t)
dups = [k for k, v in collections.Counter(labs).items() if v > 1]
print("labels:", len(labs), "dups:", dups)

refs = set(re.findall(r"\\ref\{([^}]+)\}", t))
missing = sorted(r for r in refs if r not in set(labs))
print("refs:", len(refs), "unresolved-in-file (may be cross-chapter):", missing)

# brace balance
depth = 0
for ch in t:
    if ch == "{":
        depth += 1
    elif ch == "}":
        depth -= 1
print("brace depth at EOF:", depth)
