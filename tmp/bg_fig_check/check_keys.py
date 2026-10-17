import re

ch = open(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\15-background.tex", encoding="utf-8").read()
bib = open(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal\thesis.bib", encoding="utf-8", errors="replace").read()

keys = set()
for g in re.findall(r"\\cite[tp]?(?:\[[^\]]*\])?\{([^}]*)\}", ch):
    keys.update(k.strip() for k in g.split(","))
have = set(re.findall(r"@\w+\{([^,]+),", bib))
miss = sorted(keys - have)
print("total distinct cited:", len(keys))
print("missing from thesis.bib:", miss)

# slop phrases
for p in ["delve", "crucial role", "important to note", "Moreover",
          "Furthermore", "Additionally", "landscape"]:
    hits = [i + 1 for i, line in enumerate(ch.splitlines()) if p.lower() in line.lower()]
    if hits:
        print("SLOP:", p, "lines", hits)

# em-dash occurrences
em = [i + 1 for i, line in enumerate(ch.splitlines()) if "---" in line]
print("em-dash lines:", em)
