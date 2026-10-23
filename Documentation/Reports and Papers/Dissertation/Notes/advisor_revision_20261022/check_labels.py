import re
import sys

base = open(sys.argv[1], encoding="utf-8-sig").read()
new = open(sys.argv[2], encoding="utf-8-sig").read()

bl = set(re.findall(r"\\label\{([^}]+)\}", base))
nl = set(re.findall(r"\\label\{([^}]+)\}", new))
missing = bl - nl
added = nl - bl
print("labels in base:", len(bl), " in new:", len(nl))
print("MISSING from new:", sorted(missing) if missing else "none")
print("ADDED in new:", sorted(added) if added else "none")

# refs used in new that don't exist in new
br = set(re.findall(r"\\ref\{([^}]+)\}", new))
undef = br - nl
print("undefined refs in new:", sorted(undef) if undef else "none")

# external refs: refs used elsewhere pointing at this chapter's labels
# (checked separately per chapter)
figs_b = re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", base)
figs_n = re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", new)
print("figures base:", len(figs_b), "new:", len(figs_n))
print("figure files missing in new:", sorted(set(figs_b) - set(figs_n)) or "none")
