"""Survey the Shevtsova 2026 fulltext cache: keyword coverage for the
connectome-relevant sections."""
import io
import re
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
t = open("shevtsova_2026_fulltext.txt", encoding="utf-8").read()
print("chars:", len(t))
for kw in ("V1 ", "V2b", "V2a", "V0", "V3 ", "laminar", "conductance",
           "RG-E", "RG-F", "half-center", "NaP", "persistent",
           "locomotor", "SCI", "reorganization"):
    print(f"  {kw!r}: {len(re.findall(re.escape(kw), t))}")
m = re.search(r"(?i)abstract.{0,900}", t)
if m:
    print("\nABSTRACT-ISH:\n", m.group(0)[:800])
