import re
from pathlib import Path

t = Path(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\30-results.tex").read_text(encoding="utf8")
for m in re.finditer(r"FLAG-BEN", t):
    print(f"--- offset {m.start()} ---")
    print(t[max(0, m.start() - 350):m.start() + 650].replace("\n", " | "))
    print()
i = t.find("XiUsed")
while i > 0:
    print("=== XiUsed mention ===")
    print(t[max(0, i - 300):i + 500].replace("\n", " | "))
    i = t.find("XiUsed", i + 1)
