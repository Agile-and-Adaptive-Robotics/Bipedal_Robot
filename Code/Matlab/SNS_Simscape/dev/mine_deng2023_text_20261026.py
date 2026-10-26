# Mine the Deng 2023 fulltext (kaiyu2023_fulltext.txt, single-line dump)
# for the model specification: neuron/synapse values, pathways, muscles,
# and anything per-synapse. Output: structured findings + context quotes.
import re
from pathlib import Path

P = Path(r"D:\Github\Bipedal_Robot\Code\Matlab\SNS_Simscape"
         r"\dev\kaiyu2023_fulltext.txt")
s = P.read_text(encoding="utf-8", errors="replace")
print("chars:", len(s))
# normalize for context extraction
flat = re.sub(r"\s+", " ", s)

PATTERNS = [
    ("title", r".{0,80}(?:Title|TITLE).{0,200}"),
    ("conductance_uS", r"[^.]*(?:\u00b5S|uS|micro[Ss]iemens|conductance)[^.]*\."),
    ("nA_values", r"[^.]*\b\d+(?:\.\d+)?\s*nA[^.]*\."),
    ("synap", r"[^.]*[Ss]ynap[^.]*\."),
    ("muscles", r"[^.]*\b(?:ILPSO|GMAX|HAMS|BFSH|VAS|RF|GAS|SOL|TA|VA|BF)\b[^.]*\."),
    ("reflex", r"[^.]*[Rr]eflex[^.]*\."),
    ("timestep", r"[^.]*(?:time ?step|dt ?=|0\.\d+ ?ms)[^.]*\."),
    ("optimi", r"[^.]*[Oo]ptimi[^.]*\."),
    ("model_names", r"[^.]*\b(?:AnimatLab|Mujoco|MuJoCo|SIMULINK|Simulink)[^.]*\."),
    ("freq_period", r"[^.]*\b(?:Hz|burst|rhythm|alternat)[^.]*\."),
]
seen = set()
out = []
for name, pat in PATTERNS:
    ms = re.findall(pat, flat)
    out.append("### " + name + " (" + str(len(ms)) + " hits)")
    k = 0
    for m in ms:
        mm = m.strip()[:300]
        if mm in seen:
            continue
        seen.add(mm)
        k += 1
        if k <= 8:
            out.append("  - " + mm)
    out.append("")
Path(r"D:\Github\Bipedal_Robot\Code\Matlab\SNS_Simscape\dev"
     r"\deng2023_mine_findings.txt").write_text(
    "\n".join(out), encoding="utf-8")
print("\n".join(out[:60]))
