"""Correct the wrong 'Deng original: none' afferents claim in WIRING_RULINGS_20261003.md."""
from pathlib import Path

p = Path(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\WIRING_RULINGS_20261003.md")
t = p.read_text(encoding="utf8")

old_row = "| 1 | Afferents onto RG/PF half-centers | Deng original: none; Ben master rules + Shinohara eq10/11: direct; Shevtsova: no afferents | Ben's drawing + Shinohara (kept, F4) |"
new_row = "| 1 | Afferents onto RG/PF half-centers | Deng original: PRESENT (Ben, 2026-10-03 screenshot verification: afferents DO reach the RG/PF half-centers in the original model; the earlier 'none' was an audit error in reports_20260924/deng_feedback_wiring_audit.md); Ben master rules + Shinohara eq10/11: direct; Shevtsova: no afferents in her abstraction | Ben's drawing + Shinohara + Deng original all agree: direct afferents to HCs are correct (not a deviation) |"

assert old_row in t, "old row not found"
t = t.replace(old_row, new_row)

# Also check for the same wrong claim elsewhere in the doc
old_phrase = "Current-stack deviation: ii/ia_f_central knobs connect DIRECTLY to RG HCs"
if old_phrase in t:
    t = t.replace(old_phrase,
        "ii/ia_f_central knobs connect DIRECTLY to RG HCs (consistent with the Deng original, which does have afferents reaching the HCs per Ben's 2026-10-03 verification)")

p.write_text(t, encoding="utf8", newline="\n")
print("corrected")
