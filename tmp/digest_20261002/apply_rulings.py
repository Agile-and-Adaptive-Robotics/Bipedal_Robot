"""Apply Ben's 2026-10-03 late-night rulings: flexor Xi3 = 0.1583 canonical (placeholder
resolved), extensor 09-25 margin restated as +5.02% everywhere."""
from pathlib import Path

PF = Path(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal\chapters")

def swap(fname, old, new):
    p = PF / fname
    t = p.read_text(encoding="utf8")
    assert old in t, f"anchor missing in {fname}: {old[:80]}"
    assert t.count(old) == 1, f"anchor not unique in {fname}"
    p.write_text(t.replace(old, new), encoding="utf8", newline="\n")
    print(f"{fname}: applied")

r = "30-results.tex"
swap(r,
     "together with $\\chi_{3} = 0.1583$ loaded from the extensor refit.",
     "together with $\\chi_{3} = 0.1583$ loaded from the extensor refit and adopted as the flexor value of record; the preflight value of 0.6209 that appears once in the 2026-09-20 campaign log is superseded.")
swap(r,
     "% FLAG-BEN 2026-10-03: flexor Xi3 canonicality unresolved. The 2026-09-20 campaign log (Dig_out/Opt_flxr77_20260920.log) uses Xi3 = 0.6209 in the flexor preflight, while every final flexor result file, including the pulley file, carries Xi3 = 0.1583 (the extensor value). No source states which is canonical for the flexor or whether Xi3 is inert in the flexor objective. TODO(Ben): adjudicate, then adjust the visible placeholder and the preceding sentence if needed.",
     "% [ZCODE 2026-10-03: flexor Xi3 adjudicated by Ben 2026-10-03: 0.1583 (the value carried by every final flexor result file, Bifemsh incl. pulley) is canonical; the 0.6209 in the Opt_flxr77_20260920.log preflight is superseded.]")
t = (PF / r).read_text(encoding="utf8")
old_vis = "\n\\noindent\\textit{[FLAG-BEN, FLEXOR $\\chi_{3}$ VALUE: the 2026-09-20 campaign log quotes 0.6209 in the flexor preflight while every final flexor result file carries 0.1583; which value is canonical for the flexor, and whether $\\chi_{3}$ acts in the flexor objective at all, is unresolved. TODO(Ben): adjudicate and restate.]}\n"
assert old_vis in t, "visible placeholder block not found verbatim"
(PF / r).write_text(t.replace(old_vis, "\n"), encoding="utf8", newline="\n")
print("30-results: visible Xi3 placeholder removed")
swap(r, "improved to $+5.03$\\,\\% by the 2026-09-25 re-run of record",
        "improved to $+5.02$\\,\\% by the 2026-09-25 re-run of record")
swap("94-AppendixC.tex", "and $+5.03$\\,\\% (extensor, 2026-09-25)",
                        "and $+5.02$\\,\\% (extensor, 2026-09-25)")
for f in ("30-results.tex", "40-discussion.tex", "94-AppendixC.tex", "20-methods.tex"):
    t = (PF / f).read_text(encoding="utf8")
    if "5.03" in t:
        print(f"WARNING: '5.03' still present in {f}")
print("done")
