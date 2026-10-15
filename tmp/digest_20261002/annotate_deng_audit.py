"""Annotate deng_feedback_wiring_audit.md with Ben's 2026-10-03 corrections."""
from pathlib import Path

p = Path(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20260924\deng_feedback_wiring_audit.md")
t = p.read_text(encoding="utf8")

annotation = """
---

## CORRECTION (Ben Bolen, 2026-10-03 — screenshot verification)

**The claim above that "nothing sensory touches an RG half-center" is WRONG.**
Ben opened the original Deng model in AnimatLab and confirmed that afferent
connections DO reach the RG/PF half-centers. The data-mining script below
either looked at the wrong variant (the contact-driven BilateralRG transform
rather than the pristine W2L original) or misclassified the afferent-to-HC
connections as "laminated IN drive" because of how AnimatLab groups the
connection labels.

**Implication for the current stack:** the `ii/ia_f_central` knobs that
connect afferents directly to the RG half-centers in `build_network.py`
are CONSISTENT with the Deng original, not a deviation from it. The
WIRING_RULINGS_20261003.md conflict table has been corrected to match.

**Ben's annotation on the feedback character:** "The feedback shown here is
worked for their model, although in some respects it is opposite of what
other literature claims for the type Ia and type II feedback. Namely, type
Ia and II feedbacks are usually described (e.g. Shinohara) as agonistically
excitatory, whereas here they are antagonistically inhibitory."

**Ben's doctrine on gains:** "Any gains I put in the models were based on
certain models in literature and should only be taken as a potential
starting point for optimizing/tuning, not the final word."

**Ben's forward-looking note on the connectome editor:** the neural
subsystem types visible in the Deng original (afferent-to-HC projections,
layered INs with mixed signs, the contact-driven architecture) are the
kinds of subsystem he would eventually like to see as named primitives in
the connectome block editor, but he suspects abstracting at that level may
be too large an ask for current AI tooling.
"""

t = t.rstrip() + "\n" + annotation
p.write_text(t, encoding="utf8", newline="\n")
print("annotated")
