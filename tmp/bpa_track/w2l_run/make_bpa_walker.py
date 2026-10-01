r"""BPA-ify the w2l walker body: swap the 12 Hill-type <muscle> actuators of
w2l_mjcf_fixed.xml for ZERO-GAIN <general> actuators on the SAME tendon
routes, so the force comes from Ben's BPA model (bpa_muscle.BPAMuscle via
bpa_mujoco.BPAMuscleSystem) instead of MuJoCo's Hill muscle.

This is the add_bpa_to_mjcf.py convention (gaintype fixed / gainprm 0 /
biastype none) applied to a body whose site routes ALREADY exist as
<tendon><spatial> entries - add_bpa_to_mjcf.add_bpa() would append DUPLICATE
routes, so the swap is done here instead. Everything else in the XML
(freejoint, per-joint damping, toe springs, spawn keyframe) is untouched.

Usage:  python make_bpa_walker.py [xml_in] [xml_out]
        (defaults: w2l_mjcf_fixed.xml -> w2l_mjcf_bpa.xml, same folder)
"""
from __future__ import annotations

import io
import os
import re
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

HERE = os.path.dirname(os.path.abspath(__file__))

MUSCLE_RE = re.compile(
    r'<muscle name="([^"]+)" tendon="([^"]+)".*?/>')
GENERAL = ('<general name="{n}" tendon="{t}" gaintype="fixed" '
           'gainprm="0 0 0" biastype="none" biasprm="0 0 0" '
           'ctrlrange="0 1"/>')


def make_bpa_walker(xml_in: str, xml_out: str) -> int:
    txt = open(xml_in, encoding="utf-8").read()
    names = []

    def repl(m):
        n, t = m.group(1), m.group(2)
        names.append(n)
        return GENERAL.format(n=n, t=t)

    out = MUSCLE_RE.sub(repl, txt)
    n_in = len(MUSCLE_RE.findall(txt))
    assert len(names) == n_in == 12, f"expected 12 muscles, found {n_in}"
    assert "<muscle " not in out, "unconverted <muscle> left"
    hdr = ("<!-- BPA VARIANT (make_bpa_walker.py): the 12 Hill-type <muscle>\n"
           "     actuators are replaced by zero-gain <general> actuators on the\n"
           "     SAME tendon routes; force comes from bpa_muscle.BPAMuscle via\n"
           "     bpa_mujoco.BPAMuscleSystem (qfrc_applied), pressure-driven\n"
           "     (activation in [0,1] -> 0..620 kPa). Body geometry, damping,\n"
           "     freejoint, toe springs and spawn keyframe are identical to\n"
           "     the input file. -->\n")
    open(xml_out, "w", encoding="utf-8", newline="\n").write(hdr + out)
    print(f"wrote {xml_out}: {len(names)} zero-gain general actuators:")
    print("  " + ", ".join(sorted(names)))
    return 0


if __name__ == "__main__":
    xml_in = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        HERE, "w2l_mjcf_fixed.xml")
    xml_out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
        HERE, "w2l_mjcf_bpa.xml")
    sys.exit(make_bpa_walker(xml_in, xml_out))
