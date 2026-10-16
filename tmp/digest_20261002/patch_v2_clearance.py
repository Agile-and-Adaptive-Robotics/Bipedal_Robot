"""Patch buildRoutingContext.m for Ben's 2026-10-03 clearance scheme:
- Bone: prevent direct centerline penetration only (no inflated BPA radius + margin).
- BPA-BPA: prevent direct centerline overlap only; close-but-not-touching is a penalty.
- The anatomy penalty (off-bone distance beyond d0) handles 'too far from bone'.
"""
import re
from pathlib import Path

p = Path(r"D:\Github\Bipedal_Robot\tmp\routing27_v2_src\buildRoutingContext.m")
t = p.read_text(encoding="utf8")

# 1. Bone clearance: from inflated radius + margin to centerline-only
old_bone = "ctx.boneRequiredClear = ctx.bpaRadius + boneMargin;"
new_bone = ("ctx.boneRequiredClear = 0.001;  % centerline penetration guard only (Ben, 2026-10-03:\n"
            "                                     % 'no 25 mm clearance requirements'; the anatomy\n"
            "                                     % penalty handles off-bone distance)")
assert old_bone in t
t = t.replace(old_bone, new_bone)

# 2. BPA-BPA spacing: from diameter + margin to centerline overlap guard only
old_spacing = "ctx.spacingRequired = (bpaDia_mm * 1e-3) + spacingMargin;"
new_spacing = ("ctx.spacingRequired = 0.001;  % centerline overlap guard only (Ben, 2026-10-03:\n"
               "                                       % 'constrain against direct penetration of\n"
               "                                       % bone or BPAs, but otherwise just penalize')\n"
               "ctx.spacingComfort = (bpaDia_mm * 1e-3) + spacingMargin;  % soft threshold for the\n"
               "                                                              % too-close penalty")
assert old_spacing in t
t = t.replace(old_spacing, new_spacing)

# 3. Keep the anatD0 threshold but reduce it to reflect 'closest to bone' preference
old_anat = "ctx.anatD0 = getOrDefault(opts, 'anatD0_mm', 20) * 1e-3;"
new_anat = "ctx.anatD0 = getOrDefault(opts, 'anatD0_mm', 15) * 1e-3;  % Ben: start closest to bone"
assert old_anat in t
t = t.replace(old_anat, new_anat)

p.write_text(t, encoding="utf8", newline="\n")
print("buildRoutingContext.m patched")
