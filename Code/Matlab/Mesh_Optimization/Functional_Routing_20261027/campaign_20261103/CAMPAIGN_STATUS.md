# Routing campaign paused

The source in this directory includes unfinished mechanics and packaging experiments. Optimization is paused at Ben's direction. The earlier README describes the initial campaign.

The tracked `audit_20261105/ROUTING_AUDIT_20261105.md` documents the repaired baseline. Eight of 54 cases passed its implemented constraints, and zero passed both those constraints and all sampled signed torque targets. Later local whole-leg packaging reports recorded only one of 27 rigid groups and two of 27 compliant groups clearing the modeled checks. The split-vasti redraw recorded zero of 29 rigid groups and one of 29 compliant groups clearing. Robot CAD registration remained unresolved, with approximately 87 mm knee and 89 mm ankle mismatch.

The ten `run_20261105_*` and `run_20261106_*` output directories remain local and are excluded from Git. They retain checkpoints, failed candidates, clearance reports, and figures. Source changes call the new `fr_*.m` helpers, so retain the source and configuration changes together when preserving this work.

Slide 57 in `Results/KompareXXVII/Slide57` uses a separate survey schematic. Its README documents the display transforms and limitations. Its generator requires the selected `erector_rerun/erector_compliant.json` and `erector_rerun/erector_rigid.json` records, which remain visible to Git.

Implementation checks and schematic exports establish their stated scopes. Complete physical reeving, torque continuity, assembled clearance, and robot CAD registration still require verification before a routing candidate can serve as a validated robot design.
