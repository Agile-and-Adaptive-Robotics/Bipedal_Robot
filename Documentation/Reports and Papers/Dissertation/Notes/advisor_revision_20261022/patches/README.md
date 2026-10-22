# Patches — two formats

Everything diffs against the **pristine baseline** in `../baseline/chapters/`
(the dated copy of the live source captured 2026-10-22, advisor tracked edits
already applied). Anchor identity: the `baseline/` folder IS the anchor.

## Full-chapter patches (convenience)

`10-introduction.patch`, `15-background.patch`, `30-torque.patch`,
`50-database.patch`, `60-sns.patch`, `70-overall.patch`, `80-conclusion.patch`,
`94-AppendixC.patch` — standard context diffs. Apply from the REPO ROOT (`D:\Github\Bipedal_Robot`), quoting the whole path because of the spaces:

    git apply "Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/patches/50-database.patch"

(The header paths assume the repo root; alternatively copy the chapter file over
from `../proposals/` — the proposals are the same content.)

## Discrete sub-patches (apply any subset)

`10-introduction/`, `15-background/`, `30-torque/` — zero-context patches, each
covering ONE logical fix family, each verified to apply cleanly to the
unmodified baseline:

- 10-introduction: 01-gaps-and-objectives, 02-project-connection-and-person,
  03-organization-section-fixes (INTRO-31 + terminology)
- 15-background: 01-gait2392-naming-kinematics, 02-masters-paragraph,
  03-compneuro-cpg-move, 04-remove-graphical-abstract (needs your sign-off),
  05-first-person-and-opener
- 30-torque: 01-physical-validation-sweep, 02-publication-note

Apply from the REPO ROOT with (zero-context diffs need the flag):

    git apply --unidiff-zero "Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/patches/30-torque/01-physical-validation-sweep.patch"

Every sub-patch of a chapter is a disjoint hunk of the same full-chapter diff,
so applying ALL of them reproduces the proposal exactly (content-identical;
line endings follow your git autocrlf). SNS and SADb are delivered only as
whole-chapter patches because their advisor-mandated restructuring moves
material across sections and cannot be split into independent hunks.

`70-overall`, `80-conclusion`, `94-AppendixC` are minimal first-person sweeps
(2026-10-22 ruling: "I" replaces "the author" everywhere).
