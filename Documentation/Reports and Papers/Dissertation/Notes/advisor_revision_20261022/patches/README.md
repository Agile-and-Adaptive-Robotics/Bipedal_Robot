# Patches — how to apply them (corrected 2026-10-22 after review)

Everything diffs against the **pristine baseline** in `../baseline/chapters/`
(the dated copy of the live source captured 2026-10-22, advisor tracked edits
already applied). The `baseline/` folder IS the anchor.

## Before applying anything: the baseline check (do not skip)

These patches anchor to the 2026-10-22 snapshot. If the live chapter has moved
since (for example a newer Overleaf sync), applying a stale patch can silently
clobber newer work. Pre-flight, from `ProofFinal/`:

    git --no-pager diff --no-index "<...advisor_revision_20261022>/baseline/chapters/60-sns.tex" "chapters/60-sns.tex"

- **No output** → the live chapter still equals the anchor; apply normally.
- **Differences** → do NOT apply. Regenerate the patch against the newer live
  file first (diff the newer live file against the proposal), or re-derive the
  edits by hand. (`git` = GitHub Desktop's bundled git on this machine.)

## Discrete sub-patches (apply any subset)

`10-introduction/`, `15-background/`, `30-torque/` — zero-context patches, one
logical fix family per file, each verified to apply cleanly to the pristine
baseline. Their headers target `chapters/<file>.tex`, so **apply from the
`ProofFinal/` directory** (that is where `chapters/` lives — not the repo root):

    cd "D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal"
    git apply --unidiff-zero "...\advisor_revision_20261022\patches\30-torque\01-physical-validation-sweep.patch"

(zero-context diffs require `--unidiff-zero`.) Every sub-patch of a chapter is
a disjoint hunk of the same chapter diff; applying all of them reproduces the
proposal content exactly (line endings follow your git autocrlf).

## Whole-chapter patches (50, 60, and the 70/80/94 person sweeps)

`50-database.patch`, `60-sns.patch`, `70-overall.patch`, `80-conclusion.patch`,
`94-AppendixC.patch` are **review diffs** (baseline vs. proposal, generated
inside the advisor_revision folder). Do not `git apply` them in place — their
a/b paths point at the snapshot folders, and applying them there would edit the
anchor. To adopt one: read the .patch, then copy the corresponding file from
`../proposals/` over the live chapter in `ProofFinal/chapters/` (same content,
one step). SNS and SADb cannot be split into discrete hunks because their
advisor-mandated restructuring moves material across sections.

## Inventory

- 10-introduction: 01-gaps-and-objectives, 02-project-connection-and-person,
  03-organization-section-fixes (INTRO-31 + terminology)
- 15-background: 01-gait2392-naming-kinematics, 02-masters-paragraph,
  03-compneuro-cpg-move, 04-remove-graphical-abstract (needs Ben's sign-off),
  05-first-person-and-opener
- 30-torque: 01-physical-validation-sweep, 02-publication-note
- Whole-chapter: 50-database, 60-sns, 70-overall, 80-conclusion, 94-AppendixC
  (the last three are the minimal first-person sweeps, 2026-10-22 ruling)
