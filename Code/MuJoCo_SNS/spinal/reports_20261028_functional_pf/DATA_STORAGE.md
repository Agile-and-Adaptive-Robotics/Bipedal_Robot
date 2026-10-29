# Local simulation evidence

Raw trial configurations/results, NPZ trajectories, logs, databases, and
generated source snapshots stay on easteregg2 under this directory on `D:`.
The local `.gitignore` keeps these byproducts out of subsequent commits.
The family-audit directory uses the same policy.

Code and concise Markdown reports remain visible to Git. The small
`routing_audit.json` and `routing_samples.npz` stay tracked because the
experimental controller loads them at runtime. The curated final summary,
training specification, and validation summary also remain visible.

183 previously tracked generated files, totaling 510,428,214 bytes, were
removed from the index with `git rm --cached`; every working copy was
verified to remain on disk at its original size. These appear as staged
removals in GitHub Desktop. No files were deleted, no commit was made, and
existing history was not rewritten.

The local preservation manifest is
`D:\GitHub\Bipedal_Robot\tmp\simulation_artifact_untracking_manifest.json`.
Earlier frontal-axis measurements are superseded by the corrected review
and its local evidence under `literature_training_review/`.
