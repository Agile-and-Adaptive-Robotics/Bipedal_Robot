# EB475WS4 slide-deck pickup

## Current state

On November 2, GitHub branch `KneeTestSetup_BenBo_stw` was checked directly: its head is `338da208032e5a59ea8dfe934a93b0e949e65594` (`slide v21`). The subsequent plot-view correction and this handoff are local changes awaiting Ben's commit and push. Agents have not committed or pushed.

Read this file, HANDOFF.md, DIRECTIVE_COMPLIANCE.md, and `../current_deck.json`. The compact `handoff_context/` folder preserves review reports, slide titles and speaker notes, old-to-new numbering, source maps, solver provenance, and completion checks that otherwise live in ignored build folders. These are explicit records of the work, not a copy of the entire chat history. Some historical reports use source slide numbers; use slide_map.csv for final positions.

The current v21 PowerPoint SHA-256 is `49a25c6145961bc260f28f7b657887638a1e53674f456e5a4a0c15cbc1eaefad`. Its only change after the main delivery is the shared comparison image on slides 8 and 33. Both surface views use azimuth 142.5°, elevation 30°. Data, colors, and shared scales are unchanged. Keep v21; Ben has not requested v22. Preserve the rose-green C3NS template, original small logo/footer, honeycomb title, plot colors, native media, and hidden question slides. SADb stays at the end of questions. The talk has 77 main slides, a 38.5-minute baseline, and an absolute 45-minute limit; timing has not been rehearsed.

## Transfer before rebuilding

1. Ben commits and pushes the pending v21 scripts, HANDOFF.md, this file, handoff_context/, and ../current_deck.json. Pull that branch on EB475WS4. The snapshot handoff_context/ is about 314 kB before this pickup note.
2. Copy the files listed in `handoff_context/LOCAL_ASSET_TRANSFER.csv` from easteregg2 to the same paths relative to Defense on EB475WS4. This inventory covers 482 files, 421,086,495 bytes (about 402 MiB), including v20/v21 outputs, v19/v20/v21 assets, and v18 submitted-source/figure references. These files remain outside Git. No whole-repository transfer or simulation rerun is necessary to assemble the saved deck. The inventory is a path/size list, not a checksum manifest; verify the current PowerPoint against the hash above.
3. For merely editing the finished presentation in PowerPoint, the current PPTX is sufficient. For source-driven rebuilding, transfer the full inventory. Preserve older helper/source trees: package_current.py depends on the tracked v18 transplant_v18.py helper.
4. Copy the submitted dissertation PDF, advisor's annotated PDF, and Hunt's original presentation from the external paths listed in HANDOFF.md if absent on EB475WS4. The verifier expects the submitted PDF at that path and checks its recorded hash. New neural figures and selected evidence are already staged under v21/assets/neural; the entire external neural campaign is needed only for new analysis or simulation.

## Start on EB475WS4

Check hostname and follow its AGENTS.md machine bullet. Discover that machine's bundled artifact runtime and installed presentation skill; authoring scripts contain easteregg2 runtime paths that may need updating. Use EB475WS4's Python environment, not D:/Anaconda. PowerPoint and MATLAB should run locally there. Do not launch heavy searches or simulations as part of deck work.

The source-driven builder imports v20 and applies the v21 modules. Reuse the saved assets. The compact handoff snapshots do not replace the build directory: run the documented build sequence to recreate generated state. Do not run prepare_pipeline.py over the improved scripts. Do not blindly rerun finish_delivery.py, whose assertions describe the original final-render pass; account for the subsequent in-place correction and preserve the latest pointer. For future plot regeneration, the corrected build_legacy_bflh.m is copied into handoff_context/ as well as the transferred assets tree.

Render and inspect changed slides, preserve unrelated package content, and keep only the last two delivered revisions. Read the committee-perspective reviews as agent audits, not approval by the actual committee. The new guided human-model stepping, historical EB475WS4 standing schematics, and submitted dissertation results are distinct sources; HANDOFF.md records their scope. Do not modify the submitted dissertation or rerun scientific work merely to revise slides.

No deck-building or simulation job remains running for this handoff on easteregg2. The files have been inventoried here; their arrival on EB475WS4 has not been verified.
