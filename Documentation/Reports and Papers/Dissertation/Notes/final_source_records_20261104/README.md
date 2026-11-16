# Final dissertation revision workspace

Session 1, easteregg2. Created 4 November 2026 from ProofFinal. Main entry: main.tex; chapter sources: chapters/. Shared figure trees remain ../Figures. Root support/bibliography files were copied; only the Final bibliography received the verified Rybak 2024 entry.

Initial revision changes chapters 60-sns.tex, 70-overall.tex, and 95-AppendixD.tex. Chapter 6 now follows motivation -> knowledge gap -> human-model workflow -> graphical/CAD workflow -> controller design/training/execution -> W2L functional demonstrations -> discussion. Removed bold Verification lead-ins and retained every existing figure inclusion and label. Existing result values are retained; the new incremental-speed table uses verified external benchmarks, with whole-run timing stated separately.

Chapter 7 and Appendix D now explicitly describe walker_io. The historical human zero-drive claim is corrected: the additional left/right hand-drive scale reached zero while shared DRIVE remained 2.0773 nA. Evidence: Code/MuJoCo_SNS/spinal/wean_tune_20261016.py, header lines 11-14 and run_eval --drive BASE_DRIVE. This changes interpretation, not the archived numeric results.

The newly prepared four-command/central-afferent controllers are not described as trained or physically validated results. Their construction records remain D:/Users/Ben Bolen/ControllerPreparation_20261104. Postural premotor replacement remains pending. No MATLAB jobs, archived experiments, ProofFinal source, or live Overleaf files were changed.

SOURCE_CONSERVATION.json records preserved labels/figure inclusions. Published chapters and figures are unchanged. This is a revision draft, not an uploaded replacement. See BUILD_REVIEW.md for compilation and per-directive checks after verification.
