# Chapter 6 and Appendix D revision - 30 October 2026

## Implemented

- Sections 6.1-6.3 now present the AnimatLab result, controller-transfer problem, and need for scripted/headless and graphical tools before introducing MATLAB as a solution. API, SDK, GUI, and headless execution are briefly defined. No blanket claim that Python lacks graphical programming is made.
- The Methods and Verification section follows the eight rows of the updated Table 6.1. The existing seven tool rows are retained; execution speed and actuator implementations are added last.
- A short native AnimatLab phase-1 recording moves from Appendix D into the chapter opening. It remains explicitly distinct from the original shared-RG course-project air-stepping example. Actual recorded RG/PF/MN channels and joint angles are retained without interpolation or smoothing.
- The old cross-engine trace figure now overlays antagonists and implementations on two axes over 2-5 s. The sensory-feedback figure now aligns RG, hip, knee, ankle, and heel traces over 12-16 s. PF/MN channels were not recorded in that run and are not manufactured.
- The reduced Deng-derived PF-stimulation plot aligns RG/PF/MN signals, joint kinematics, and applied current. It is a translated suspended-model diagnostic, not a reproduction of measured rat gait.
- D1 now says neural subsystems, with bold (A)/(B)/(C) before their descriptions. Old D6 and D7 become one coupled graded/spiking comparison over 10-14 s; their labels both point to the same retained figure. Filtered spiking readouts are identified correctly. The old D8 startup figure/table are archived in full and no longer presented as current reflex evidence.
- Hunt's four neural-training operations are combined with the five physical evaluation stages. The final stage is refinement/support reduction, and unpassed stages remain identified as unpassed.
- Runtime statements identify their workloads: neural-only 11.23x versus the preceding preallocated kernel; 16 s coupled pilots 69.51 to 15.86 s and 113.74 to 34.59 s; latest isolated compliant-muscle bank 16 s in 3.55 s. The approximate 240x and 8x claims could not be verified and are not inserted.
- The chapter ends with Ben's confirmed SNS-Simscape CAD/BPA virtual-harness-assisted standing example. It does not borrow the separate ideal-tension MuJoCo CAD model's 80-percent assistance or other measurements.

## Preservation and explicit archival

No simulation, optimization, model fitting, or experimental rerun occurred. Figures use saved data; file hashes, channels, windows, and font profiles are recorded with their generators in Figures/60-sns/Components/sns_redesign_20261030. Source arrays and scientific metrics were not changed. Old graphic files are in the respective Deprecated folders. Complete earlier startup source blocks and pre-edit sources are retained locally in D:/Bipedal_Robot_Local_Artifacts/sns_rewrite_20261030.

Published Chapter 20, the torque-paper working sources, both downloaded Overleaf ZIPs, and Zotero are outside this pass. No live upload, commit, or push.

## Remaining work

The pending-only v5 audit lists remaining work with current source lines. In particular, the available phase-1 model contains an initial horizontal-force stimulus, not a documented timed platform drop. The original course-project antagonist/frequency recording and a verified drop response remain pending. The later instrumented recording is shown without claiming either.

## Final verification

Clean build: 298 pages, zero errors or undefined references. Five placed figure minima: 10.120, 10.205, 10.246, 10.110, and 10.120 pt; no visual clipping or label collisions. Both original downloaded ZIPs and protected Chapter 20 retain their checked hashes. Latest PDF and audit v5 links are in FINAL_DELIVERY_LINKS.md.
