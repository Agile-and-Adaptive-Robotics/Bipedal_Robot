# Submission audit — pending only, v3

30 October 2026. This supersedes the seven-page v2 action list. Rejected proposals and verified completed fixes are omitted. Source line numbers below refer to the current local ProofFinal files. Yellow highlighting is reserved for the items still open; approval of an action does not mean that action is complete.

## Dissertation: decisions or unfinished presentation work

1. **SNS diagrams: preview before insertion.** The requested enlargement is shown in `output/pdf/SNS_Diagram_Preview_20261030.pdf`; it has not replaced the active figures. Ben's request for neural-pathway schematics concerns the **Gait2392-Based Models** subsection, not the knee-demonstration subsection. Readable, configuration-specific Gait2392 circuits remain to be selected. Locations: §6.4.2, `60-sns.tex:122`; §6.4.3, `60-sns.tex:149` and `156`. (D03–D05; comments 3–6.)

2. **Earlier knee comparison: clarify its relation to the later experiment.** The first curve is a graded/non-spiking baseline; the second is its spiking twin. The retained five-second comparison is separate from the later settled, pulsed graded-neuron test. Its table still says “matched operating point” and “settle.” Proposed wording is “Earlier five-second startup comparison of the graded baseline and spiking twin; statistics describe the final two-second window.” No new experiment is requested. Location: §D.4.2, `95-AppendixD.tex:252`. (D03; comment 3.)

3. **Compact configuration map: presentation decision.** A small table could distinguish Gait2392 human biological-muscle tests, the BPA robot counterpart, the Li–Deng/W2L walker, and the CAD-derived robot, giving controller, actuator law, support, and demonstrated outcome. These are not interchangeable models or matched support conditions. Latest fetched Hunt-training results concern the Gait2392 human model, not AH001. The CAD record shows supported unloading/lifting gains but does not establish regular alternating footfalls or progression. Locations: §6.4.2, `60-sns.tex:122`; §D.3.1, `95-AppendixD.tex:133`; §D.5.4, `95-AppendixD.tex:359`. (D04; comment 5.)

4. **Route-geometry figure readability: approved, unfinished.** Retain direct-route and cam-route identities. If multiple views are used, show an initial bone-overlay view for orientation; add frame/axis units to the main pose explanation and improve remaining dense contact/transmission labels. Do not combine the two configurations as though their routes were identical. Locations: §4.3.3, `30-torque.tex:715`; §C.6, `94-AppendixC.tex:325` and `333`; §C.8, `94-AppendixC.tex:406`. (D09; comment 10.)

5. **Torque authorship/publication note: author confirmation.** The dissertation names seven authors, whereas the companion manuscript names four. The agent has not chosen authors or contributor roles. The published 13-page Morrow (2020) paper has no separate author-contribution statement; this does not establish each person's contribution. Retain the legacy placement work in the dissertation, as Ben directed. Locations: Chapter 4 opening, `30-torque.tex:5`; paper front matter, `Main.tex:48` and `69`. (D14/D16; comments 15 and 17.)

6. **Unit formatting: expressly deferred.** Handwritten units remain in Appendix C, with isolated instances in Appendix B and the paper supplement. No formatting sweep was made. Representative location: §C.4, `94-AppendixC.tex:240`. Ben asked that these remain unchanged for now.

## Companion paper: keep separate from dissertation submission

7. **Small plot labels: approved, unfinished.** The remaining six-panel pinned-extensor plot and supplementary route comparison need readable final-size labels. Preserve every trace and update a reproducible generator. Locations: §3.2, `Main.tex:441`; supplement §1.3, `frontiers_SupplementaryMaterial.tex:171`. (P02; comment 20.)

8. **Identification settings table: explanation before approval.** The proposal is to replace the long settings paragraph with two rows, one for flexor identification and one for extensor identification, recording test assignments, bounds, population/generation settings, and the shared stiffness constraint. It changes presentation, not optimization or FVU. No table has been inserted. Location: §2.8, `Main.tex:382`. (P04; comment 22.)

9. **Abstract hardware sentence: wording before approval.** The abstract currently singles out extensor hardware despite physical results for both routes. Proposed replacement: “At the measured configurations, both redesigned routes met or exceeded the torque magnitude of their selected human-muscle references.” It does not claim every pressure or unmeasured angle was tested. No abstract edit has been made. Locations: abstract, `Main.tex:81`; §4.3.2, `Main.tex:602`. (P05; comment 23.)

The two Chapter 5 presentation credits are implemented. Source citations for the three reproduced background images are retained, without a new license search. Historical scorer expansion, a new SADb release/snapshot statement, and moving hardware results into another section remain rejected. Zotero is untouched and its cleanup is paused by Ben. No broad paper replacement or live upload is included in this pass.
