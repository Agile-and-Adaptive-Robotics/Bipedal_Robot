# Chapter6 and AppendixB/D: second PDF comments

These dispositions apply to the19 comments in USER_COMMENTS_19.json, not the earlier47-comment numbering. Ben authorized local source changes in this pass. Only60-sns.tex,93-AppendixB.tex,95-AppendixD.tex, and the separately authorized pipeline.pdf are changed. The existing corrected motor PowerPoint is verified without modification. No simulation, model builder, fit, optimization, data, or figure generator was run.

| Directive | Change and evidence |
|---|---|
| Framing context | Chapter6 starts with controller transfer and hand tuning, air stepping rather than walking, and the need to test alternatives on Gait2392/BPA models. AnimatLab is identified as usable but dated; MATLAB is explicitly a later route. 60-sns.tex:4; 60-sns.tex:14; 60-sns.tex:29. |
| Table6.1 order | Methods and Results each have seven subsections in table-row order: W2L; Gait2392; blocks; CAD assembly; bridge; curriculum; literature/spiking. Existing results and figure blocks were moved intact. 60-sns.tex:72,145; 60-sns.tex:137,229. |
| C8: repeated capability framing | Removed repeated prose inventories and the repeated gap/workflow recap; the capability table remains exact. All recorded result statements and limitations remain. 60-sns.tex:37; 60-sns.tex:241. |
| C10: clarify/condense AppendixB | Shortened the toolbox, mirror implementation, burst-termination, and conversion instructions; retained every equation, numeric token, parameter/calibration/gate/variant table, graphic, and label. 93-AppendixB.tex:70; 93-AppendixB.tex:86; 93-AppendixB.tex:166. |
| C10: correct FigureB1 letters | Native artwork has three row labels(A1),(B1),(C1), with hip left and knee right. Caption/body now use those actual labels instead of six nonexistent A–F panel labels. The asset is unchanged. 93-AppendixB.tex:196; 93-AppendixB.tex:190. |
| C12: orphaned material | This meant existing content without useful explicit references, not missing figures or sections. All six identified AppendixD figures now have explanatory local calls and direct main-text calls (inventory below). No content was recreated or removed. |
| C14: c1/V3 correction | The primary Shinohara network explicitly contains C1 and V3 commissural interneuron populations. The prior implementation-only description was unsupported and is withdrawn. A brief function gloss now cites the existing primary-source key. 93-AppendixB.tex:207. |
| C16: conversion evidence | Main text directly calls the normalized force traces and absolute-statistics table. The comparison concerns the repaired MyoConverter-derived human model at matched poses/activations, with muscle-specific exceptions retained. It does not test dynamic tendon equivalence or isolate converter errors from later model repairs. 60-sns.tex:173; 95-AppendixD.tex:409. |

## Six figure references

| Existing supporting figure | Main text | AppendixD |
|---|---|---|
| fig:animatlab_phase1_preliminary | 60-sns.tex:12 | 95-AppendixD.tex:27 |
| fig:mujoco_model_views | 60-sns.tex:90 | 95-AppendixD.tex:36 |
| fig:app_knee_transient | 60-sns.tex:188 | 95-AppendixD.tex:66 |
| fig:layer_stimulus_validation | 60-sns.tex:235 | 95-AppendixD.tex:85 |
| fig:app_w2l_perturbations | 60-sns.tex:160 | 95-AppendixD.tex:316 |
| fig:force-compare | 60-sns.tex:173 | 95-AppendixD.tex:409 |

## Primary-source and implementation checks

[Shinohara et al.(2025), neural-controller description](https://pmc.ncbi.nlm.nih.gov/articles/PMC12642437/) identifies C1-mediated inhibition between flexor rhythm generators and a V3 pathway through contralateral extensor interneurons. The saved primary full-text extraction independently confirms this at Code/MuJoCo_SNS/spinal/shinohara_2025_biophiv687930_fulltext.txt:294–298. The local W2L builder specifies the corresponding inhibitory c1 and excitatory V3 connections at w2l_mujoco/build_w2l_split_net.py:119–126. The text describes literature-supported modeled populations; it does not invent a new genetic identity for C1.

The executable body of Code/MuJoCo_SNS/spinal/_muscle_force_compare.py:78–105 loads OpenSim static-optimization forces/activations, assigns each stored pose and activation, calls mj_forward, and records actuator forces. The plotting block uses the matched-activation arrays and normalizes each curve to its own peak. No trajectory velocities or time integration of activation/fiber dynamics are supplied. Its opening docstring and some historical legend language overstate the compliant-tendon comparison; the manuscript retains the accurate pose-wise scope. The script was inspected, not run.

## Editable motor and pipeline footer

C7 is already satisfied by Figures/15-background/Components/actuator_technology_comparison.pptx. Its one slide contains native editable shapes/connectors and zero raster pictures. Six8-point stator squares have centers33points from the common motor-ring center, and match the active PDF within0.003points. Both file hashes match the earlier verified assets. See asset_evidence/MOTOR_EDITABLE_VERIFIED.json and motor_current.png. No rebuild.

The active pipeline.pdf contained the real footer “Computer-in-the-loop: test embedded control against the simulated body / before transfer to the physical robot.” Only those two text lines are removed, followed by cropping their empty strip (468×303 original MediaBox;468×260 CropBox). All vector drawing records and remaining text/font positions are identical, and the visible diagram has zero changed raster channels at2× resolution. The substantive idea remains for the root to retain in Chapter7 body/caption. No design_workflow.pdf change: its alleged second caption was an unsupported diagnosis. See asset_evidence/PIPELINE_FOOTER_CONSERVATION.json and native before/after images.

## Conservation and verification boundary

44/44 source checks pass. All recorded numeric tokens, equations, table blocks, graphic inclusions and sizing, old labels, and negative evidence are retained. Chapter6 figure blocks are byte-identical after reordering; AppendixB figure blocks differ only in the approved FigureB1 caption. AppendixD figures/tables and all image assets are unchanged. Before copies, exact diffs, hashes, and checks are in before_sns/,*.tex.diff,SNS_CHANGED_FILES.json,andSNS_CONSERVATION.json.

- 60-sns.tex: 4419→3968 whitespace-delimited source words. Final SHA-256 3e86db84fcf3f8379a332d268288a036b0dc18735f84bed3ca14ca7ee7a2c6e8.
- 93-AppendixB.tex: 2680→2435 whitespace-delimited source words. Final SHA-256 86ee57048a9268ae06d47d11f4b83cbbcebbc61a842878b9fdc99409685776f4.
- 95-AppendixD.tex: 5513→5579 whitespace-delimited source words. Final SHA-256 b9ffbe47f68153413b92c4036e11ca3765c6145c7c1696f8cd166d98c9c73a6f.

Fresh integrated dissertation build and page-level layout/cross-reference inspection remain the root’s final verification step; no clean-build claim is made in this source-only record.

Recommended Chapter7 Aim4 future-work scope (root-owned): model-specific tests of alternative literature controllers; dynamic force comparisons with matched velocity/activation histories and tendon assumptions; tonic-drive/full-export discrepancies; spiking-network retuning; and reduction of the position-control harness as well as constant weight assistance. Retain DiRusso and spiking failures and the separate force-prediction exceptions. Do not present those future studies as completed.
