# Human 5x5 grading surfaces index (2026-10-23)

One line per actuator: the two graded joint coordinates (with their
range-of-motion spans in degrees), the grid, the in-domain artifact-flag
count, and the 9x9 viewing export. The single grading source of truth is
`human_surfaces_5x5_grading_20261023.csv` (this folder; 25 cells per
actuator, direct OpenSim evaluations at the nodes; knee flexion NEGATIVE,
hip adduction POSITIVE). The 9x9 exports are bilinear interpolations from
the 1-deg OpenSim samples, for VIEWING only; percentages are computed on
the 5x5 grading grid, never on the 9x9. All 9x9 exports live in
`tmp/paper_bridge/capability_table_20261022/human_surfaces_9x9_20261023/`.

| id | Actuator | Graded coordinates (deg) | Grid | Flagged graded cells | 9x9 export |
|---|---|---|---|---|---|
| 1 | Erector Spinae | lumbar_extension [-20, 15] x lumbar_bending [-15, 15]; pinned at 0: lumbar_rotation | 5x5 grading (25 graded of 25); 9x9 view | 0 of 25 | `human_surfaces_9x9_20261023/act01_ercspn_9x9_20261023.csv` |
| 2 | Internal Oblique | lumbar_extension [-20, 15] x lumbar_bending [-15, 15] | 5x5 grading (25 graded of 25); 9x9 view | 0 of 25 | `human_surfaces_9x9_20261023/act02_intobl_9x9_20261023.csv` |
| 3 | External Oblique | lumbar_extension [-20, 15] x lumbar_bending [-15, 15] | 5x5 grading (25 graded of 25); 9x9 view | 0 of 25 | `human_surfaces_9x9_20261023/act03_extobl_9x9_20261023.csv` |
| 4 | Bicep Femoris, Long Head | knee_angle_r [-120, 10] x hip_flexion_r [-15, 110] | 5x5 grading (25 graded of 25); 9x9 view | 5 of 25 | `human_surfaces_9x9_20261023/act04_bifemlh_9x9_20261023.csv` |
| 5 | Semimembranosus | knee_angle_r [-120, 10] x hip_flexion_r [-15, 110] | 5x5 grading (25 graded of 25); 9x9 view | 13 of 25 | `human_surfaces_9x9_20261023/act05_semimem_9x9_20261023.csv` |
| 6 | Sartorius | knee_angle_r [-120, 10] x hip_flexion_r [-15, 110] | 5x5 grading (25 graded of 25); 9x9 view | 0 of 25 | `human_surfaces_9x9_20261023/act06_sar_9x9_20261023.csv` |
| 7 | Tensor Fasciae Latae | knee_angle_r [-120, 10] x hip_flexion_r [-15, 110]; pinned at 0: hip_adduction_r | 5x5 grading (25 graded of 25); 9x9 view | 11 of 25 | `human_surfaces_9x9_20261023/act07_tfl_9x9_20261023.csv` |
| 8 | Rectus Femoris | knee_angle_r [-120, 10] x hip_flexion_r [-15, 110] | 5x5 grading (25 graded of 25); 9x9 view | 13 of 25 | `human_surfaces_9x9_20261023/act08_rect_fem_9x9_20261023.csv` |
| 9 | Gracilis | knee_angle_r [-120, 10] x hip_flexion_r [-15, 110] | 5x5 grading (25 graded of 25); 9x9 view | 5 of 25 | `human_surfaces_9x9_20261023/act09_grac_9x9_20261023.csv` |
| 10 | Bicep Femoris, Short Head | knee_angle_r [-120, 10] x hip_flexion_r [-15, 110] | 5x5 grading (25 graded of 25); 9x9 view | 0 of 25 | `human_surfaces_9x9_20261023/act10_bifemsh_9x9_20261023.csv` |
| 11 | Vastus Intermedius | knee_angle_r [-120, 10] x hip_flexion_r [-15, 110] | 5x5 grading (25 graded of 25); 9x9 view | 10 of 25 | `human_surfaces_9x9_20261023/act11_vas_int_9x9_20261023.csv` |
| 12 | Soleus | ankle_angle_r [-50, 20] x knee_angle_r [-120, 10] | 5x5 grading (25 graded of 25); 9x9 view | 10 of 25 | `human_surfaces_9x9_20261023/act12_soleus_9x9_20261023.csv` |
| 13 | Medial Gastroocnemius | ankle_angle_r [-50, 20] x knee_angle_r [-120, 10] | 5x5 grading (25 graded of 25); 9x9 view | 16 of 25 | `human_surfaces_9x9_20261023/act13_med_gas_9x9_20261023.csv` |
| 14 | Lateral Gastrocenemius | ankle_angle_r [-50, 20] x knee_angle_r [-120, 10] | 5x5 grading (25 graded of 25); 9x9 view | 16 of 25 | `human_surfaces_9x9_20261023/act14_lat_gas_9x9_20261023.csv` |
| 15 | Tibialis Posterior | subtalar_angle_r [-25, 35] x ankle_angle_r [-50, 20] | 5x5 grading (25 graded of 25); 9x9 view | 1 of 25 | `human_surfaces_9x9_20261023/act15_tib_post_9x9_20261023.csv` |
| 16 | Tibialis Posterior (duplicate name) | ankle_angle_r [-50, 20] x subtalar_angle_r [-25, 35] | 5x5 grading (25 graded of 25); 9x9 view | 0 of 25 | `human_surfaces_9x9_20261023/act16_tib_ant_9x9_20261023.csv` |
| 17 | Peroneus Brevis | subtalar_angle_r [-25, 35] x ankle_angle_r [-50, 20] | 5x5 grading (25 graded of 25); 9x9 view | 3 of 25 | `human_surfaces_9x9_20261023/act17_per_brev_9x9_20261023.csv` |
| 18 | Peroneus Longus | subtalar_angle_r [-25, 35] x ankle_angle_r [-50, 20] | 5x5 grading (25 graded of 25); 9x9 view | 3 of 25 | `human_surfaces_9x9_20261023/act18_per_long_9x9_20261023.csv` |
| 19 | Peroneus Tertius | subtalar_angle_r [-25, 35] x ankle_angle_r [-50, 20] | 5x5 grading (25 graded of 25); 9x9 view | 0 of 25 | `human_surfaces_9x9_20261023/act19_per_tert_9x9_20261023.csv` |
| 20 | Flexor Digitorum Longus | mtp_angle_r [-40, 40] x ankle_angle_r [-50, 20]; pinned at 0: subtalar_angle_r | 5x5 grading (25 graded of 25); 9x9 view | 0 of 25 | `human_surfaces_9x9_20261023/act20_flex_dig_9x9_20261023.csv` |
| 21 | Flexor Hallucis Longus | mtp_angle_r [-40, 40] x ankle_angle_r [-50, 20]; pinned at 0: subtalar_angle_r | 5x5 grading (25 graded of 25); 9x9 view | 0 of 25 | `human_surfaces_9x9_20261023/act21_flex_hal_9x9_20261023.csv` |
| 22 | Extensor Digitorum Longus | mtp_angle_r [-40, 40] x ankle_angle_r [-50, 20]; pinned at 0: subtalar_angle_r | 5x5 grading (25 graded of 25); 9x9 view | 0 of 25 | `human_surfaces_9x9_20261023/act22_ext_dig_9x9_20261023.csv` |
| 23 | Extensor Hallucis Longus | mtp_angle_r [-40, 40] x ankle_angle_r [-50, 20]; pinned at 0: subtalar_angle_r | 5x5 grading (25 graded of 25); 9x9 view | 0 of 25 | `human_surfaces_9x9_20261023/act23_ext_hal_9x9_20261023.csv` |
| 24 | Gluteus Maximus | hip_flexion_r [-15, 110] x hip_adduction_r [-45, 20]; envelope-coupled adduction axis; pinned at 0: hip_rotation_r | 5x5 grading (15 graded of 25); 9x9 view | 3 of 15 | `human_surfaces_9x9_20261023/act24_glut_max1_9x9_20261023.csv` |
| 25 | Adductor Magnus | hip_adduction_r [-45, 20] x hip_flexion_r [-15, 110]; envelope-coupled adduction axis | 5x5 grading (15 graded of 25); 9x9 view | 3 of 15 | `human_surfaces_9x9_20261023/act25_add_mag1_9x9_20261023.csv` |
| 26 | Iliacus | hip_flexion_r [-15, 110] x hip_adduction_r [-45, 20]; envelope-coupled adduction axis | 5x5 grading (15 graded of 25); 9x9 view | 7 of 15 | `human_surfaces_9x9_20261023/act26_iliacus_9x9_20261023.csv` |
| 27 | Psoas X | hip_flexion_r [-15, 110] x hip_adduction_r [-45, 20]; envelope-coupled adduction axis | 5x5 grading (15 graded of 25); 9x9 view | 7 of 15 | `human_surfaces_9x9_20261023/act27_psoas_9x9_20261023.csv` |

Totals: 635 graded cells, 40 envelope-excluded cells
(acts 24-27, 10 each), 126 flagged graded cells. Flag layer,
hinge-chain, gate, and provenance detail per actuator:
`human_surfaces_manifest_20261023.json` (this folder). Pinned DOFs: actuator
1 lumbar_rotation; 7 hip_adduction_r; 20-23 subtalar_angle_r; 24
hip_rotation_r.
