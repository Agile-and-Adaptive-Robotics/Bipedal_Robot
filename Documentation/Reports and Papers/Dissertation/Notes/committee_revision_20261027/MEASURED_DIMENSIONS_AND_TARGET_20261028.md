# Final-test measured dimensions and human-target wording

The updated Results_table_20mm.xlsx supplies separate physical dimensions for the final test pairs. This supersedes the earlier report that those measurement cells were empty. A snapshot and cell-level extraction are preserved in Archive/final_test_dimensions_20261028.

| Actuator | Resting length (mm) | Tendon length (mm) |
|---|---:|---:|
| Flexor 1 | 561 | 28 |
| Flexor 2 | 560 | 30 |
| Extensor 1 | 828 | 92 |
| Extensor 2 | 826 | 85 |

Sources: FlxTest20mm_56cm and ExtTest20mm_2, C2/D2 for resting length and C4/D4 for tendon length. Resting values are stored in meters. Ben confirmed that extensor entries 0.92/0.85 mean 92/85 mm; that same worksheet conversion is applied to flexor entries 0.28/0.30 as 28/30 mm. The conversion and source workbook SHA256 are recorded in measured_dimensions.json. The workbook itself was read only.

Table C7 now contains these physical measurements. Its explanatory paragraph keeps the original prediction inputs separate: shared flexor resting length 561.49 mm and effective tendon parameter 5.00 mm after the 20 mm fitting adjustment, and extensor resting length 827.74 mm with approximate tendon parameter 92.49 mm. The main chapter's cross-reference now describes measured dimensions and shared prediction inputs accurately. The paper includes the measured dimensions in its physical-test description. Both physical-test captions identify 562 mm as the nominal flexor plot label, avoiding its interpretation as the separate installed lengths.

Ben also supplied a close-up of the extensor measurement near -86.5 degrees. The blue measured point is above the black human target. Both manuscripts retain the statement that all nine extensor measurements meet or exceed the human reference. The pointwise sentence now explicitly says that the -86.5-degree measurement exceeds the human reference and is 2.06 N m below the pressure-specific prediction. That prediction-error comparison must not be interpreted as failure to meet the human target. The supplied plot is preserved as user_extensor_target_detail.png.

No measured torque, pressure, target, model coefficient, predicted curve, or recorded goodness-of-fit value was changed. No model was reevaluated. Graphics conservation verifies all 90 dissertation graphic files and 19 paper/supplement graphic files against the prior deliveries. This pass changes metadata, captions, a cross-reference, and explanatory wording only. The earlier 11-versus-21 mm biomimetic flexor tendon and 520-versus-518 mm extensor resting-length metadata differences remain separately documented; they were not silently reconciled.

Verification: fresh dissertation build passes five steps with no errors or undefined references, 321 pages. Paper and supplement compile successfully with no undefined references, 33 and 10 pages. Inherited local typesetting/font messages remain. Final live-source readbacks match local 30-torque.tex, 94-AppendixC.tex, and Main.tex. Live compile status and screenshots are stored with the pass evidence. Root delivery PDFs and ZIPs were updated and checked. Build products are archived outside ProofFinal. No tracked changes accepted/rejected, no commits, and no pushes.

Final live verification: dissertation 322 pages, Errors 0 / Warnings 0; paper 33 pages, Errors 0 / Warnings 19 (existing font warnings). The live paper was left at the physical-test paragraph that explicitly states all nine measurements meet/exceed the human reference and that the -86.5-degree point remains above that target.
