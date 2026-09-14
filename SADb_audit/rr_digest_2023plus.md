# RR-style digest — recent papers (2023–2026) relevant to SADb — built 2026-09-12

Method: Ben's Research Rabbit collection was observed live (selected paper: Shevtsova
2026 eLife "Linking spinal circuit reorganization to recovery after thoracic SCI";
recommendations incl. Mari 2023, Audet 2023, Rybak 2024, Harrison 2023). Those plus the
lab's model family (Rybak 2013, Rybak 2024/25, Shevtsova 2025/26, V3/commissural lines)
were used as OpenAlex seeds; all works 2023+ citing any seed were pulled (137 candidates),
deduped against the 483-DOI SADb corpus, preprint/eLife-version twins collapsed, and the
shortlist below verified via Crossref. Full candidate dump: `batch3\rr_digest_candidates.json`.

NOT yet written to Airtable/Zotero (session lost its Airtable route + Zotero key).
Airtable import file: `batch3\rr_digest_shortlist.json` (Name/Author/Year/DOI fields ready).
Zotero: `python zotero_add_digest.py --key <KEY> --post` (payload `batch3\zotero_digest_items.json`).

## Shortlist (19)

| # | DOI | Year | First author | Title (short) | Why relevant |
|---|-----|------|-----------|---------------|--------------|
| 1 | 10.7554/elife.98841 | 2024 | Rybak | Operation regimes of spinal circuits controlling locomotion… | RR-recommended; the lab's model family, final eLife version |
| 2 | 10.7554/elife.103504 | 2025 | Rybak | Operation of spinal sensorimotor circuits controlling phase… | same family, sensory-phase transitions |
| 3 | 10.1152/jn.00104.2024 | 2024 | Harnie | Forelimb movements contribute to hindlimb cutaneous reflexes | cutaneous reflex network during locomotion |
| 4 | 10.1113/jp286151 | 2024 | Mari | Changes in intra/interlimb reflexes from hindlimb cutaneous afferents | RR-recommended (Mari 2023) follow-up; cutaneous gating |
| 5 | 10.1113/jp286808 | 2024 | Mari | Sister paper: forelimb cutaneous | same series |
| 6 | 10.1016/j.neunet.2024.106422 | 2024 | Zhu | Spinal circuit model, asymmetric cervical-lumbar layout | MODEL |
| 7 | 10.1098/rsos.240207 | 2024 | Molkov | Sensory feedback and central neuronal interactions | Rybak-lab synthesis of sensory+central control |
| 8 | 10.1152/jn.00248.2023 | 2023 | Refy | Dynamic spinal reflex adaptation during locomotor adaptation | reflex plasticity in task |
| 9 | 10.3389/fncir.2023.1235181 | 2023 | Chacon | Lumbar V3 interneurons provide direct excitatory input | extends the V3 (Zhang 2022) line already in SADb |
| 10 | 10.1016/j.cub.2023.07.014 | 2023 | Laflamme | Distinct roles of spinal commissural interneurons | experimental kin of Rybak 2013 left-right model |
| 11 | 10.1371/journal.pcbi.1012101 | 2025 | Pazzaglia | Balancing central control and sensory feedback | CPG+feedback model |
| 12 | 10.1371/journal.pcbi.1013494 | 2025 | Severini | Physiologically inspired hybrid CPG/reflex controller | robot control — directly the lab's application |
| 13 | 10.1016/j.cub.2025.09.030 | 2025 | Toscano | A spinal circuit for skilled locomotion | spinal circuit, recent |
| 14 | 10.1016/j.expneurol.2023.114496 | 2023 | Danner | Spinal control of locomotion before and after SCI | lab collaborator review |
| 15 | 10.1523/jneurosci.2015-22.2023 | 2023 | Giorgi | Excitatory and inhibitory descending commissural interneurons | commissural circuitry |
| 16 | 10.1007/s00422-023-00970-z | 2023 | Kohler | BCM rule lets a spinal cord model learn rhythmic | learning in CPG models |
| 17 | 10.1038/s42003-024-06843-w | 2024 | Arai | Interlimb coordination is not strictly controlled during walking | coordination control theory |
| 18 | 10.1152/jn.00331.2025 | 2025 | Yassine | Speed-dependent locomotor adjustments following staggered… | Frigon-group spinal-cat line |
| 19 | 10.1101/2025.11.11.687930 | 2025 | Shinohara | Adaptive interlimb coordination to sudden ground loss (cat CPG NMS model) | already in Ben's Zotero (C5LYIBNK) — Airtable row still needed |

## Observed RR items resolved
- Zhang 2022 V3 → already in SADb (batch 2, KIG9DEKJ). ✓ no action
- Shevtsova 2026 eLife 10.7554/elife.107480 → shortlist seed; its final version row
  appeared in candidates as 10.7554/elife.107480 (2025). Add to Airtable alongside #19
  (not double-counted — one row; it IS in the candidates file).
  NOTE: item 19 and this item are different papers; both were in Ben's tabs.
- Harrison 2023 podokinetic path integration → judged OUT of SADb scope (attention/path
  integration, not afferent circuitry). Skipped.
