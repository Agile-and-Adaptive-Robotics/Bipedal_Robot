# Literature rules conflict audit (2026-09-30)

Cross-comparison of rybak_rules.json, shevtsova_rules.json,
shinohara_rules.json vs ben_rules_20260924.json (the four
connectome rule sets in spinal\). Normalized homolog labels;
143 distinct (src,dst) families across the sets.

## SIGN CONFLICTS (8)

| src | dst | per-file |
|---|---|---|
| ia | ia | rybak:inh(g=0.1) [rybak2006_iainmut]; rybak:inh(g=0.1) [rybak2006_iainmut]; ben:exc(g=1.0) [ia_recip]; ben:inh(g=0.5) [] |
| ia | mn | rybak:inh(g=0.6) [rybak2006_recip]; rybak:inh(g=0.6) [rybak2006_recip]; ben:exc(g=2.0) [ia_homo (0 hops shown)]; ben:inh(g=0.5) [ia_recip hops=1] |
| ib | ib | ben:exc(g=1.0) [ib_auto]; ben:exc(g=1.0) [ib_rev]; ben:inh(g=0.5) []; ben:exc(g=0.5) [] |
| ib | mn | ben:inh(g=0.5) [ib_auto hops=1]; ben:exc(g=0.5) [ib_rev hops=1]; ben:inh(g=0.5) [] |
| ii | mn | ben:exc(g=0.5) [ii_exc hops=1]; ben:inh(g=0.5) [ii_inh hops=1] |
| rc | rc | ben:exc(g=1.0) [rc]; ben:exc(g=0.5) []; ben:inh(g=0.5) []; ben:inh(g=0.5) []; ben:inh(g=0.5) []; ben:inh(g=0.5) [] |
| rg-e | pf e | rybak:exc(g=0.0075) [rybak2006_rgpf]; rybak:inh(g=0.05) [rybak2006_pfgate]; shinohara:exc(g=0.7) [shin2025_a1] |
| rg-e | rg-e | rybak:exc(g=0.0125) [rybak2006_rgrec]; rybak:inh(g=0.115) [rybak2006_lam] |

## GAIN DISAGREEMENTS >2x, same sign (6)

| src | dst | per-file |
|---|---|---|
| pf e | mn | rybak:exc(g=0.5); shinohara:exc(g=0.07); shinohara:exc(g=0.12); shinohara:exc(g=0.08); shinohara:exc(g=0.12) |
| pf f | mn | rybak:exc(g=0.5); shinohara:exc(g=0.09); shinohara:exc(g=0.1); shinohara:exc(g=0.08) |
| rg f | pf f | rybak:exc(g=0.0075); shinohara:exc(g=0.7) |
| rg-e | v3 | shevtsova:exc(g=0.35); shevtsova:exc(g=0.35); shinohara:exc(g=0.5); ben:exc(g=1.0) |
| supra | rg f | rybak:exc(g=1.0); shinohara:exc(g=0.02) |
| supra | rg-e | rybak:exc(g=1.0); shinohara:exc(g=0.15) |

## Agreements (same sign, compatible gains; 0)

