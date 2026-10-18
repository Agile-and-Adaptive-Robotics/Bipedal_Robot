# 80/20 15-Series Connection Design Loads — Reference Sheet

Compiled 2026-10-16 (ZCode session) after 80/20 removed all quantitative connection-load
data from their current materials. Original question: design load for 1515-series bracket
connections using 2 vs. 4 vs. 8 of the 5/16-18 fasteners.

## Why this data is hard to find now

Verified 2026-10-15 by searching, downloading, and text-scanning: 80/20 Product Catalog 26
(876 pp), Catalog 23 fasteners section, 2006-era catalog sections, the live 8020.net
bracket product/category pages, and Wayback-archived product pages back to 2019. None
publish load values for brackets or internal fasteners — current language is qualitative
only ("strong connection", "tighten to the required torque"). The old online
**Fastener Application Test** tool (referenced by Catalog 23, fastener-methods intro, and
the 2023 80/20 blog "T-Slot Framing Fasteners: Internal vs. External") is gone from the
live site and not recoverable from the Wayback Machine.

The numbers below survived in **ExtrusionBuilder News (EBN) Spring 2015**, published by
F&L Industrial Solutions — an authorized 80/20 distributor — presenting 80/20 fastener
test results. Archived copy: [`8020_fastener_strength_EBN_Spring2015.pdf`](8020_fastener_strength_EBN_Spring2015.pdf)
(next to this file). Mirror it was recovered from:
<https://www.greenstello.com/9dd65ba1477d4ba99dbc23002bfabca7/ebn-spring15.pdf>

## Tested connection strength — 1515-Lite profile, 6 in. sample (lbf)

| Connection method | Direct force | Cantilevered | Torsional | Avg |
|---|---:|---:|---:|---:|
| Single anchor fastener | 950 | 625 | 540 | 705 |
| Double anchor | 1,200 | 700 | 1,150 | 1,017 |
| End fastener (5/16-18) | 1,000 | 500 | 680 | 727 |
| Joining plate | 225 | 200 | 1,000 | 475 |
| 90° joining plate | 250 | 200 | 1,120 | 523 |
| **Corner bracket (2-hole, 2× 5/16-18)** | **375** | **225** | **500** | 367 |
| Corner gusset (2-hole) | 375 | 750 | 500 | 542 |
| Hidden corner (inside-inside) | 50 | 50 | 50 | 50 |
| Hidden corner (inside-outside) | 240 | 220 | 240 | 233 |

Test-article note: conducted on 1515-**Lite** (thinner slot walls than standard 1515) —
standard 1515 is at least this good. Per the source, results are similar for most bar
types. These are ultimate/failure-style test values, not allowable stresses — apply your
own factor of safety.

## Scaling by fastener count (2 vs. 4 vs. 8 of the 5/16-18)

80/20 never published 4-hole or 8-hole bracket values. Governing mode for direct and
cantilevered loads is **slip of the bracket against the slot face**, so capacity scales
~linearly with bolts per interface:

| Config | Bolts per leg | Direct-force estimate | Basis |
|---|---:|---:|---|
| 2-hole (4302-style) | 1 | 375 lbf (tested) | EBN Spring 2015 |
| 4-hole (2 per leg) | 2 | ≈ 700–750 lbf | 2× slip capacity |
| 8-hole (4304-style) | 4 | ≈ 1,400–1,500 lbf | 4× slip capacity |

First-principles check: 5/16-18 button-head screw at 80/20's specified torque
(10–15 ft-lbf, part 3320 assembly) → ≈ 2,300 lbf preload; at μ ≈ 0.10–0.15 on anodized
aluminum → 230–350 lbf slip resistance per screw per interface, consistent with the
tested 375 lbf for one-bolt-per-leg.

**Caveats:**
- Torsional stiffness does NOT scale linearly with bolt count (the gusset's 750 lbf
  cantilever value comes from its stiffness, not extra bolts).
- At high bolt counts the failure mode can shift from slip to T-slot wall bearing/tear-out.
- Torque must actually be reached — under-torqued connections get a fraction of these values.
- For a quotable primary source, 80/20 tech support will email current connection test
  data on request (they retired the public Fastener Application Test).

## Rig application note (knee test setup)

BPA routes are a few hundred lbf each → a 4-bolt bracket (~700 lbf direct) has reasonable
margin; 8-bolt is comfortable. Anchor fasteners (strongest per dollar/hole) reach
950–1,200 lbf where the joint geometry allows the counterbore machining.
