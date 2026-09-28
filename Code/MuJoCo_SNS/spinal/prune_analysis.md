# Prune-matrix analysis (easteregg2, 2026-09-28)

## s3k  (79 cells, 0 errors)

| cut | AIR d | AIRDEAFF d | WALK d | STAND d | sway d (m) | push fell | push sway d (m) | verdict |
|---|---|---|---|---|---|---|---|---|
| noaff | -1.9 | +0.0 | -81.4 | -15.2 | -0.020 | 0/4 | -0.020 | **KEEP (walk, stand)** |
| interleg | -0.4 | +0.0 | +11.5 | +0.0 | +0.000 | 0/4 | -0.000 | **PRUNE** |
| contact | -0.4 | +0.0 | +2.5 | +0.0 | +0.000 | 0/4 | -0.000 | **PRUNE** |
| ia | +1.4 | -0.0 | -133.6 | -7.9 | -0.035 | 0/4 | -0.043 | **KEEP (walk)** |
| ii | -2.0 | +0.0 | -68.8 | -8.7 | -0.033 | 0/4 | -0.031 | **KEEP (walk)** |
| ib | -0.4 | -0.0 | +6.2 | +0.0 | +0.000 | 0/4 | -0.000 | **PRUNE** |
| renshaw | -0.4 | -0.0 | -80.8 | -2.8 | +0.007 | 0/4 | +0.007 | **KEEP (walk)** |
| rgweak | -0.4 | +0.0 | +0.0 | +0.0 | +0.000 | 0/4 | +0.000 | **PRUNE** |
| combo | +0.0 | n/a | +10.3 | +0.0 | +0.000 | 0/4 | -0.000 | **PRUNE** |

FULL reference: air -8.50842134773583 (rises 0, period Nones), walk kine -186.41420708638907 (duty 1.0, kz 0.86), stand score 31.5 (sway 0.128, fell False)

## w2lvar  (72 cells, 0 errors)

| cut | AIR d | AIRDEAFF d | WALK d | STAND d | sway d (m) | push fell | push sway d (m) | verdict |
|---|---|---|---|---|---|---|---|---|
| noaff | -48.3 | +0.0 | -80.0 | -20.8 | +0.034 | 0/4 | +0.047 | **KEEP (air, walk, stand, stand, push)** |
| interleg | -45.7 | +0.0 | -56.7 | -30.9 | +0.052 | 0/4 | +0.107 | **KEEP (air, walk, stand, stand, push)** |
| contact | -48.3 | +0.0 | -76.7 | -4.9 | +0.006 | 0/4 | +0.006 | **KEEP (air, walk)** |
| ia | -48.3 | -0.1 | -129.1 | -19.7 | +0.033 | 0/4 | +0.047 | **KEEP (air, walk, stand, stand, push)** |
| ii | -48.3 | +0.0 | -6.1 | -3.6 | +0.004 | 0/4 | +0.006 | **KEEP (air, walk)** |
| ib | -48.4 | -0.0 | -88.2 | +6.5 | -0.000 | 0/4 | +0.022 | **KEEP (air, walk, push)** |
| renshaw | -48.3 | +0.0 | +2.1 | +6.7 | -0.000 | 0/4 | +0.011 | **KEEP (air, push)** |
| rgweak | -48.3 | -0.0 | -135.5 | -4.9 | +0.006 | 0/4 | +0.006 | **KEEP (air, walk)** |

FULL reference: air 37.981074306639684 (rises 4, period 2.44s), walk kine -144.46440311078322 (duty 0.5, kz 0.79), stand score 25.6 (sway 0.098, fell False)

## syn6  (64 cells, 0 errors)

| cut | AIR d | AIRDEAFF d | WALK d | STAND d | sway d (m) | push fell | push sway d (m) | verdict |
|---|---|---|---|---|---|---|---|---|
| noaff | -46.8 | +0.0 | +9.8 | -18.1 | +0.048 | 0/4 | -0.027 | **KEEP (air, stand, stand)** |
| interleg | -46.8 | +0.0 | -12.2 | -0.1 | +0.000 | 0/4 | +0.001 | **KEEP (air, walk)** |
| contact | -46.7 | +0.0 | -11.4 | -19.1 | +0.050 | 0/4 | +0.010 | **KEEP (air, walk, stand, stand)** |
| ia | -46.7 | +0.0 | -4.2 | -0.2 | +0.000 | 0/4 | +0.002 | **KEEP (air, walk)** |
| ii | -46.7 | +0.0 | -4.2 | -0.2 | +0.000 | 0/4 | +0.002 | **KEEP (air, walk)** |
| ib | -46.7 | +0.0 | -4.2 | -0.2 | +0.000 | 0/4 | +0.002 | **KEEP (air, walk)** |
| rgweak | -46.9 | +0.0 | -119.6 | -0.4 | +0.000 | 0/4 | +0.001 | **KEEP (air, walk)** |

FULL reference: air 36.43985845560516 (rises 4, period 2.45s), walk kine -200.37955834081458 (duty 0.42, kz 0.84), stand score 37.5 (sway 0.088, fell False)

