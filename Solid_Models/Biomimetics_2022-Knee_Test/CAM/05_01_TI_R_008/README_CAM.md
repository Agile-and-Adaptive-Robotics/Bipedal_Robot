# CAM package — 05_01_TI_R_008 (right tibia plate)

Built 2026-10-16 by ZCode. Everything here came from your actual part + drawing:
the drawing was exported to PDF/DXF straight out of SolidWorks, the part to
STEP/STL, and every coordinate in the G-code was computed from the extracted
geometry (not read off the drawing by eye), then cross-checked against the
drawing dimensions. The generation and checking scripts are included so this can
be regenerated in minutes when the part revs to 009.

## Files

| File | What it is |
|---|---|
| `05_01_TI_R_008_OP10_TOP_3AXIS.nc` | 3-axis program, first setup (part top up) |
| `05_01_TI_R_008_OP20_BOTTOM_3AXIS.nc` | 3-axis program, second setup (flipped) |
| `05_01_TI_R_008_2AXIS_XY.nc` | 2-axis (XY CNC, manual Z) program for the American mill |
| `05_01_TI_R_008_drawing.pdf/.png/.dxf` | Your drawing, exported 2026-10-16 |
| `05_01_TI_R_008.step/.STL` | Part exports — usable in any CAM seat (Fusion, Mastercam, SW CAM) |
| `part_geometry.json` | The extracted geometry + drawing cross-check (readable) |
| `gen_cam.py` | Generator: part geometry in, all three .nc files out |
| `verify_nc.py`, `wall_guard.py` | The two automated checkers (see Verification) |

## The part, as extracted

Rectangular stepped plate, 37.83 × 67.90 × 36.43 mm, made of four prisms plus
four holes (all mm, in the STL frame where x 0..37.83 across the width, y 0..67.90
along the length, z 0 at the rail bottom):

- **Top plate** — full 37.83 × 67.90 footprint, z 24.000 → 36.425 (12.43 thick)
- **Pad A** — y 0 → 30.00, z 13.475 → 24.000
- **Pad B** — y 48.65 → 67.90, z 14.615 → 24.000
- **Base rail** — y 0 → 5.00, z 0 → 13.475
- **Holes** (x from the left edge): Ø7.00 THRU-plate+padA @ (18.915, 17.500),
  22.95 deep; Ø6.60 THRU-plate+padB @ (18.915, 57.400), 21.82 deep; Ø6.60 thru
  plate only @ (28.214, 39.325) and (10.414, 39.325), 12.43 deep.

Every extracted number ties to a printed drawing dimension: 22.95 = plate top to
pad A bottom, 21.81 = to pad B bottom, 18.65 = middle gap width, 19.25 = pad B
depth, 10.50 = pad A step, 17.80 = the two middle holes apart, 28.58 = top edge
to the middle hole line. All four hole diameters extract to within 0.002 mm of
the callouts.

## Assumptions (change these when you know the machine)

You said you had no machine details, so I picked conservative defaults and wrote
them all down:

1. **Material 6061-T6 aluminum** — the drawing title block has no material. If
   it's actually steel, cut the feeds in the files roughly in half and the RPM by
   ~60% (formulas below).
2. **Imperial, Fanuc/HAAS-style RS-274** (G20, G54, G81/G83, T# M06, M08 flood).
   Runs on Haas, Fanuc, LinuxCNC/PathPilot, Centroid, and ProtoTRAK in EIA mode.
3. **Stock 1.500 × 2.750 × 1.500 in** bar — 0.005 in/side over the part in X,
   0.039 in/side in Y, 0.066 in on the face. Saw to ~2.80 length.
4. **Hole drills**: the holes are metric (6.60 / 7.00 mm). For an imperial tool
   crib: **letter G (0.261 in)** for the three 6.60 holes and **letter J
   (0.277 in)** for the 7.00 hole — that's 6.63/7.04 mm, normal drill-oversize
   clearance. If you have metric drills, use 6.60/7.00 directly.
5. **1/2 in 3-flute carbide end mill** for everything else; a 1/2 in 90° spot
   drill; a 1/2 in 90° chamfer tool (optional, hole mouths only).

Coordinate frame in ALL programs: **G54 X0 Y0 = center of the part footprint;
Z0 = the finished face that is up for that setup.** Long direction along Y.

## OP10 — top side (3-axis)

Fixture: stock attached to a **sacrificial plate** (hot-melt glue or
double-sided tape; bolting through the middle-gap region also works since that
area is removed later). Do NOT put parallels under the footprint — the profile
cuts 0.020 in into the sacrifice.

Order: face top (two 0.033 passes, stock is 0.066 over) → rough outer profile
(0.150 steps, 0.010 wall stock) → finish profile full depth (tool center exactly
one radius off the wall, cuts to Z −1.454) → spot 4 holes → peck-drill (G83,
0.150 pecks) → chamfer hole mouths. The finish contour plunges in air 0.15
outside the corner, then leads in.

## OP20 — bottom side (3-axis)

Flip part bottom-up on the sacrificial plate. Z0 = the surface now on top (the
former rail-bottom plane). Re-find X-Y off the two machined long side walls
(indicate or use an edge finder on the finished walls; the part is symmetric
about X0).

Three open pockets, zig-zag roughing at 0.200 axial / 0.200 stepover with 0.010
left on walls and floor, then a finish pass per pocket:
**P1** beside the rail (y 5→30 mm, 13.475 deep — removes the overhang material
under pad A), **P2** the middle gap (y 30→48.65, 24.000 deep), **P3** beyond
pad B (y 48.65→67.9, 14.615 deep). The floors of P1/P2/P3 ARE the pad/plate
bottom faces; the pocket walls at y = 5, 30, 48.65 are the finished step faces.

## 2-axis program (your American mill)

Same geometry, structured for XY-CNC + manual quill: the file positions X-Y and
feeds, and every depth is a **[SET Z] checkpoint** — run the contour, raise
clear, drop the quill 0.100, re-run. Stages: 1) hole positions (spot + peck to
the depth stop; depths listed per hole in the header), 2) outer profile
(15 passes to −1.454), 3) flip, then the three pockets as lane passes with their
own schedules (P1 6 passes, P2 10, P3 6). If your 2-axis control *does* have a
powered Z, the 3-axis OP10/OP20 files are the better starting point.

## Speeds & feeds (and how to rescale)

| Tool | S | F | basis |
|---|---|---|---|
| 1/2 3FL carbide EM | 5200 | 45 in/min (12 plunge) | ~680 SFM, 0.0029 in/tooth |
| 1/2 spot drill | 3400 | 5 | ~445 SFM |
| Letter G / J drill | 4200 / 4000 | 8 | ~285 SFM, 0.002 in/rev |

RPM = SFM × 3.82 / diameter(in). If the mill tops out lower (a manual
Bridgeport does ~2500–3000), scale RPM and feed proportionally — chip load is
what matters. Holes ≥3×D get G83 pecking; the two deep holes are ~3.4×D.

## Verification (what was actually checked)

`verify_nc.py`: parses the G-code, simulates positions, and checks Z floors per
op, XY envelope of all cutting moves, hole-cycle coordinates vs the extracted
hole positions (exact match to 0.001), and per-section monotonic depth
schedules. **PASS on all three files.**

`wall_guard.py`: geometric keep-out check — no cutting point below facing depth
comes within tool-radius (±0.006) of the part walls (OP10), and no OP20 pocket
point violates the rail/pad/plate bodies that exist at the depth being cut
(with 0.0015 rounding tolerance). **PASS — 84 + 49 cutting points checked.**

Honest note: the generator's first draft had the finish-profile offset in the
wrong units (0.020 in instead of 0.250 in radius — it would have cut 0.23 in
into the walls) and facing lanes that missed the top edge strip. Both were
caught in review before any chip, and both checkers now guard exactly those
failure modes. Still: **air-cut or single-block the first run** — this G-code
has never run on a machine.

## When you know the machine

Tell me (or edit the obvious spots): controller dialect oddities (G83 vs G73,
tool-call format, whether T# M06 or M06 T#), max RPM, tool-holder taper + tool
length offsets, and how you zero (edge finder vs probe). If you'd rather use a
real CAM seat, the STEP file drops straight into Fusion 360 / Mastercam /
SOLIDWORKS CAM (the full suite incl. CAM is installed on this machine) — this
README's op order mirrors what you'd rebuild there.

To regenerate for a new revision: update the geometry constants at the top of
`gen_cam.py` (or point it at a new STL extraction), run it, then run both
checkers. Everything re-derives.
