# ZCode pickup: validated defense deck v18

Run `hostname` first, then read repository `AGENTS.md`. These paths are for easteregg2. The current deck pointer is `Defense/current_deck.json`.

## Deliverables

```text
D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\C3NS_Revision_20261031_v18\
├── output\
│   ├── Bolen_PhD_Defense_RoseGreen_20261031_v18.pptx
│   └── Bolen_PhD_Defense_RoseGreen_20261031_v18.pdf
├── revise_v18.mjs
├── transplant_v18.py
├── finalize.mjs
├── render_review.ps1
├── check_video_playback.ps1
├── make_contact_sheets.py
├── submitted_figs\
├── submitted_source\
└── build\
    ├── validation_delivery.json
    ├── transplant_report.json
    ├── native_video_playback.json
    ├── native_text_geometry.json
    ├── render\slide_001.png ... slide_123.png
    └── contact\contact_01.jpg ... contact_11.jpg
```

PPTX SHA-256: `d7e0b9521e94ca391ff8470fb565ff9823acbcf2d0487ff3a5c122546fde738c`.
It is 112,036,176 bytes, so it must remain Git-ignored or outside ordinary Git history.

The deck has 70 main slides, 53 hidden question slides, and a planned 33:51 main-talk duration. All 123 slides were rendered at 1600 × 900 and inspected; the PDF has 123 pages. Native charts and embedded workbooks remain on owner slides 16, 26, 31, 44, 97, and 98. Both slide-52 AnimatLab movies autoplay, loop, and remain muted. Slide 74 retains the newer Root-fixed GIF.

Slide 44 uses the committee-submitted dissertation values: torque coverage 37.35% → 65.81% and active-force coverage 93.00% → 99.31%.

The prior validated deck retained under the two-revision rule is:

```text
D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\C3NS_Revision_20261031_v17\output\Bolen_PhD_Defense_RoseGreen_20261031_v17.pptx
```

Revision 16's final tree was removed after v18 passed. Older v12–v15 folders contain no final PPTX and were not recursively destroyed.

## Sources of record

```text
D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Bolen_Dissertation.pdf
D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Bolen_Dissertation.zip
C:\Users\Ben Bolen\Documents\ME603 Dissertation\Defense\Bolen_PhD_Defense_C3NS_20261030_v13_review_AJH.pdf
C:\Users\Ben Bolen\Downloads\19_02_15 - INRC.pptx
D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\Movement and Control of Biomimetic Humanoid Robots.pptx
C:\Users\Ben Bolen\Documents\Github\ZCode_Skills\assertion-evidence\SKILL.md
```

The assertion-evidence skill now includes the advisor-linked projected-figure guidance from Rougier, Droettboom, and Bourne (2014): design one message per figure, simplify paper figures for projection, directly label evidence, use selective attention cues rather than decorative boxes, keep accessible redundant encodings, subordinate branding, and inspect at delivery scale.

## Local media

```text
C:\Users\Ben Bolen\Videos\Screen Recordings\frozen-root treadmill.mp4
C:\Users\Ben Bolen\Videos\Screen Recordings\free_drop.mp4
C:\Users\Ben Bolen\Videos\Screen Recordings\w2l_bilateralRGfreeze.gif
C:\Users\Ben Bolen\Downloads\2021-02-17-20-58-24.mp4
```

The deck is self-contained, but the external files are useful for recutting. These recordings predate the corrected body/joint geometry; do not present them as corrected-geometry evidence.

## Continuation rules

Use v18 as the baseline for v19. Preserve the C3NS rose-green master, title honeycomb, small original logo/footer, original plot colors, hidden question-slide state, native charts/workbooks, and native media. Do not flatten through PDF or screenshots. Render and inspect every slide, export a full PDF including hidden questions, run package/layout/overflow validation, and test media playback in PowerPoint. Retain only v18 and the new validated revision.

The full W2L/AnimatLab/MuJoCo continuation tree remains in `Defense/ZCODE_PICKUP_20261031.md`.
