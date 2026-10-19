# ChatGPT handoff — defense-deck external review (Bipedal_Robot, Oct 2026)

This file briefs any AI assistant (ChatGPT or otherwise) working on this project while
the primary assistant (ZCode) is unavailable. Read it fully before doing anything. The
repo-wide `AGENTS.md` next to this file has additional standing context — both are plain
markdown. Previous handoff entries were erased on Ben's explicit instruction (2026-10-18);
their content is preserved in this file's git history if ever needed.

## 2026-10-18 — REVIEW THE DISSERTATION DEFENSE DECK (ChatGPT task, from Ben)

**Machine: easteregg2.** Branch `KneeTestSetup_BenBo_stw`, working tree left clean for
Ben to commit via GitHub Desktop — do NOT commit or push anything yourself.

### What to review

The defense deck:
`Documentation/Reports and Papers/Dissertation/Defense/Bolen_Dissertation_Defense.pptx`
(54 slides; defense is 2026-11-06, dissertation submitted 2026-10-21). Fresh per-slide
PNG renders live in `Documentation/Reports and Papers/Dissertation/Defense/deck_builder/
deck_render/slide_01.png … slide_54.png` — use those for judging appearance without
opening PowerPoint.

### Task, in Ben's words and order

1. **Evaluate the deck's adherence to the assertion-evidence approach** — the
   Alley/Penn State method. The references of record are
   https://www.assertion-evidence.com/ (start here; its Templates and Tutorial links
   point to craftscicom.org — `ae_templates.html`, `ae_tutorial.html`, and the official
   template files under `/uploads/5/6/1/4/56145985/`), plus
   https://www.writing.engr.psu.edu/teaching_presentations.html and
   https://ulrikeboehm.org/assertion-evidence-presentations/. The local standard,
   including the ground-truth geometry extracted from the official 16:9 template, is
   codified in the `assertion-evidence` skill
   (`C:\Users\Ben Bolen\Documents\GitHub\ZCode_Skills\assertion-evidence\SKILL.md`) —
   read it first. Key proportion rules on this deck's 13.333×7.5 in canvas: assertion
   headline 34 pt bold near edge-to-edge at the top (≤ 2 lines, a full-sentence claim
   with the result, ideally the number); visual evidence filling most of the slide
   below; at most 1–2 short call-outs; no bullet lists; prose lives in speaker notes;
   nothing below ~11 pt.
2. **Add or revise figures** where the evidence is weak: a slide whose visual does not
   carry its assertion, a chart too dense to project, or a claim with no on-slide
   visual. Figure sources live in `Documentation/Reports and Papers/Dissertation/
   Figures/<chapter>\` trees and the deck's own `deck_builder/deck_assets/`. The
   pipeline is Python (`deck_builder/build_deck.py` rebuilds from `slides_a-d.py` +
   `deck_lib.py`; `render_deck.py` re-exports PNGs; `audit_ae.py` is the programmatic
   A-E audit; `convert_figs.py` regenerates figure assets from the chapter PDFs).
   python-pptx runs under `D:/Anaconda/python.exe`. Edit the builder scripts, not the
   pptx by hand, so the deck stays regenerable. Keep every number exactly as cited —
   each slide's source line and the speaker notes carry the provenance; if a figure
   you add needs a number, take it from the cited chapter table or leave a TODO for
   Ben rather than inventing one.
3. **Write your findings to `CHATGPT_REPORT.md`** (repo root, same folder as this
   file): per-slide findings (slide number, what fails the A-E bar, what you changed
   or recommend), then a summary of the deck-level verdict. Ben and ZCode read that
   report on the next pass.

### Constraints

- The deck asserts numbers from the dissertation's chapters and appendices — never
  change a number; if one looks wrong, flag it in the report instead.
- "Model", never "plant"; Oxford comma; no em-dash-heavy prose; no dates or campaign
  provenance on slides; no emoji.
- Working tree stays clean for Ben's GitHub Desktop commit — no commits, no pushes.
- This repo is public: no credentials, license keys, or machine-access details in
  anything you write.

## Next step for whoever picks this up (ZCode or ChatGPT)

Read `CHATGPT_REPORT.md` when it appears, verify its claims against the rendered
slides (`deck_builder/deck_render/`), and fold accepted findings into
`deck_builder/slides_a-d.py` + `deck_lib.py`, then rebuild + re-render + re-audit
(`build_deck.py`, `render_deck.py`, `audit_ae.py`) before Ben's next look.
