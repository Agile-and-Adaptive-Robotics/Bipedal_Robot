# source_sync.md — Overleaf 2026-10-15 pull → ProofFinal reconciliation

**Task:** unzip `Bolen_Dissertation.zip` (2026-10-15 12:12 Overleaf pull), diff its `.tex`
files against `ProofFinal/` (chapters + main.tex), port Overleaf-only changes so ProofFinal
carries the Overleaf substance. Overleaf zip = canonical input; ProofFinal = working copy.

**Run date:** 2026-10-15. **Scratch:** `tmp/disscampaign_20261015/overleaf_pull/` (unzipped)
+ `tmp/disscampaign_20261015/diffs/` (per-file unified diffs, PF→OL so `+` = Overleaf).

**Method:** `Expand-Archive` of the zip; per-file `git diff --no-index --unified=3`
(GitHub Desktop bundled git 2.53.0.windows.4 — no git on PATH) for all 16 `.tex` + 2 `.bib`
+ `.cls` pairs; ports applied with exact-string edits; post-port verification by
whitespace-and-graphicspath-normalized equality (script `normalized_check.ps1`, result
`ALL_MATCH=True` across all 15 pairs) plus a residual-line audit (`audit_residual.ps1`)
showing the only remaining diff lines are `\graphicspath` and blank/whitespace-only lines.

## Key fact discovered: the zip carries tracked changes ALREADY APPLIED

The zip's `.tex` content is the all-tracked-changes-accepted state of the live Overleaf
project, even though `CHATGPT_REPORT.md` lists T01–T15 as *pending* (none accepted or
rejected in the review UI). Evidence: the zip contains "first **improvement** of the torque
tools" (T02's replacement), "I used these tools to convert" (T03+T04's replacement), and the
rewritten introduction (T01's old intro deleted) — while ProofFinal contained the
pre-advisor wording. Consequence: porting the zip substance into ProofFinal applies
T01–T15 content locally; the *pending* status lives only in Overleaf's review system.

## Introduction reconciliation (the ask's known divergence)

The ask stated "local 10-introduction.tex is 28 lines while Overleaf anchors comments at
lines 8/22". Measured now, **ProofFinal's intro is 50 lines / 51 LF (14,529 B) and already
carries the NEW introduction** — it is byte-identical to the zip's intro (49 LF, 14,521 B)
except for (a) the `\graphicspath` line and (b) a trailing blank line + EOF newline. The
"28 lines" premise does not match the current file; possibly it described an older local
state or a non-blank line count. Verified facts:

- Overleaf comment anchors C01 (line 8, "Artificial") and C02 (line 22, "Second") land on
  exactly the same text in both copies (intro lines are aligned 1:1 up to line 46; PF's two
  extra trailing linefeeds are at EOF only). Zip line 8 = "Artificial muscles are a key
  enabling technology…"; line 22 = "Second, actuator force alone does not determine…".
- T01 (bbolen83, 27 Sep, pending) deletes the OLD introduction (old Motivation / Project /
  Problem Statement / Objectives / Organization sections — full text preserved in
  `CHATGPT_REPORT.md` T01). That old text exists in **neither** the zip nor ProofFinal;
  both carry the replacement introduction (e.g. "The author's experience with
  \qty{10}{\mm} actuators made this gap concrete…", real `\ref{sec:robot_model}`-style
  targets instead of the old broken `\ref{sec}` stubs).
- **Action taken: none needed.** ProofFinal already carries the Overleaf substance for the
  introduction. T01's pending status is an Overleaf-review action for Ben (see Flags).

## Ports applied to ProofFinal (Overleaf-only substance)

### `ProofFinal/chapters/15-background.tex` (8 edits)
| Line (post-port) | Change | Origin |
|---|---|---|
| 38 | "first **requirement** of the torque tools" → "first **improvement** of the torque tools" | tracked change T02 (ajh26, Oct 14) |
| 44 | `festo_bpa_cad.png` `height=0.68\textheight` → `0.45\textheight` | direct Overleaf edit (figure shrink) |
| 92 | `00_GraphicalAbstract.png` `width=0.94\textwidth` → `0.75\textwidth` | direct Overleaf edit (figure shrink) |
| 110 | `aarl_bpa_biped.png` `width=\linewidth` → `width=\linewidth,height=0.52\textheight,keepaspectratio` | direct Overleaf edit (height cap) |
| 134 | "Hand placing muscle origin and insertions was a painstaking, iterative task." → "Hand placement of muscle origins and insertions was a painstaking and iterative task." | direct Overleaf edit |
| 138 | "understanding of **the** mammalian locomotor **CPG comes**" → "understanding of mammalian locomotor **CPGs come**" | direct Overleaf edit |
| 158 | "these tools were used to convert" → "**I** used these tools to convert" | tracked changes T03+T04 (ajh26) |
| 170 | "animals carry real moral cost" → "animals **also** carry real moral cost" | tracked change T05 (ajh26) |

### `ProofFinal/chapters/20-methods.tex` (2 edits)
| Line (post-port) | Change | Origin |
|---|---|---|
| 389 | "balance problem encountered by a **walking** humanoid" → "encountered by a humanoid" | tracked change T06 (ajh26) |
| 410 | "increasing foot compliance and conserving energy" → "simplifying the system" | tracked change T07 (ajh26) |

### `ProofFinal/chapters/92-AppendixA.tex` (1 edit)
| Line (post-port) | Change | Origin |
|---|---|---|
| 154 | `\begin{figure}[h]` → `\begin{figure}[htbp]` (muscle_force_compare.png figure) | direct Overleaf edit (float spec) |

### `ProofFinal/thesis.bib` (1 deletion)
- Removed the unused `@article{steele_experimental_2023, …}` entry (was lines 1931–1941).
  Verified `steele_experimental_2023` is cited **nowhere** in either tree's `.tex` files
  (grep of both trees returned zero hits; the cited Steele works are
  `steele_development_2017`, `steele_biomimetic_2018`, `steele_biomimetic_2018-1`), so the
  Overleaf deletion is a safe unused-entry cleanup.

## T01–T15 disposition (ledger: `CHATGPT_REPORT.md`, 2026-10-15)

| ID | File | Kind | Disposition |
|---|---|---|---|
| T01 | 10-introduction | Deleted (old intro) | Already reflected in both copies before this sync; nothing to port. Still *pending* on Overleaf. |
| T02 | 15-background | Changed | **Ported now** (requirement→improvement). |
| T03+T04 | 15-background | Changed+Added | **Ported now** ("these tools were"→"I"+…). |
| T05 | 15-background | Added | **Ported now** ("also"). |
| T06 | 20-methods | Deleted (" walking") | **Ported now**. |
| T07 | 20-methods | Changed | **Ported now** (→"simplifying the system"). |
| T08 | 20-methods | Deleted (grouped fragments) | Already reflected in both copies (PF 20-methods text post-dates the Sep-27 edits; diff showed no such content). |
| T09–T12 | 30-results / 40-discussion / 50-futurework / 60-conclusion | Deleted (grouped fragments) | Already reflected in both copies (only graphicspath/whitespace diffs remain in those files). |
| T13 | 92-AppendixA | Deleted fragments | Already reflected (", as they exist…"/"Fcompressive" absent from both). |
| T14 | 94-AppendixC | Deleted ("Regenerating the results files updates this table.") | Note: this sentence is still present in BOTH copies at the end of the adopted-values paragraph (94-AppendixC.tex:223 area) — see Flags. |
| T15 | main.tex | Deleted ("Creating"/"with") | Already reflected (no such content in either). |

## Divergences intentionally NOT ported (environment/cosmetic)

1. **`\graphicspath` in main.tex + all 11 chapter files:** Overleaf `figs/…` vs local
   `../Figures/…`. The Overleaf project keeps figures under a project-root `figs/` tree;
   ProofFinal builds against `Documentation/Reports and Papers/Dissertation/Figures/`
   (per-chapter subfolders — all four resized/referenced figure files verified to exist
   there). Porting `figs/` would break local builds. Keep the split; re-point only if the
   folder layouts are ever unified.
2. **Trailing whitespace / EOF newlines:** 10-introduction (PF ends with blank line + EOF
   newline; OL file has no final newline), 20-methods & 40-discussion (OL +1 trailing blank
   line), 93-AppendixB (OL +2 trailing blank lines), 94-AppendixC (OL has a trailing space
   after the ZCODE Xi-audit comment and after `\end{figure}`, +1 blank line), thesis.bib
   blank-line count (matched exactly during the port). Zero LaTeX effect; ProofFinal's
   endings kept.
3. **Zip-only artifacts:** `figs/` (69 images), `psutex.cls`/`Frontiers-Harvard.bst`
   (byte-identical to ProofFinal's), `BolenFrontiers22.bib` (byte-identical) — nothing to do.

**Byte-identical already (0-byte diffs):** 01-copyright, 02-abstract, 04-dedication,
06-ack, 90-bib, psutex.cls, BolenFrontiers22.bib.

## Verification run

- Post-port `git diff --no-index` per file → residual `+`/`-` lines are ONLY `\graphicspath`
  lines and blank/whitespace-only lines (`audit_residual.ps1` output, this session).
- Whitespace-and-graphicspath-normalized comparison of all 15 file pairs (main.tex, 11
  chapters present in both, thesis.bib, BolenFrontiers22.bib, psutex.cls):
  `ALL_MATCH=True` (`normalized_check.ps1`, this session). ProofFinal now equals the
  Overleaf zip in substance.
- Local LaTeX compile: **not run** (out of this task's scope; the ask named diff+port only,
  and the build folder is ProofFinal/build_wf). The ported edits are wording, figure-size
  options, a float specifier, and a bib-entry deletion — no new packages, refs, or files.

## Files touched this session

Modified (4): `ProofFinal/chapters/15-background.tex`, `ProofFinal/chapters/20-methods.tex`,
`ProofFinal/chapters/92-AppendixA.tex`, `ProofFinal/thesis.bib`.
Created (1): this log (folder `Notes/revision_campaign_20261015/` created).
Scratch (disposable): `tmp/disscampaign_20261015/` (unzip + diffs + check scripts).
Not touched: no `upload/` mirror exists anywhere under `Documentation/.../Dissertation/`
(searched this session — the folder is absent, so nothing to avoid or modify there);
no code, data, result, figure, or Overleaf-side files were touched, and no git operations
were run.

## Flags for Ben

1. **T01–T15 are still PENDING on Overleaf**, but the zip (and now ProofFinal) already
   carry their content applied. If you REJECT any tracked change in the Overleaf review UI,
   the Overleaf file will revert and re-diverge from ProofFinal — re-run this sync after
   any reject.
2. **T14 looks un-applied on Overleaf itself:** its deleted text ("Regenerating the results
   files updates this table.") is still present in the live zip's 94-AppendixC.tex (end of
   the adopted-values paragraph, ~line 223) — i.e., this grouped deletion entry does not
   match the current file content. Not ported either way (content identical in both
   copies); ruling on that sentence is yours.
3. **The ask's "local intro = 28 lines" premise did not hold** (actual: 50 lines, already
   the new introduction, byte-equal to Overleaf modulo path + EOF newline). If a 28-line
   introduction exists somewhere else you expected to be the working copy, point me at it —
   it is not `ProofFinal/chapters/10-introduction.tex`.
4. **Open advisor threads C01/C02 (intro) and C03+ (background onward) are untouched** —
   they are review-system comments, not file content; addressing them is the revision
   campaign's next step, not this sync's.
5. The three 15-background figure resizes and the 92-AppendixA `[h]`→`[htbp]` change were
   untracked direct edits on Overleaf (author not recorded in the ledger); ported as-is on
   the zip-is-canonical rule.
