# Style contract — Ben's dissertation voice

Campaign: post-defense revision, 2026-10-15. Governs **all new and revised prose** in
`Documentation/Reports and Papers/Dissertation/` (working set: `ProofFinal/chapters/*.tex`).
Untouched existing text is grandfathered, but a revision that touches a paragraph makes the
whole paragraph subject to this contract. Paths below are relative to the repo root
(`D:\Github\Bipedal_Robot` on this machine; same relative layout on the other clones).

Every rule here is grounded in the existing chapters; where a rule and an existing passage
disagree, flag it in the revision report instead of silently rewriting grandfathered text.

## 1. Register and audience

- Audience: a mechanical-engineering faculty committee. Professional academic prose,
  technically dense, direct declarative sentences. Say the thing; support it; move on.
- Exemplars of the voice — calibrate new prose against these:
  - Framing/argument: `ProofFinal/chapters/10-introduction.tex:6-24`.
  - Methods rationale (why a correction is formulated, not derived):
    `ProofFinal/chapters/20-methods.tex:161-163`.
  - Results interpretation with honest hedging: `ProofFinal/chapters/40-discussion.tex:9-17`.
- Hedge only where uncertainty is real ("It remains to be seen how $\epsilon_{620}$ can be
  known a priori", `40-discussion.tex:17`). Never hedge to soften a claim the data support.

## 2. Person and tense

- **SUPERSEDED 2026-10-22 (Ben):** first person "I/my/me" replaces "the author"
  everywhere in the dissertation chapters (the advisor's active-review request).
  "We/our" remains for work done with collaborators. Never "the researcher" as a
  first-person dodge. Same date: self-plagiarism gating against Ben's own papers
  (published or in-preparation) is REMOVED — chapters may follow his own papers
  closely, and his own paper is cited wherever one exists for the content.
- Past tense for procedures and data collection (`20-methods.tex:80,131`); present tense
  for model behavior, equations, and interpretation (`20-methods.tex:106`; `40-discussion.tex:15`).

## 3. Punctuation

- **Oxford comma, always**, in every list of three or more items
  (`10-introduction.tex:6`; `20-methods.tex:94`).
- **Em-dash policy.** The existing text uses spaced `---` parentheticals rarely (one in
  `10-introduction.tex:24`). New prose: at most one em-dash per ~500 words, never to splice
  independent clauses habitually; prefer a comma, parenthesis, or colon. Em-dashes as
  connective padding are an AI-slop tell (rubric Gate 2).
- En-dash `--` for numeric and conceptual ranges is correct and exempt
  (`\SIrange[range-phrase=--]{112}{518}{\mm}`, "pressure--contraction pairs",
  `20-methods.tex:80`). Do not "fix" these to hyphens or em-dashes.

## 4. Spelling

- American English. Verified: zero hits for `behaviour|fibre|modelled|colour|centre` across
  `ProofFinal/chapters/*.tex`; "fiber" is used (`10-introduction.tex:8`). Keep it that way;
  do not introduce British forms via quoted material without `[sic]` or paraphrase.

## 5. Numbers and units

- Units go through siunitx: `\qty{620}{\kPa}`, `\qtylist{10;20}{\mm}`,
  `\SIrange{...}`. No bare `10mm` in prose. This is existing practice
  (`20-methods.tex:80,104,133`) — keep it.
- **Every number needs a source.** A number in dissertation prose is one of: a measured
  value (traceable to a `Testing_Data/` file or results table), a computed value (traceable
  to the generating script/mat), or a cited value (traceable to the reference, at the cited
  page for datasheets — see the page-numbered Festo citations in `40-discussion.tex:13`).
  If you cannot name the source, the number does not go in. The supervisor enforces this
  (rubric Gate 3).
- Never "harmonize" two disagreeing numbers by editing one to match the other. Discrepancies
  get investigated or footnoted, not silently resolved.

## 6. Citations

- natbib `\citet`/`\citep` as the existing chapters do. Factual claims about others' work
  carry a citation at the claim, not paragraphs later.
- Text adapted from Ben's own two papers follows the adaptation rules in the rubric
  (Gate 1) and keeps the publication note at `40-discussion.tex:5` truthful.

## 7. Banned AI-slop patterns

Never appear in new or revised text (regexes and tiers in the rubric, Gate 2). Headline
bans: "delve", "it is important to note", "plays a crucial role", "landscape" as metaphor,
"moreover/furthermore/additionally" chains, em-dash padding, uniform paragraph rhythm.
Tier A is a hard fail; Tier B is judged in context.

## 8. Rhythm and paragraphing

- Vary paragraph length and opening structure. Long analytical paragraphs are fine when the
  content earns them (`20-methods.tex:163` is one long paragraph and correctly so); a run of
  same-shaped 3–4-sentence paragraphs is not.
- Sentence length varies too. If three consecutive sentences share the same
  subject-verb-object skeleton, restructure one.

## 9. Reviewer-comment etiquette

- **Never paste reviewer or committee text into the dissertation.** Comments are responded
  to by changing the document, not by quoting the comment, and not with meta-commentary
  ("as requested by the committee…"). The dissertation reads as if it always knew.
- The supervisor judges whether an edit actually answers a comment, against the comment text
  (rubric Gate 4).

## 10. Change discipline

- An edit answers its comment and stops there. Drive-by rewording of untouched sentences is
  scope creep and gets flagged.
- Math edits additionally follow `mr-nomenclature.md` (same folder).
- Every revision batch reports: comment quoted, files touched with line ranges, numbers
  verified with sources, gates run.
