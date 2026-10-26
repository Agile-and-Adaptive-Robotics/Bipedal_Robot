# Local prose revision D — 23 October 2026

**18 targeted edits across Introduction, Background, Database, and SNS are complete locally. Overleaf is unchanged.** Ben clarified that ChatGPT Astra is his intended final prose model. No model was launched, no provider was configured, and no scheduled monitoring was restarted.

Working copy: [Revision_D_local.pdf](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/Revision_D_local.pdf>). Actionable code/model handoff: [ZCODE_HANDOFF.md](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/ZCODE_HANDOFF.md>). Exact before/after changes: [changes.diff](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/review/changes.diff>); per-edit record: [changes.json](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/review/changes.json>).

## Source and scope

The five captured Overleaf chapter snapshots supplied the reviewed source. Uncaptured chapters, main.tex, bibliography, and initial assets came from the existing Revision_C copy solely to provide local build support. This is **not a fully reconciled current Overleaf project** and is not an upload bundle. ZCode must merge the targeted changes into fresh source, not replace live chapters wholesale. ProofFinal and Revision_C were not edited by this pass.

The existing flexor panel Flx20mm_meetsHuman.pdf was missing from the Revision_C support assets. It was copied unchanged from tmp/overleaf_live_20261023 for the build. No figure was generated or replaced. The database pivot caption's “animal preparations” was explicitly changed to “experimental preparations”; all other figure blocks were preserved, and all figure inclusions, labels, and citation commands remain conserved.

## Changes and self-check

| Edit | File and final source line | Change / reason | Local self-score |
|---|---|---|---|
| D-01 | [10-introduction.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/10-introduction.tex:12>) | Correct subject/verb agreement and tense. | 3/4; acceptance unverified |
| D-02 | [10-introduction.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/10-introduction.tex:12>) | Make the dispersed findings the explicit plural subject. | 3/4; acceptance unverified |
| D-03 | [10-introduction.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/10-introduction.tex:12>) | Remove duplicate text and misspelling; scope the missing combination to this project instead of asserting an unverified field-wide first. | 3/4; acceptance unverified |
| D-04 | [50-database.tex:10](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/50-database.tex:10>) | Define archived by its actual app-only status; source sadb_app.html archive control/header. | 3/4; acceptance unverified |
| D-05 | [50-database.tex:10](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/50-database.tex:10>) | Explain citation direction; verified build_citation_graph.py:8 and build_app.py:320-322,725-744. | 3/4; acceptance unverified |
| D-06 | [50-database.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/50-database.tex:12>) | Explain record types and relationships in reading order; preserve existing example and explicitly avoid claiming complete links. | 2/4; acceptance unverified |
| D-07 | [50-database.tex:25](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/50-database.tex:25>) | Explain the purpose of controlled vocabulary without inventing category definitions. | 2/4; acceptance unverified |
| D-08 | [50-database.tex:40](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/50-database.tex:40>) | Use the existing preparation category without implying that every category is an animal species. | 2/4; acceptance unverified |
| D-09 | [50-database.tex:52](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/50-database.tex:52>) | Answer the explicit future-work request without claiming all requested links are absent. | 2/4; acceptance unverified |
| D-10 | [60-sns.tex:14](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/60-sns.tex:14>) | Replace abstract missing-step language with observed outcome, limitation, and specific requirement; the existing 92-actuator value remains later in this chapter. | 2/4; acceptance unverified |
| D-11 | [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/60-sns.tex:27>) | Keep the advisor five-part structure but distinguish tasks from completed systems, remove tested-it ellipsis, and define graphical construction. | 2/4; acceptance unverified |
| D-12 | [60-sns.tex:54](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/60-sns.tex:54>) | Remove literal ellipsis from table cell; preserve all capabilities and limits. | 3/4; acceptance unverified |
| D-13 | [60-sns.tex:61](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/60-sns.tex:61>) | Define checks rather than implying unspecified audit discoveries; actual defects and repairs remain a ZCode evidence request. | 2/4; acceptance unverified |
| D-14 | [60-sns.tex:63](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/60-sns.tex:63>) | Define BPA interface inputs, force output, feedback, and limits at first mention without changing model functionality. | 2/4; acceptance unverified |
| D-15 | [15-background.tex:224](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/15-background.tex:224>) | Complete a sentence fragment with the existing scoped attribution. | 3/4; acceptance unverified |
| D-16 | [15-background.tex:226](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/15-background.tex:226>) | Remove unsupported overnight throughput/statistical-power assertions; explicitly logged, not silently harmonized. No numerical results or figures removed. | 2/4; acceptance unverified |
| D-17 | [15-background.tex:230](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/15-background.tex:230>) | Clarify repeatability versus biological replication; retain Three Rs definition/citations and the existing figure. | 2/4; acceptance unverified |
| D-18 | [50-database.tex:8](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/50-database.tex:8>) | State the project-specific evidence gap in a direct sentence. | 2/4; acceptance unverified |

Scores of 3 concern a directly checked local clarification; scores of 2 retain an unresolved contextual, reader-understanding, or implementation-evidence dependency. No score of 4 was assigned. Earlier 0/4 clarity failures remain in the historical matrix; these local revisions do not erase them. Context IDs in changes.json link nearby matrix entries but do not claim every nearby short comment is addressed.

## Verification

- Fresh pdflatex → bibtex → pdflatex ×2: **257 pages, zero errors, zero undefined references, zero missing figures.** Full logs and build_result.json are in review/ and the named build directory.
- Figure paths/counts, label lists, and citation command lists are unchanged in all four edited files. Bibliography was copied byte-for-byte and never re-encoded.
- All numerical measured/model results are preserved. In SNS, duplicated mentions of the 92-actuator count and Gait2392 name were shortened; their substantive values remain. Literal list ordinals 1–5 became an enumerate list. Numeric-token differences are logged in conservation.json; no result value was silently changed.
- Background's unsupported claims about thousands of overnight combinations, hundreds of cycles per condition, unmatched throughput, and statistical power were explicitly replaced by a qualitative statement of repeatable simulation tests. These were unverified performance claims, not benchmark results. The assertion that identical models are independent biological subjects was explicitly removed. Old wording is retained in changes.json/diff.
- Citation directions checked against SADb_audit/build_citation_graph.py:8 and SADb_audit/app/build_app.py:320–322,725–744. App-only archive status checked against the app header/archive control. These are local code checks, not claims of an exhaustive live-base schema audit.
- Database opening and SNS task-list page renders were inspected: text was legible, with no clipping on those sampled pages. **26 overfull boxes remain in the full build.** The document is not declared layout-clean.

## Rubric disposition and unresolved work

Gate 1: no new external prose copied; existing citation commands and publication notes preserved. Gate 2: changed prose checked for banned boilerplate, unexplained labels, and reviewer-comment leakage; no new reviewer text or agent comments placed in .tex. Gate 3: numerical results conserved, scoped claims checked as above; detailed implementation claims still require ZCode's artifact audit. Gate 4: direct clarifications made; broad clarity concerns remain provisional until Ben/advisor reads the result. **No overall scientific-validation GO is claimed.**

ZCode's handoff specifies evidence to retrieve for actual conversion defects, five-system integration status, bilateral organization/authorship, source-model replication, earlier-library comparison, and database links. It requires existing evidence before new runs, source-backed edge tables before optimization, and one baseline plus one published stimulus/deletion protocol before expanding a replication campaign. GLM-5.3 / Flash are Ben's preferred execution models there. ChatGPT Astra is reserved for bounded final prose polishing after evidence is settled.

All source edits are local and reviewable. No Overleaf changes, uploads, replies, accept/reject actions, resolutions, commits, or pushes were performed.

---

