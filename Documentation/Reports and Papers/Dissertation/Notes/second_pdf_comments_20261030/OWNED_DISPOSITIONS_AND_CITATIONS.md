# Owned-source disposition and citation provenance

Only source files `10-introduction.tex`, `15-background.tex`, `50-database.tex`, and `96-AppendixE.tex` were edited by this owner (Chapters 1, 2, 5, and Appendix E). Root integrates the bibliography and clean builds. No figure/data files, the published muscle chapter, or database records were edited.

| Directive | Current source evidence | Disposition |
|---|---|---|
| Expanded goals | `10-introduction.tex:6`; `10-introduction.tex:8` | Biology informs autonomous robotic movement in unstructured terrain; controllable robots return questions to biomechanics/neuroscience; prosthetic/orthotic design is a further translational goal. Energy economy is explicitly a potential motivation, with no claim of measured net savings. |
| Knowledge gaps and four contributions | `10-introduction.tex:10`; `10-introduction.tex:12`; `10-introduction.tex:16` | State tested-Festo force-prediction and robot-specific transmission gaps, followed by sensory-evidence access and controller-integration gaps. Preserve the four ordered objectives: BPA characterization, isometric torque, SADb, and neuromechanical tools. Remove the introduction's incorrect BPA-humanoid classification and use musculoskeletal arms. |
| Primary experimental examples | `15-background.tex:158`; `15-background.tex:189`; `15-background.tex:191` | Add the isolated stick-insect preparation, cat group I stimulation, and distinct hip-position phase-transition experiments, retaining preparation and step-phase qualifications. |
| Primary model/posture examples | `50-database.tex:76`; `50-database.tex:78`; `50-database.tex:80` | Explain cat postural-synergy evidence, cat split-belt/model comparisons, reflex and CPG neuromechanical models, and the salamander robot's testable biological prediction. Each example states a finding rather than using an unexplained citation train. |
| Three-review provenance | `50-database.tex:6`; `50-database.tex:21` | Remove the undocumented statement that three particular reviews guided or inspired the collection. Retain the documented Zotero/Airtable origins and backward/forward citation-search process. Reviews remain cited for their actual background role in Chapter 2. |
| C3 | `96-AppendixE.tex:12` | Remove only the claim that 38 retained records lack grounding text. Retain 550 DOIs and every other count. No database, Zotero, or Airtable record was deleted or modified. |
| C4 | `50-database.tex:44`; `50-database.tex:46` | Explain spring attraction/repulsion, distinguish network coordinates from year/citation axes, define the selected focus paper, and identify the illustrated Grillner (1985) focus-paper neighborhood. Preserve 57/13/555 and all callouts. |
| C11/C18 | `50-database.tex:74`; `96-AppendixE.tex:20` | Call Table E1 for afferent classes/receptors/signals and Table E2 for phase-dependent pathways and their source studies. Both existing table bodies are unchanged. |
| Authorized responsibility sentence | `50-database.tex:67` | Use “I made the final inclusion, classification, and removal decisions.” All preceding ratings/comparison wording remains byte-identical; no formal or documented assessment method is invented. |
| C19 | `30-torque.tex:342–348` (read-only) | The passage describes lumbosacral actuator-placement searches, Cartesian perturbations, numerical weights, and the later knee mesh search. “Spine anatomy” was a misleading audit label. No source deletion or change to Chapter 4 was made by this owner. |
| Concrete layout repairs | `10-introduction.tex`; `15-background.tex`; `50-database.tex` | Move the exact Appendix-roadmap paragraph before Figure 1.1. Remove only the section barriers before Sensory Feedback, Computational Neuroscience, and Curation Methods to prevent sparse continuation pages. Retain the mechanical-before-CPG barrier and every graphic dimension/content. |

## Primary sources selected from the actual AARL library

All 7 new references also match records in the existing 1,393-entry SADb export. Six match by DOI; Ijspeert (2007) matches by title because its stored DOI is blank. No ranking by citation counts or historical inspiration is claimed.

| New key | AARL item | Primary evidence | Requested-author connection |
|---|---|---|---|
| `grillner_initiation_1978` | `HCRQPHCM` | [10.1016/0006-8993(78)90973-3](https://pubmed.ncbi.nlm.nih.gov/274169/) | Grillner |
| `mcvea_role_2005` | `HBAEC9DM` | [10.1152/jn.00511.2005](https://pubmed.ncbi.nlm.nih.gov/16093331/) | McVea, Pearson |
| `guertin_ankle_1995` | `J347PRR7` | [10.1113/jphysiol.1995.sp020871](https://physoc.onlinelibrary.wiley.com/doi/10.1113/jphysiol.1995.sp020871) | Guertin, Perreault |
| `buschges_rhythmic_1995` | `CK2NF8VD` | [10.1242/jeb.198.2.435](https://pubmed.ncbi.nlm.nih.gov/9318078/) | Buschges |
| `latash_organization_2020` | `ZXDCXIBI` | [10.3389/fnins.2020.598888](https://www.frontiersin.org/journals/neuroscience/articles/10.3389/fnins.2020.598888/full) | Frigon, Rybak |
| `ijspeert_swimming_2007` | `WNZGBKM2` | [10.1126/science.1138353](https://www.epfl.ch/labs/biorob/research/amphibious/salamandra/) | Ijspeert |
| `ting_limited_2005` | `NF3C56HN` | [10.1152/jn.00681.2004](https://pubmed.ncbi.nlm.nih.gov/15342720/) | Ting |

Existing primary citations supply Akay's mouse feedback experiments, Geyer and Herr's reflex model, Rybak's two-level model, and Deng (2022) (coauthors include Szczecinski, Chiel, and Quinn). The uncertain name “Gerr” is covered cautiously by the existing Geyer/Herr citation; no author identity was invented in bibliographic metadata.

The spring-layout explanation is supported by `SADb_audit/build_citation_graph.py:144` and the [official NetworkX spring-layout documentation](https://networkx.org/documentation/stable/reference/generated/networkx.drawing.layout.spring_layout.html).

## Preservation and verification

29/29 source-conservation checks pass. They include the untouched force-model comparison/coefficients, complete Chapter 2 ethics/film section, four objectives, all graphic calls/dimensions, the two Appendix E table bodies, corpus counts, and excluded AI comparison. The only removed count is the expressly withdrawn 38-record interpretation. Detailed checks, source hashes, before copies, and diffs accompany this ledger.

Initial r1 visual inspection identified sparse continuation pages at PDF pages 30, 48, 53, and 143; the exact paragraph relocation and three local barrier removals address those defects. Actual stable r2 visual verification passes at PDF pages 28, 46, 49, and 139; Figure 5.1 remains early on PDF 140, the third page of Chapter 5. The final PDF has 297 pages and SHA256 `af0b89130342a9d8383b523924782125473da1d931a6c79d569925aee25db9c4`. These are direct checks of the compiled PDF.

## Read-only curation flag retained for later work

Ting (2005) record `rec6Vhp6wlxjLUCWZ` currently has no animal tag and its note calls it a computational-model study. The primary article reports cat postural experiments; the revised manuscript accurately describes that experiment. This metadata discrepancy was reported, and the database record remains unchanged. The Ijspeert (2007) record `recyFbFcjkWSlqlMU` lacks its DOI in the export; the DOI is verified for the manuscript citation without changing the app.
