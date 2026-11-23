# Reviewing dissertation prose against Ben’s editing history

Ben’s September machine-teaching ledger paired each judgment with the words that earned it. Its two positive examples show an author telling readers where a model will be used and why an optimization is being introduced. Ben’s September 12 grade of B− concerned the coherence of the combined papers. His October 22 clarification removed credit for sections that drew no comments. [September ledger preserved at 5bf8e762](https://github.com/Agile-and-Adaptive-Robotics/Bipedal_Robot/blob/5bf8e7620441b40161ddd53f241f2b2d4c54552a/Documentation/Reports%20and%20Papers/Dissertation/summary.md) · [October scoring instructions](https://github.com/Agile-and-Adaptive-Robotics/Bipedal_Robot/blob/5bf8e7620441b40161ddd53f241f2b2d4c54552a/Documentation/Reports%20and%20Papers/Dissertation/Notes/advisor_revision_20261022/ROUND4_NOTES.md)

The October 23 review recorded 0/4 for three **previous attempts at revision**. The advisor still could not follow the database organization (M-065) or two passages on the SNS work (M-106 and M-122). The other 134 cards were unassessed at capture. The October 30 audit counted five distinct assistant mistakes with a weighted cost of 115, including missing a tracked deletion and attributing a task to the wrong person. Those weights and the older campaign rewards retain their original accounting. [October readiness rules and matrix](../../../../../CHATGPT_REPORT.md) · [October 30 audit rubric](https://github.com/Agile-and-Adaptive-Robotics/Bipedal_Robot/blob/f99106cfce46cb64c86a66d6429c9de668586cbf/Documentation/Reports%20and%20Papers/Dissertation/Notes/publication_audit_20261030/REVIEW_RUBRIC.md)

## Use on a new chapter draft

Ben reserves the ethical argument in Chapter 2, “Testing Neuromechanical Hypotheses on Robots: The Ethical Case,” for his own judgment. The advisor-style reader skips that section. Use Ben’s exact wording for his position and verify quotations, citations, and factual descriptions independently when he requests an edit.

Run the following command from the repository root on easteregg2. It uses the preserved pre-November advisor signals and the dated November ledgers. Supplying chapter paths after the options selects different TeX files.

```powershell
& 'D:\Anaconda\python.exe' 'Documentation/Reports and Papers/Dissertation/Notes/advisor_feedback_history/review.py' --output 'tmp/dissertation_history_review.json'
```

The output contains the current passages and nearby historical issues found by word and phrase similarity. It also gives actual advisor examples from different months for each recurring reader question and includes the [praised examples, Ben’s rulings, and earlier editorial minuses](EDITORIAL_EXAMPLES.json). Read the marked words, the reply or tracked edit, and the source. A topic match can point to an unrelated criticism. A reader can recognize a recurring presentation problem even when the subject and vocabulary have changed.

Read the surrounding section, then answer the historical question from the paragraph. For a mechanism, identify what acts on what. For a contribution, identify who built it and how a colleague uses it. For an experiment, name the change, controls, readout, and warranted conclusion. Preserve passages that already give a concrete answer. For a recurrence, record the exact draft, the linked historical examples, the evidence, the revised passage, and the independent reader’s answer. The [November application](../chapters_6_7_20261118/PASSAGE_REVIEW.json) includes technical corrections, improvements in presentation, a consequential consistency edit, and passages retained after reading.

To verify an existing review, add the following argument to the command. It checks the exact text even if the passage has moved. A later modification makes the previous local judgment stale. A newly captured criticism anchored to the same passage reopens the judgment even if the wording is identical. New feedback with ambiguous anchors is listed for inspection.

```powershell
--decisions 'Documentation/Reports and Papers/Dissertation/Notes/chapters_6_7_20261118/PASSAGE_REVIEW.json'
```

## How to interpret a score

Ben’s November 19 toolchain corrections reopen the names and uses of the bipedal controllers, wording for suspended stepping, the joint-design heading, and slide-to-manuscript figure choices. Check the actual user action in each interface and the configuration used in each experiment. The five author examples travel with each new historical review. [Latest passage judgments and evidence](../toolchain_sadb_20261119/PASSAGE_REVIEW.json)

| Evidence | Accounting |
|---|---|
| Attributable advisor tracked change | −1 against the original wording for each distinct insertion, deletion, or replacement. |
| Attributable advisor comment | historically 0 for pure praise, −1 for wording or clarification, −2 for a substantive revision, or −3 for a factual or interpretive error. Mixed praise and a request receive the severity of the request. |
| Several edits or messages on one concern | Preserve every distinct signal and link the underlying issue. Recapturing an event earns no further deduction. |
| Present draft resembles a historical criticism | A reviewer may assign an inspection priority and explain the parallel. This adds zero observed advisor events. Silence leaves a passage unlabeled. |
| A specific revision has been checked | Use the October 0–4 scale. A repeated failure reopens the prior judgment. A 3 requires source and contextual verification. A 4 requires Ben’s or the advisor’s confirmation. |

The records preserve authors, event dates, exact wording when available, and qualified anchors. User directives, replies from Ben, other readers, and defense-slide feedback retain their own provenance. The November ledger preserves the original strings and individual scores for the three November 17 edits and two comments, and all 36 November 18 Introduction edits. [Historical extraction and exclusions](HISTORICAL_SIGNALS.json) · [November advisor ledger](../advisor_scoring_20261118/CLASSIFIER_RESULTS.json)

For subsequent captures, preserve a stable event identity and link the source before extending the ledger. The checker reads dated `advisor_scoring_*/CLASSIFIER_RESULTS.json` captures, ignores repeated records with the same identity, and flags incompatible recaptures. New issue categories need their editorial tags in the capture or the scorer. The `--before YYYY-MM-DD` option restricts evidence to earlier event dates when testing whether an older record anticipates later feedback.

## Figures, tables, and conclusions

The November 20 Chapter 5 review reopened the result-first judgment. The prior check retrieved examples paragraph by paragraph but missed that an entire Results section opened with inventory. The checker now extracts the first prose under a Results heading and its subheadings, alongside three attributed Chapter 5 comments and the advisor’s before-and-after Introduction edit `INTRO-T36`. Read each opening with its example and write down what a researcher can learn. The historical torque correction asks the same question about a physical design and its human reference.

Add `--check-results-openings` to flag openings that lead with a record inventory, instructions, or a roadmap for contextual review. The report also lists headings without opening prose. A flag requires inspection, and an unflagged opening still needs a reader’s judgment. The regression test uses the rejected Chapter 5 lead, its procedural example, and a quantitative torque result that answers a research question. [November 20 review and readback](../submission_review_20261120/REVIEW.md)

The November 19 Overleaf review added six Chapter 6 comments and five tracked edits from the advisor (−15), and one Introduction comment and three tracked edits (−5). Across the November 17–19 captures, 56 distinct advisor events total −62. Three November 19 edits by Ben are recorded separately. The reviewer now returns the first-person phrases in each passage and distinguishes ownership, narration, and tool capability for contextual reading. The two filler comments and the correction from “our” to “my” belong to different editorial judgments. [Dated review table](../ch6_advisor_20261119/REVIEW_TABLE.html) · [New advisor events](../advisor_scoring_20261119/CLASSIFIER_RESULTS.json)

Ben’s November 18 readback added four author examples to the teaching ledger. For a stimulus and response, compare their time coordinates in the printed figure. For a model table, check row order against first introduction. At the end of a chapter, name the experiment a future researcher can run. For Chapter 6, start with the AnimatLab hand-tuning problem and identify the specific timed implementation when discussing the gain in speed. The model retrieves these examples on every pass. [Author feedback and independent reread](../ch6_sequence_20261118/REVIEW_TABLE.html)

An October revision spelled out Walker Two-Layer at first use. A later rewrite lost the definition, and Ben raised the same problem in November. For each new chapter draft, read the model and program names in printed order, including figure and table labels. Ask what each name denotes and what distinguishes each variant. The Chapter 6 check records the source line and printed page for 36 terms, including W2L, W2Lvar, BiW2Lvar+, W2L108, s3k, and syn6. [Revision and evidence](../ch6_first_use_20261118/FIRST_USE_REVIEW.md) · [Chapter 6 first-use check](../ch6_first_use_20261118/verify_first_use.py)


## November 20 afternoon review

At Ben’s direction, explicit praise of an identifiable revision now earns +1. The earlier mixed comment on the Motivation paragraph and its deduction remain in the historical ledger. Hunt’s subsequent “Much better” earns +1 on the revised paragraph. Ben’s positive assessment of the amount of Chapter 5 retained is recorded separately. No passage acquires advisor approval solely because it remained in the document.

The reviewer now flags three questions across prose in methods, results, and appendices. Does a statement about checking name the correction or measured agreement? Does a possible import identify who performed it and which program received the model? Does an experiment introduce its scientific purpose before numerical settings? The flags direct a reader to the surrounding evidence. On five documented before-and-after pairs, the checks caught four original problems; the unexplained return to McCrea still required a reader. Seven short, separately worded examples exercise these distinctions in the regression suite. [Attributed scoring events](../advisor_scoring_20261120/CLASSIFIER_RESULTS.json) · [Before-and-after cases, missed example, and remaining flags](../advisor_feedback_20261120/VERIFICATION.json)
