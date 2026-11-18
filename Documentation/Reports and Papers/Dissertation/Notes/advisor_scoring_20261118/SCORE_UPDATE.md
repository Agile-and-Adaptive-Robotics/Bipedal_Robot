# Introduction feedback and cumulative score, November 18, 2026

Alex Hunt (`ajh26`) made 36 tracked edits in Chapter 1, Section 1.1, Motivation, between 8:57 and 9:21 a.m. Pacific. The original wording receives −36. The three tracked edits and two comment messages scored on November 17 contributed −6, bringing the cumulative observed score to **−42**. The [single review table](REVIEW_TABLE.html) includes the exact original and proposed text, attribution, location, individual deduction, editorial interpretation, and disposition.

| Introduction line | What the changes address | Replacements | Insertions | Deletions | Score |
|---:|---|---:|---:|---:|---:|
| 10 | Why the Festo force model needed resting-length data | 3 | 2 | 0 | −5 |
| 12 | How force reaches the knee and what the route redesign achieved | 8 | 5 | 3 | −16 |
| 14 | Why neural-controller design needs the SADb and how readers use it | 10 | 2 | 3 | −15 |

The new entries contain 21 replacements, nine insertions, and six deletions. A replacement counts once. The comma deletion and the single-letter replacement each count once. The two Chapter 4 corrections were Ben’s edits. His decision to retain line 451 is recorded in the review table.

### Editorial findings

- **BPA force model.** Specify the missing relation between resting length and output. Connect the bench measurements to the force model they made possible.
- **Knee torque.** Explain how a muscle pulling on a limb creates a joint moment before introducing brackets, tendons, and bending. Connect the observed prediction error to the rig, the redesigned routes, and the human torque criterion.
- **Neural control and the SADb.** Start with the human-like controller the robot needs. The dispersed evidence motivates the database. Name who built it and what a reader can find or do with it.

### Evidence and scoring method

The chapter had no entries in the November 17 Review overview. The November 18 overview had 423 entries across the project, including 36 in `10-introduction.tex`. All 36 were authored by `ajh26`. Each change was read in full in the open editor. Replaying its exact strings and offsets against the November 17 text reproduced all three visible proposed paragraphs. The [Introduction capture](LIVE_INTRODUCTION_EDITS.json) preserves the original and proposed strings, including leading and trailing spaces.

The [cumulative classifier ledger](CLASSIFIER_RESULTS.json) records the model identity, prompt, rubric, individual judgments, original and proposed paragraphs, and links among related edits. The language-model classifications group the 36 edits into six editorial issues. Each deduction follows Ben’s fixed per-edit rule and is charged to the original draft. The 130 previously flagged unmarked passages remain predictions. Publication and review decisions for the tracked changes remain with Ben. The [scoring script](update_scores.py) reproduces the ledger and table and checks the live-text anchors.

