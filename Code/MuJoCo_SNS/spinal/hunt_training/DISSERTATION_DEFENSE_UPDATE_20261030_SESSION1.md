# Session 1: dissertation and defense handoff, 30 October 2026

Machine: easteregg2. Research record: `TRAINING_RECORD.md` in this folder.
This note supplies material to the separate dissertation session; it does not
revise manuscript sources or the defense deck.

## Hunt source correction supplied by Ben

Hunt, A. J. (2019, February 15). *Real-time control of dynamic walking in a
dog-like robot*. Presented at the 2nd Intel Neuromorphic Research Community
(INRC2) Meeting.

Short slide credit: A. J. Hunt, INRC2 Meeting, Intel Labs (Feb 15, 2019).

Treat this as a professional-meeting lecture/presentation, not a journal
article. Local presentation identified by Ben:
`C:\Users\Ben Bolen\Downloads\19_02_15 - INRC.pptx`.
Use the actual title above instead of a placeholder. Distinguish the four
training operations from our five physical validation stages.

## Results suitable for a bounded methods demonstration

| Validation | Current evidence | Essential qualification |
| --- | --- | --- |
| 1: deafferented air | Pass; 3/3 knee cycles, period 2.467 s | Fully supported, no ground contact |
| 2: afferented air | Same saved weights pass; 7/7 cycles, period 1.214 s | Fully supported; minimum ankle-origin width 0.1086 m |
| 3: standing | Same weights pass assisted standing | About 55% BW on feet and 45% vertical spring support; joint restraints remain |
| 4: supported ground stepping | Not passed | Bilateral contact events alone do not establish alternating steps |
| 5: refinement/support reduction | Withheld | Requires stage-4 acceptance |

These are the Gait2392 human model with the experimental four-cell RF/hamstring
PF controller; they are not new AH001 robot results. Sartorius and gracilis
remain separate. Removing the added Peterka/Ting-inspired module also passes
assisted standing, so no demonstrated benefit from that module is claimed.
Six scalar perturbation checks pass; only one of two simultaneous synaptic
jitters passes. Do not describe this as broadly robust or independent walking.

The selective stance-recruitment experiment reaches about 89% BW on the feet,
but fails standing balance (22.19° tilt against a 15° limit). This is a useful
negative result/diagnostic, not a support-reduction success.

## Defense visual

Candidate folder:
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\runs\tonic_abduction_20261030_s1_v1\0p02`.

Its `air_pair\preview\stage2_air_neural.gif` and `.png` replay actual recorded
states at 1× speed with an explicit full-support caption. `CHECKPOINT.json`
contains hashes and assessment provenance. No corrected poses or teacher
trajectory is used in this animation.

Suggested slide statement: "A single trained configuration produced alternating
air motion with and without afferent feedback and passed assisted standing.
Transfer to alternating load-bearing ground steps remains unresolved."

For the dissertation session, these fit a tools/curriculum demonstration or
an appendix diagnostic table. Do not replace published or frozen results or
introduce a walking-success claim. Raw outputs are ignored local artifacts.

## Contact correction in progress

The legacy runner copied whole-foot normal force into both heel and toe
channels. The experimental evaluator now offers separate hindfoot/forefoot
regions based on actual contact locations, plus optional transient hindfoot
loading-event excitation of contralateral RG-F. This is a regional virtual
sensor partition, not calibrated heel/toe pads. Stage-4 comparisons and
force-accounting checks are recorded separately before any promotion.
