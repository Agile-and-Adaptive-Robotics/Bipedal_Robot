# Curation flags for Ben (2026-09-28 v2 campaign — hyperlinked edition)

164 flagged papers. Every paper title below is a clickable link
straight to its Airtable record (opens in your browser, logged in). DOI links
open the publisher page. Categories are auto-tagged from the flag wording —
skimming 'Data concerns' first is worthwhile. Regenerate after future batches:
`myo python build_flags_doc.py`. No Airtable API calls are used.

## Data concerns (possible duplicate / wrong metadata) (4)

- [Carlson-Kuhta 1998 — Forms of forward quadrupedal locomotion. II. A comparison of posture, hindlimb kinematics,](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rec2EuBpkiSUE9ZtV) · [doi](https://doi.org/10.1152/jn.1998.79.4.1687)  
  Abstract is identical to Smith 1998 'Forms III' (downslope) despite the title reading 'II ... upslope and level walking' - likely a metadata/abstract swap; replace with the true upslope abstract before relying on upslope-specific conclusions.
- [Molkov 2023 — Sensory Feedback and Central Neuronal Interactions in Mouse Locomotion](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recFqSV9WOXZsTiV9) · [doi](https://doi.org/10.1101/2023.10.31.564886)  
  probable duplicate of recCjTyA6Mz9tHw0G (published version)
- [Pearson and Kg 1993 — Common Principles of Motor Control in Vertebrates and Invertebrates](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recATBxX2a8F6fw7l) · [doi](https://doi.org/10.1146/annurev.ne.16.030193.001405)  
  ABSTRACT MISMATCH: the queue's abstract is default-mode-network text, not Pearson's motor-control review; all content fields left empty per the grounding rule. Re-fetch the abstract from DOI 10.1146/annurev.ne.16.030193.001405 and re-curate. is_review was inferred from the Annual Reviews title/venue only.
- [Pratt 1991 — Functionally complex muscles of the cat hindlimb. I. Patterns of activation across sartori](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recKW2KMXXHAHCpKi) · [doi](https://doi.org/10.1007/bf00229406)  
  Title says 'across sartorius' but the abstract describes biceps femoris compartments - title/abstract mismatch; verify which record is correct.

## New-pathway / vocabulary proposals (92)

- [Barbeau 1999 — Tapping into spinal circuits to restore motor function.](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rec3wynmDwk03ghmF) · [doi](https://doi.org/10.1016/s0165-0173(99)00008-9)  
  Chick embryo rhythmogenesis is discussed but 'Bird/Chick' is absent from the animals vocabulary - mapped to 'Vertebrates'; consider adding.
- [Beres-Jones and Harkema 2004 — The human spinal cord interprets velocity-dependent afferent input during stepping](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rec7PY2fLq3auXrSs) · [doi](https://doi.org/10.1093/brain/awh252)  
  The velocity-dependent afferent class is not named in the abstract; afferents left empty.
- [Bosco and Poppele 2001 — Proprioception from a spinocerebellar perspective.](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recRarTZSTWhAEBWn) · [doi](https://doi.org/10.1152/physrev.2001.81.2.539)  
  Abstract does not name the species (DSCT work is classically cat); animals left empty per the grounding rule.
- [Britz 2015 — A genetically defined asymmetry underlies the inhibitory control of flexor–extensor locomo](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recpqM3xOZzQDCuEm) · [doi](https://doi.org/10.7554/elife.04718)  
  Abstract never names the species (V1/V2b ablation in awake behaving animals); animals left empty rather than assumed.
- [Brooks 2015 — Learning to expect the unexpected: rapid updating in primate cerebellum during voluntary s](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recMGok8WlXgIUDE5) · [doi](https://doi.org/10.1038/nn.4077)  
  Primate species not specified in the abstract; animals recorded as the closest vocabulary term Mammals.
- [Bui 2013 — Circuits for Grasping: Spinal dI3 Interneurons Mediate Cutaneous Control of Motor Behavior](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/reclNWiYPtVRdFpb1) · [doi](https://doi.org/10.1016/j.neuron.2013.02.007)  
  Cutaneous -> dI3 -> motoneuron grip circuit fits no current Feedback vocabulary (grasp control rather than locomotor phase switching); feedback left empty.
- [Burke 1999 — The use of state-dependent modulation of spinal reflexes as a tool to investigate the orga](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rec5TrEW6v3HlB4bx) · [doi](https://doi.org/10.1007/s002210050847)  
  Abstract names no species (classic lesion/l-DOPA/FRA experiments implied), so animals left empty; the phase-gated oligosynaptic reflex pathways described have no one-to-one feedback vocabulary entries.
- [Burke 2001 — Patterns of locomotor drive to motoneurons and last-order interneurons: clues to the struc](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recXedSz4OBEuh3EX) · [doi](https://doi.org/10.1152/jn.2001.86.1.447)  
  Cutaneous reflex pathway modulation to FDL is described, but the abstract does not specify a stance-modification or flexor-excitation form — no exact vocabulary fit, so feedback left empty.
- [Buschmann 2015 — Controlling legs for locomotion—insights from robotics and neurobiology](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recMHdyJNF0qGg4Wp) · [doi](https://doi.org/10.1088/1748-3190/10/4/041001)  
  Abstract names no specific animal preparations (only common principles across species); animals left empty.
- [Chacon 2023 — Lumbar V3 interneurons provide direct excitatory synaptic input onto thoracic sympathetic ](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recEBU9QYruCYKihU) · [doi](https://doi.org/10.3389/fncir.2023.1235181)  
  Abstract does not state the species (optogenetic spinal cord and slice preparations) — animals left empty pending full text. The V3-to-SPN excitatory propriospinal pathway is not an afferent feedback pathway and has no vocabulary entry.
- [Chopek 2018 — Sub-populations of Spinal V3 Interneurons Form Focal Modules of Layered Pre-motor Microcir](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recQA3AJ46DeMxTIb) · [doi](https://doi.org/10.1016/j.celrep.2018.08.095)  
  Motoneuron-to-V3 recurrent glutamatergic excitation is a recurrent-excitation motif absent from the Feedback vocabulary (existing recurrent/Renshaw entries are inhibitory).
- [Council 2014 — Deadbeat control with (almost) no sensing in a hybrid model of legged locomotion](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rectZDCJfjdHtvPJI) · [doi](https://doi.org/10.1109/icamechs.2014.6911592)  
  Geometric/feedforward stabilization without sensing — conceptually adjacent to 'Biomechanically mediated preflexive feedback' but explicitly not feedback; feedback vocabulary left empty.
- [Dewolf 2022 — Left-right locomotor coordination in human neonates](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recRAV1TB6KDUwAPV) · [doi](https://doi.org/10.1523/jneurosci.0612-22.2022)  
  Neonatal load (midstance extensor activation) and hip-position (release-triggered swing) feedback rules are described without afferent types; candidate vocabulary additions for Ben.
- [Dietz 1994 — HUMAN NEURONAL INTERLIMB COORDINATION DURING SPLIT-BELT LOCOMOTION](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recOicfoYh5TcIN93) · [doi](https://doi.org/10.1007/bf00227344)  
  Abstract attributes the adaptation to 'proprioceptive feedback' without specifying afferent types or a specific pathway, so afferents and feedback were left empty.
- [Donelan and Pearson 2004 — Contribution of sensory feedback to ongoing ankle extensor activity during the stance phas](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recUDhhSVfIJ5cJg9) · [doi](https://doi.org/10.1139/y04-043)  
  Abstract attributes the extensor enhancement to 'likely candidates' (II in humans, Ib in humans and cats) — 'type II excitatory' and 'Ib excitatory' are tagged on that candidate basis, not a demonstrated receptor attribution.
- [Duysens 1990 — Gating and reversal of reflexes in ankle muscles during human walking](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recRsTA1RaPWpafjV) · [doi](https://doi.org/10.1007/bf00231254)  
  Ankle nerve (tibial/sural) stimulation at perception threshold; the abstract does not specify the afferent class (likely cutaneous/FRA-class), so 'Cutaneous stance modification' was not asserted and feedback was left empty.
- [Duysens 2002 — A walking robot called human: lessons to be learned from neural control of locomotion.](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recQTXZrLF3QI8qUo) · [doi](https://doi.org/10.1016/s0021-9290(01)00187-7)  
  Load-triggered extensor activation and unload-gated flexor initiation are described functionally without afferent types (Ib-like in the cat literature) and fit no exact vocabulary entry; feedback left empty rather than asserting Ib.
- [Duysens 2006 — How deletions in a model could help explain deletions in the laboratory.](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/reci17p4EWu9R0BYg) · [doi](https://doi.org/10.1152/jn.00888.2005)  
  Abstract missing - substantive fields left empty per grounding rule; appears to be a brief Duysens commentary (likely on modeling burst deletions); Ben should decide whether to fetch the PDF and how to classify it.
- [Ekeberg 1998 — A neuro-mechanical model of legged locomotion: single leg control](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/receswW0WAOzeYQcm) · [doi](https://doi.org/10.1007/s004220050468)  
  Sensory entrainment of the phase generator and phase-gated fast feedback have no exact feedback vocabulary entry.
- [Engberg and Lundberg 1969 — An electromyographic analysis of muscular activity in the hindlimb of the cat during unres](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rechFy16bSaCN7BJR) · [doi](https://doi.org/10.1111/j.1748-1716.1969.tb04415.x)  
  Ia reflex regulation is discussed as a possibility, not demonstrated; feedback left empty.
- [Frigon 2010 — Effects of ankle and hip muscle afferent inputs on rhythm generation during fictive locomo](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rec7NCgPDniURHbDa) · [doi](https://doi.org/10.1152/jn.01028.2009)  
  Nerve stimuli given at group I/group II strength without Ia-vs-Ib dissociation (afferents tagged Ia, Ib, and II to reflect the fiber classes engaged); demonstrated transitions - plantaris group I prolongs extension and terminates flexion, sartorius group II shortens flexion - have no exact feedback vocabulary entry, consider adding 'group I swing to stance'.
- [Frigon and Rossignol 2006 — Experiments and models of sensorimotor interactions during locomotion](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recYNXjJEH4sr5u1v) · [doi](https://doi.org/10.1007/s00422-006-0129-x)  
  Abstract names no species (likely cat/human content); animals left empty per grounding rule. The generic state-/phase-dependent gating it describes has no vocabulary entry, so feedback left empty.
- [Fuchs 2012 — Proprioceptive feedback reinforces centrally generated stepping patterns in the cockroach](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recue3mOOqQApmwgP) · [doi](https://doi.org/10.1242/jeb.067488)  
  Inter-leg proprioceptive reinforcement of the central pattern in cockroach — sensory organs and afferent types not named in the abstract; no exact feedback vocabulary fit.
- [Gao 2001 — Whisker Deafferentation and Rodent Whisking Patterns: Behavioral Evidence for a Central Pa](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recq863es8E4t9Vjx) · [doi](https://doi.org/10.1523/jneurosci.21-14-05374.2001)  
  Whisking rhythm survives bilateral deafferentation — conceptually 'Fictive locomotion without sensory feedback' but for whisking rather than locomotion; feedback left empty pending Ben's call.
- [Gossard 1989 — Intra-axonal recordings of cutaneous primary afferents during fictive locomotion in the ca](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recPs5aq2lbEr1ozm) · [doi](https://doi.org/10.1152/jn.1989.62.5.1177)  
  CPG-driven presynaptic (PAD) control of cutaneous afferents fits no current vocabulary entry; consider 'Cutaneous presynaptic inhibition'.
- [Gossard 1990 — Phase-dependent modulation of primary afferent depolarization in single cutaneous primary ](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recGNWivbxaXIUKPy) · [doi](https://doi.org/10.1016/0006-8993(90)90334-8)  
  Abstract demonstrates CPG-driven phase modulation of presynaptic inhibition of cutaneous afferents (PAD); the feedback vocabulary has no cutaneous presynaptic inhibition term, so feedback left empty.
- [Gregor 2006 — Mechanics of slope walking in the cat: quantification of muscle load, length change, and a](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recHGgTieMzq4gPKA) · [doi](https://doi.org/10.1152/jn.01300.2004)  
  Abstract names feedback only as 'muscle length and force, and paw pad cutaneous afferents' - spindle and Golgi tendon classes are implied but unnamed, so only Cutaneous tagged; no pathway is demonstrated, so feedback left empty.
- [Grillner 2011 — Control of Locomotion in Bipeds, Tetrapods, and Fish](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recLKbcEm4la8JtDV) · [doi](https://doi.org/10.1002/cphy.cp010226)  
  Abstract is a section table of contents only; the feedback-related items (load sensitivity, hip position, group II/III gating) appear as topic headings without demonstrated pathways, so feedback was left empty rather than guessed.
- [Grillner and Rossignol 1978 — On the initiation of the swing phase of locomotion in chronic spinal cats.](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recyqCjQ3Ha4n9k1r) · [doi](https://doi.org/10.1016/0006-8993(78)90973-3)  
  Abstract names no afferent type for the hip-position signal (the vocabulary item 'Ia or II stance to swing' would require that assumption); the contralateral-phase gating also has no vocabulary entry.
- [Harischandra 2011 — Sensory feedback plays a significant role in generating walking gait and in gait transitio](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recwgCoGUA2e5cPSH) · [doi](https://doi.org/10.3389/fnbot.2011.00003)  
  Feedback pathway (late-stance stretch-receptor input aiding gait transition and essential for walking coordination) has no vocabulary entry - stretch receptors are not typed as Ia/II in the abstract.
- [Harkema 1997 — Human Lumbosacral Spinal Cord Interprets Loading During Stepping](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recewflUcft50Rh2J) · [doi](https://doi.org/10.1152/jn.1997.77.2.797)  
  Load-dependent EMG modulation reported with afferent type unspecified (presumably Ib or plantar load receptors) and no vocabulary item for load-dependent extensor enhancement; feedback left empty.
- [Harnie 2024 — Forelimb movements contribute to hindlimb cutaneous reflexes during locomotion in cats](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recp5dgFg2CMMDQAz) · [doi](https://doi.org/10.1152/jn.00104.2024)  
  Forelimb-movement gating of hindlimb cutaneous long-latency reflexes has no vocabulary entry; feedback left empty.
- [Hasegawa 2012 — Pseudo-proprioceptive motion feedback by electric stimulation](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recWr7gpW0WgDVeff) · [doi](https://doi.org/10.1109/mhs.2012.6492480)  
  Electro-tactile sensory-substitution feedback; no natural afferent pathway is demonstrated — afferents and feedback left empty. No vocabulary entry fits prosthetic substitution channels.
- [Hiebert and Pearson 1999 — Contribution of Sensory Feedback to the Generation of Extensor Activity During Walking in ](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recpozMUyMMicDWNu) · [doi](https://doi.org/10.1152/jn.1999.81.2.758)  
  Extensor burst reinforcement is attributed to unnamed 'proprioceptors' (loading/unloading and dorsal-root manipulations — classic Ib autogenic excitation territory), but no afferent type is named; afferents and feedback left empty.
- [Hinckley 2010 — Sensory Modulation of Locomotor-Like Membrane Oscillations in Hb9-Expressing Interneurons](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rec5I9UPZcLsoKZrS) · [doi](https://doi.org/10.1152/jn.00996.2009)  
  Abstract specifies only 'low-threshold, presumably muscle afferents' (no Ia/II dissociation), so afferents left empty; the demonstrated monosynaptic excitation of Hb9 interneurons and phase-dependent cycle resetting fit no exact feedback vocabulary entry.
- [Hofstoetter 2017 — Probing the Human Spinal Locomotor Circuits by Phasic Step-Induced Feedback and by Tonic E](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recuGq3xm0qslYtBg) · [doi](https://doi.org/10.2174/1381612822666161214144655)  
  No feedback vocabulary entry for 'entrainment of spinal reflex circuits by cyclic proprioceptive feedback' — candidate new term for Ben.
- [Hubbuch 2015 — Proprioceptive feedback contributes to the adaptation toward an economical gait pattern](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recuI1q2b43ma7DG2) · [doi](https://doi.org/10.1016/j.jbiomech.2015.04.024)  
  Afferent type behind the vibration-disrupted 'proprioception' is not specified in the abstract (tendon vibration typically implicates Ia) — afferents left empty.
- [Hultborn 2006 — Spinal reflexes, mechanisms and concepts: from Eccles to Lundberg and beyond.](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recyWbN0VMRB86b0W) · [doi](https://doi.org/10.1016/j.pneurobio.2006.04.001)  
  Renshaw recurrent inhibition and (untyped) presynaptic inhibition are core topics but absent from the feedback vocabulary; afferents inferred as Ia/II/Ib from 'muscle spindles and Golgi tendon organs' in the abstract.
- [Iwasaki 2006 — Sensory Feedback Mechanism Underlying Entrainment of Central Pattern Generator to Mechanic](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/reciacLFlfKR61mAv) · [doi](https://doi.org/10.1007/s00422-005-0047-3)  
  Demonstrates sensory-feedback entrainment of the CPG to mechanical resonance (positive coupling plus complete adaptation); no matching feedback vocabulary entry - candidate 'CPG resonance entrainment' term.
- [Jankowska and McCrea 1983 — Shared reflex pathways from Ib tendon organ afferents and Ia muscle spindle afferents in t](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/reciJ8qFe7RzrrVyI) · [doi](https://doi.org/10.1113/jphysiol.1983.sp014663)  
  Core finding is Ia-Ib convergence on shared premotoneuronal interneurons (common feedback system) - no vocabulary entry covers shared-pathway convergence; feedback left empty for Ben to decide (e.g. 'Ia-Ib shared interneuronal convergence').
- [Kasumacic 2012 — Vestibular-mediated synaptic inputs and pathways to sympathetic preganglionic neurons in t](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recM5HgQ3RhYKMbTq) · [doi](https://doi.org/10.1113/jphysiol.2012.234609)  
  Vestibular (VIIIth nerve) afferents and the vestibulosympathetic reflex have no entries in the afferent or Feedback vocabularies; both fields were left empty rather than forced.
- [Kim 2022 — Contribution of Afferent Feedback to Adaptive Hindlimb Walking in Cats: A Neuromusculoskel](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recbsjTuemJTPPkv5) · [doi](https://doi.org/10.3389/fbioe.2022.825149)  
  Abstract names no specific afferent types or reflex pathways (the feedback detail lives in the full text); afferents and feedback left empty.
- [Klarner 2010 — Contribution of load and length related manipulations to muscle responses during force per](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/reclGnyvIRax1xURR) · [doi](https://doi.org/10.14288/1.0071394)  
  Afferents left empty: the abstract attributes responses to 'load sensitive and length sensitive afferents' without naming Ia/Ib/II.
- [Krause 2013 — Central drive and proprioceptive control of antennal movements in the walking stick insect](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/reca494vWEYiEgpAM) · [doi](https://doi.org/10.1016/j.jphysparis.2012.06.001)  
  Antennal hair fields regulate joint working ranges, but the feedback vocabulary only contains trochanteral hair-plate items; left empty.
- [Kriellaars 1994 — Mechanical entrainment of fictive locomotion in the decerebrate cat.](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recHnYpaDZj46nxzi) · [doi](https://doi.org/10.1152/jn.1994.71.6.2074)  
  Afferents are described only as 'low-threshold, stretch-sensitive muscle afferents' (spindles) with no Ia/II split, so afferents left empty; feedback tagged 'Ia or II stance to swing' for extensor-stretch prolongation of extension.
- [Laflamme 2023 — Distinct roles of spinal commissural interneurons in transmission of contralateral sensory](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recIjP2Mi5teZyaov) · [doi](https://doi.org/10.1016/j.cub.2023.07.014)  
  Crossed-reflex roles of V0/V3 CINs are described without afferent types; no matching feedback vocabulary term, so feedback left empty.
- [Leonardis 2014 — Multisensory Feedback Can Enhance Embodiment Within an Enriched Virtual Walking Scenario](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recgRFsKPiTtwLYEQ) · [doi](https://doi.org/10.1162/pres_a_00190)  
  Tangential to locomotor sensory-feedback control (VR embodiment); vestibular feedback is not in the afferent vocabulary - candidate for exclusion.
- [Li 2019 — Flexor and Extensor Ankle Afferents Broadly Innervate Locomotor Spinal Shox2 Neurons and I](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recekls2eVExj4oy8) · [doi](https://doi.org/10.3389/fncel.2019.00452)  
  Abstract demonstrates low-threshold ankle-afferrent phase resetting of fictive locomotion to flexion but never types the afferents (Ia/II/Ib), so afferents and feedback are left empty; no vocabulary item covers untyped proprioceptive reset.
- [Mari 2024 — Changes in intra- and interlimb reflexes from forelimb cutaneous afferents after staggered](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recB7n3thv2EEgswz) · [doi](https://doi.org/10.1113/jp286808)  
  Forelimb-to-hindlimb cutaneous interlimb reflex reorganization (homolateral/diagonal responses after incomplete SCI) fits no existing feedback vocabulary entry, so feedback is left empty; propose a 'Cutaneous interlimb reflex' term if Ben wants this family captured.
- [Markin 2010 — Afferent control of locomotor CPG: insights from a simple neuromechanical model](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recnEhLZMGIbVnYto) · [doi](https://doi.org/10.1111/j.1749-6632.2010.05435.x)  
  Model paper with ensemble Ia/II/Ib feedback modulating CPG phase durations; no vocabulary entry covers generic CPG entrainment by afferent feedback, so feedback left empty.
- [Maxwell and Soteropoulos 2020 — The mammalian spinal commissural system: properties and functions.](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recbNW3rCTH0YS9k7) · [doi](https://doi.org/10.1152/jn.00347.2019)  
  Commissural (left-right coordinating) pathways are not represented in the feedback vocabulary; left empty despite the paper's direct relevance to interleg coordination.
- [McCrea 1995 — Disynaptic group I excitation of synergist ankle extensor motoneurones during fictive loco](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recBMsJS4yKuHNuo5) · [doi](https://doi.org/10.1113/jphysiol.1995.sp020897)  
  Group Ia disynaptic excitation of synergists during extension (tendon-stretch-evoked EPSPs) has no vocabulary entry; 'Ib disynaptic excitation' is used for the Ib component, and an 'Ia disynaptic excitation' term is proposed if Ben wants the Ia component captured.
- [Mendes 2013 — Quantification of gait parameters in freely walking wild type and sensory deprived Drosoph](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rechzXQZNS6xn5IoW) · [doi](https://doi.org/10.7554/elife.00231)  
  Leg sensory neurons are not typed in the abstract (no chordotonal or campaniform named), so afferents and feedback left empty; the demonstrated pathway (proprioceptive feedback for step precision, coordination-independent) has no matching vocabulary entry.
- [Merlet 2020 — Cutaneous inputs from perineal region facilitates and modulates spinal locomotor activity ](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recZxoFgqdQh2aUJO) · [doi](https://doi.org/10.1101/2020.07.29.226530)  
  No vocabulary item for cutaneous excitation of the CPG rhythm with reflex-gain suppression; feedback left empty rather than forcing the existing cutaneous items.
- [Merlet 2021 — Cutaneous inputs from perineal region facilitate spinal locomotor activity and modulate cu](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recIpNrSMGJDUWbbu) · [doi](https://doi.org/10.1002/jnr.24791)  
  Perineally driven, state-dependent modulation of cutaneous reflex transmission has no exact feedback vocabulary term; feedback left empty.
- [Molkov 2024 — Sensory feedback and central neuronal interactions in mouse locomotion](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recCjTyA6Mz9tHw0G) · [doi](https://doi.org/10.1098/rsos.240207)  
  Postural-imbalance feedback controlling swing-to-stance transitions and walking direction has no Feedback vocabulary entry — consider proposing one. Near-identical bioRxiv preprint (2023) appears in this batch as recFqSV9WOXZsTiV9 — dedupe candidate.
- [Musselman and Yang 2007 — Loading the Limb During Rhythmic Leg Movements Lengthens the Duration of Both Flexion and ](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recAggUQj90C8YTU8) · [doi](https://doi.org/10.1152/jn.00891.2006)  
  The load-related feedback that prolongs the loaded phase is not typed to an afferent class in the abstract (an Ib-family pathway by mechanism) so feedback is left empty per the grounding rule.
- [Noga 1987 — The role of Renshaw cells in locomotion: antagonism of their excitation from motor axon co](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recI3QJS6rddKdgQH) · [doi](https://doi.org/10.1007/bf00236206)  
  Renshaw recurrent inhibition (rate control, not part of the CPG) has no feedback vocabulary term; feedback left empty.
- [Norris 2011 — Constancy and Variability in the Output of a Central Pattern Generator](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recd1FHQNmrOPeUNp) · [doi](https://doi.org/10.1523/jneurosci.5072-10.2011)  
  Leech (annelid) preparation not in the anim al vocabulary; animals left empty.
- [Ogihara and Yamazaki 2001 — Generation of human bipedal locomotion by a bio-mimetic neuro-musculo-skeletal model.](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recmkSPf9bTo465ER) · [doi](https://doi.org/10.1007/pl00007977)  
  Feedback left empty: the abstract credits 'reciprocal innervation in the muscle spindles' without typing it as Ia reciprocal inhibition.
- [Onushko 2013 — Hip proprioceptors preferentially modulate reflexes of the leg in human spinal cord injury](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recqN5btwC8zeVr9w) · [doi](https://doi.org/10.1152/jn.00261.2012)  
  Hip-position-dependent gating of knee/ankle stretch reflexes has no vocabulary entry; afferent types not named; feedback left empty.
- [Owaki 2013 — Simple robot suggests physical interlimb communication is essential for quadruped walking](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rec35pluvCkgMozzQ) · [doi](https://doi.org/10.1098/rsif.2012.0669)  
  Robot's sole coordinating signal is local per-leg force feedback (mechanosensory analog) - no exact feedback vocabulary term; candidate 'local leg force feedback'.
- [Ozeri-Engelhard 2022 — Inhibitory interneurons within the deep dorsal horn integrate convergent sensory input to ](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rec8SaFrCjyBMQ3J8) · [doi](https://doi.org/10.1101/2022.05.21.492933)  
  Abstract never names the species (intersectional genetics and parvalbumin labeling suggest mouse) so animals is left empty per the grounding rule; the dPV pathway (convergent cutaneous/proprioceptive inhibitory gating via dorsal horn interneurons) may also deserve its own feedback vocabulary entry beyond 'Cutaneous stance modification'.
- [Pang and Yang 2000 — The initiation of the swing phase in human infant stepping: importance of hip position and](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rec9RDRjijNiY2b26) · [doi](https://doi.org/10.1111/j.1469-7793.2000.00389.x)  
  Abstract does not type the afferents behind the hip-position and load regulation of swing initiation (in the vocabulary this would map to 'Ia or II stance to swing' plus a load/Ib-type gate) so feedback is left empty per the grounding rule.
- [Pazzaglia 2025 — Balancing central control and sensory feedback produces adaptable and robust locomotor pat](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recMZoiH13ytYevS7) · [doi](https://doi.org/10.1371/journal.pcbi.1012101)  
  Axial stretch (proprioceptive) feedback fits no Feedback vocabulary entry - feedback left empty per the guide; afferents mapped to Mechanosensory as the closest term. Consider adding a stretch-receptor feedback vocabulary entry.
- [Perreault 1999 — Proprioceptive Control of Extensor Activity during Fictive Scratching and Weight Support C](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recSRi5Pinh7wq0vV) · [doi](https://doi.org/10.1523/jneurosci.19-24-10966.1999)  
  Abstract never names the species (fictive locomotion/scratching/weight-support preparation, almost certainly cat); animals left empty per the grounding rule. 'Extensor group I' decomposed as Ia+Ib; the excitation is called oligosynaptic, so it was not mapped to 'Ib disynaptic excitation' or 'Ib excitatory'.
- [Perret and Cabelguen 1980 — Main characteristics of the hindlimb locomotor cycle in the decorticate cat with special r](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rec6Fd1OrPbNA2Yl0) · [doi](https://doi.org/10.1016/0006-8993(80)90207-3)  
  Alpha-gamma coactivation and FRA-pattern gating of bifunctional motoneurons have no exact feedback vocabulary entries; feedback left empty.
- [Procházka 1976 — Discharges of single hindlimb afferents in the freely moving cat.](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rec4XpUFoz8LpQ10x) · [doi](https://doi.org/10.1152/jn.1976.39.5.1090)  
  Brisk stretch of ankle extensors evoked a rapid (disynaptic or trisynaptic) reflex arc in the conscious cat; no vocabulary entry covers a disynaptic/trisynaptic Ia stretch-reflex excitation, so feedback left empty - consider adding one.
- [Quevedo 2000 — Group I disynaptic excitation of cat hindlimb flexor and bifunctional motoneurones during ](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rec3BRkL2st5JzzeI) · [doi](https://doi.org/10.1111/j.1469-7793.2000.t01-1-00549.x)  
  Also shows Ia spindle afferents evoking DISYNAPTIC excitation of flexor motoneurones - 'Ia monosynaptic excitation' does not fit; candidate vocabulary term 'Ia disynaptic excitation (flexors)'.
- [Refy 2023 — Dynamic spinal reflex adaptation during locomotor adaptation](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recZJU3Tv1NeOUDiG) · [doi](https://doi.org/10.1152/jn.00248.2023)  
  Abstract does not name the reflex modality or afferent type, so afferents and feedback left empty.
- [Rossignol 2006 — Dynamic Sensorimotor Interactions in Locomotion](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recTUr9yGNL2hrV5y) · [doi](https://doi.org/10.1152/physrev.00028.2005)  
  Review abstract names no specific species (animals left empty per grounding rule) and no muscle afferent types; the extensor-proprioceptive stance timing/amplitude modulation has no exact vocabulary entry — 'Cutaneous stance modification' tagged for the foot-placement role only.
- [Rybak 2024 — Operation regimes of spinal circuits controlling locomotion and the role of supraspinal dr](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rechmhCcKA0ILQWKG) · [doi](https://doi.org/10.7554/elife.98841)  
  Slow-speed phase transitions depend on 'sensory feedback related to limb extension' with no afferent group named (Ia/II length vs Ib force), so feedback left empty - Ben may want a 'limb-extension feedback phase transition' vocabulary entry.
- [Sastra 2008 — A Biologically-Inspired Dynamic Legged Locomotion With a Modular Reconfigurable Robot](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recWPp5EyxMBGayi3) · [doi](https://doi.org/10.1115/dscc2008-2402)  
  Passive leg compliance is the stability mechanism (arguably 'Biomechanically mediated preflexive feedback'), but the abstract describes no explicit feedback pathway — feedback left empty. Pure robotics paper; no animal content.
- [Schomburg 1977 — Phase-dependent transmission in the excitatory propriospinal reflex pathway from forelimb ](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recOx68YKAQvYBUZD) · [doi](https://doi.org/10.1016/0304-3940(77)90187-2)  
  Phase-gated descending propriospinal excitation from forelimb afferents to lumbar motoneurons fits no current Feedback vocabulary entry; consider 'propriospinal interlimb excitation, phase-gated'. Afferent types not specified in the abstract.
- [Schomburg 1998 — Flexor reflex afferents reset the step cycle during fictive locomotion in the cat](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rec0Y7IJ9rNh4reiU) · [doi](https://doi.org/10.1007/s002210050522)  
  FRA-train reset of the rhythm (flexion-pattern reset interrupting extension) fits no exact feedback vocabulary term; candidate 'FRA resets step cycle'. Abstract also names joint afferents (no vocabulary entry) and extension-prolonging group I effects (Ia vs Ib unspecified).
- [Severini and Muñoz 2025 — A physiologically inspired hybrid CPG/Reflex controller for cycling simulations that gener](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recp7yrg2a8E7DgyR) · [doi](https://doi.org/10.1371/journal.pcbi.1013494)  
  The abstract does not enumerate the reflex pathways in the feedback component; afferents and feedback left empty.
- [Shaw 2015 — The significance of dynamical architecture for adaptive responses to mechanical loads duri](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recwRVXeyikVNCPwN) · [doi](https://doi.org/10.1007/s10827-014-0519-3)  
  Animal is Aplysia (mollusc), absent from the animals vocabulary; the load/proprioceptive-modulation pathway has no feedback vocabulary entry.
- [Shefchyk 1990 — Activity of interneurons within the L4 spinal segment of the cat during brainstem-evoked f](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recsA05QaEgbSWsbo) · [doi](https://doi.org/10.1007/bf00228156)  
  Group II afferent-to-midlumbar-interneuron-to-motor-nucleus pathway demonstrated, but the sign of the action on motoneurons is not stated in the abstract — no exact feedback vocabulary fit; feedback left empty.
- [Shinohara 2025 — Mechanisms of adaptive interlimb coordination to sudden ground loss: a neuromusculoskeleta](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recpHeyTiV5C1veLV) · [doi](https://doi.org/10.1101/2025.11.11.687930)  
  Afferent feedback controls fast-slow transition dynamics in the model, but the abstract names no afferent classes; feedback left empty.
- [Sinkjær 2000 — Major role for sensory feedback in soleus EMG activity in the stance phase of walking in m](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/reca4OL8XXDjk9p7L) · [doi](https://doi.org/10.1111/j.1469-7793.2000.00817.x)  
  Abstract implicates group II and/or Ib but names no synaptic pathway, so feedback left empty rather than guessing among the Ib/II vocabulary items.
- [Sra 2019 — Adding Proprioceptive Feedback to Virtual Reality Experiences Using Galvanic Vestibular St](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recawoXfKMSJgoPyx) · [doi](https://doi.org/10.1145/3290605.3300905)  
  Vestibular afferents and vestibular reflex pathways are not represented in the current afferent or feedback vocabularies; fields left empty.
- [Syed 1990 — In Vitro Reconstruction of the Respiratory Central Pattern Generator of the Mollusk <i>Lym](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recATJ2aCZl7aUXZY) · [doi](https://doi.org/10.1126/science.2218532)  
  The animal (pond snail Lymnaea stagnalis) has no exact name in the animals vocabulary, so animals is left empty; propose adding 'Mollusk' or 'Snail' if non-arthropod invertebrate CPGs should be captured.
- [Szczecinski 2021 — A computational model of insect campaniform sensilla predicts encoding of forces during wa](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/reczLE9Hbbjz1dwyz) · [doi](https://doi.org/10.1088/1748-3190/ac1ced)  
  The feedback vocabulary has only the trochanteral CS-to-MN item; this paper is a force-encoding transducer model without a specified central pathway, so feedback left empty.
- [Szczecinski and Quinn 2018 — Leg-local neural mechanisms for searching and learning enhance robotic locomotion](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rechOTg6oKWMbHpYM) · [doi](https://doi.org/10.1007/s00422-017-0726-x)  
  Abstract describes negative feedback of leg depressor force (load signals adjusting muscle activation) and a no-load search mode but never names the receptor (e.g. campaniform sensilla), so feedback left empty - candidate vocabulary entry 'leg depressor force negative feedback / load-gated searching'.
- [Tripodi 2011 — Motor antagonism exposed by spatial segregation and timing of neurogenesis](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recnBgLt46yZT88vr) · [doi](https://doi.org/10.1038/nature10538)  
  Abstract attributes extensor premotor targeting to unnamed 'proprioceptive' feedback (no afferent class); feedback left empty — developmental targeting, no vocabulary fit.
- [Von Stetina et al. 2005 — The Motor Circuit](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recrZGsaVKPx35Xdx) · [doi](https://doi.org/10.1016/S0074-7742(05)69005-8)  
  Animal preparation is C. elegans (nematode) — no matching name in the allowed animals vocabulary; animals left empty.
- [Webster-Wood 2020 — Control for multifunctionality: bioinspired control based on feeding in Aplysia californic](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rec4qaxGiMJ2k6c7e) · [doi](https://doi.org/10.1007/s00422-020-00851-9)  
  Aplysia californica (mollusk) is not in the animals vocabulary, so animals left empty - consider adding a Mollusk/Aplysia entry.
- [Whelan 1995 — Stimulation of the group I extensor afferents prolongs the stance phase in walking cats.](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recih7yOrlandW0mo) · [doi](https://doi.org/10.1007/bf00241961)  
  Effect is of mixed group I extensor input (authors could not separate Ia from Ib, though they consider Ib likely) and includes a swing-terminating action - no 'group I stance to swing' vocabulary entry, so feedback left empty for Ben to decide.
- [Wilmink and Nichols 2003 — Distribution of Heterogenic Reflexes Among the Quadriceps and Triceps Surae Muscles of the](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recjfbKv36iU9scmz) · [doi](https://doi.org/10.1152/jn.00833.2002)  
  Length and force feedback are not attributed to named afferents in the abstract (presumed spindle vs tendon organ) and 'interjoint force inhibition' has no exact vocabulary entry - feedback left empty for Ben to decide.
- [Yakovenko 2004 — Contribution of stretch reflexes to locomotor control: a modeling study](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recDGA2Yt2ZU5pYhG) · [doi](https://doi.org/10.1007/s00422-003-0449-z)  
  Abstract does not name the species of the bipedal planar model (human-like slow level walking) — animals left empty. State-gated autogenic Ia/Ib reflex has no exact vocabulary entry beyond the preflexive one.
- [Yang 1991 — Contribution of peripheral afferents to the activation of the soleus muscle during walking](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rechPjLTzHpL5JYzS) · [doi](https://doi.org/10.1007/bf00227094)  
  Abstract attributes the soleus stretch response to unspecified peripheral afferents (no Ia/spindle named), so afferents and feedback left empty per grounding rule.
- [Zhang 2008 — V3 spinal neurons establish a robust and balanced locomotor rhythm during walking.](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recpmLF4VjQIRfxJ7) · [doi](https://doi.org/10.1016/j.neuron.2008.09.027)  
  Abstract does not name the species; animals left empty rather than assumed.

## Animal-vocabulary gaps (8)

- [Elson 1992 — Identified proprioceptive afferents and motor rhythm entrainment in the crayfish walking s](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recm4urpzERKRV42b) · [doi](https://doi.org/10.1152/jn.1992.67.3.530)  
  Crayfish mapped to 'Arthropods' (no Crayfish list entry); TCMRO S/T afferent entrainment and burst-resetting pathways fit no current Feedback vocabulary (arthropod promote/remote phases rather than vertebrate stance/swing), so feedback is empty.
  # create the feedback vocabulary. Animals: Crayfish and Arthopods.
- [Jindrich 2009 — Maneuvers during legged locomotion](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recD8qqPsM8KJCYJX) · [doi](https://doi.org/10.1063/1.3143031)  
  Ostriches are discussed as the comparison biped but are not in the animals vocabulary.
  # Animals: Ostriches (and humans?)
- [LaBella 1992 — Low-threshold, short-latency cutaneous reflexes during fictive locomotion in the "semi-chr](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recLuYXbI1HoROtsh) · [doi](https://doi.org/10.1007/bf00231657)  
  Feedback mapping: short-latency excitation in extensors (triceps surae) maximal during their active period was mapped to Cutaneous stance modification, and excitation in the flexor semitendinosus to Cutaneous flexor excitation.
  # Again, I am confused by "Cutaneous". Like contact or like an external stimulus?
- [Maas 2007 — The effects of self-reinnervation of cat medial and lateral gastrocnemius muscles on hindl](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rec2RZbu3cjgkAJvv) · [doi](https://doi.org/10.1007/s00221-007-0938-8)  
  'Length feedback' from MG/LG mapped to spindle Ia/II (reinnervation also silences Ib); Cutaneous included because the abstract names it explicitly.
- [Rome 1993 — How Fish Power Swimming](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rec2LGJaT2IS5CkEM) · [doi](https://doi.org/10.1126/science.8332898)  
  Abstract says only 'fish' (no species); animals mapped to the closest allowed bucket 'Vertebrates' - consider adding 'Fish' to the vocabulary.
  # Add fish
- [Rossi-Durand 1993 — Peripheral proprioceptive modulation in crayfish walking leg by serotonin](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/reckkLAj4fWC3npA8) · [doi](https://doi.org/10.1016/0006-8993(93)91131-b)  
  Crayfish (crustacean) mapped to 'Arthropods' - no Crayfish entry in the animal list.
  # add crayfish
- [Schmitt 2002 — Dynamics and stability of legged locomotion in the horizontal plane: a test case using ins](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recN3bMEzLPbHVamF) · [doi](https://doi.org/10.1007/s00422-001-0300-3)  
  Passive-dynamic stability from leg mechanics was mapped to Biomechanically mediated preflexive feedback; the abstract does not use the term preflex itself.
  # Is preflex here like stance to swing transition?
- [Ting and Chiel 2017 — Muscle, Biomechanics, and Implications for Neural Control](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recMLbRa9efnWHXEx) · [doi](https://doi.org/10.1002/9781118873397.ch12)  
  Chapter explicitly covers both vertebrates and invertebrates; the animals vocabulary has no general invertebrate entry, so only Vertebrates was recorded.
# add invertebrates

## Judgment calls (tagging conservatism etc.) (60)

- [Akay 2004 — Signals from load sensors underlie interjoint coordination during stepping movements of th](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recMqk9O3UkRUMNMz) · [doi](https://doi.org/10.1152/jn.01271.2003)  
  The vocabulary entry names MN magnitude adjustment, but this paper demonstrates timing/switching of protractor-retractor motoneuron pools by trochanteral campaniform sensilla signals - consider a timing variant of that entry.
- [Akay 2020 — Sensory Feedback Control of Locomotor Pattern Generation in Cats and Mice.](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recfjnlS1Vy7XmC5v) · [doi](https://doi.org/10.1016/j.neuroscience.2020.05.008)  
  Abstract too thin to extract specific afferent types or pathways.
- [Akazawa 1982 — Modulation of stretch reflexes during locomotion in the mesencephalic cat](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recetENlBYcGrP7zO) · [doi](https://doi.org/10.1113/jphysiol.1982.sp014319)  
  Core finding is phase-dependent gain modulation of the Ia stretch reflex; no vocabulary item for reflex gain modulation.
  # Well can we call it "stance to swing", "swing to stance", "stance", "swing", and then add either "Ia inhibition" or "Ia excitation"? In general, we should allow for this logic algebra: "afferent" + "inhibit/excite" + "contralateral/ipsilateral/agonist/antagonist" + "joint" + "phase". In this case though it is the phase feeding back to the afferent. Do you understand?
- [Andersson and Grillner 1981 — Peripheral control of the cat's step cycle I. Phase dependent effects of ramp-movements of](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recPGqPZU3LRiZkX3) · [doi](https://doi.org/10.1111/j.1748-1716.1981.tb06867.x)  
  Hip-ramp phase reset and swing-triggering effects fit no current vocabulary entry; consider 'hip position afferent stance to swing'. The abstract does not specify the afferent type mediating the effect.
  # Exactly, see note above. 
- [Angel 1996 — Group I extensor afferents evoke disynaptic EPSPs in cat hindlimb extensor motorneurones d](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recWccr9LDvwi2WzR) · [doi](https://doi.org/10.1113/jphysiol.1996.sp021538)  
  Selective group Ia activation also evoked the disynaptic excitation; the vocabulary has no 'Ia disynaptic excitation' entry, so the finding is tagged under 'Ib disynaptic excitation' (the group I mixed-nerve result).
  # I would like to know the history of the types of afferent discovery. Did they think about load or position? Otherwise we might need to leave it 'I disynaptic excitation'. Or, it could be 'Ia disynaptic excitation' and 'Ib disynaptic excitation'. I think both exist in the stance-phase.
- [Angel 2005 — Candidate interneurones mediating group I disynaptic EPSPs in extensor motoneurones during](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/reca10selgQlFrpHA) · [doi](https://doi.org/10.1113/jphysiol.2004.076034)  
  Abstract specifies 'group I' without dissociating Ia from Ib; afferents listed as both, and the 'Ib disynaptic excitation' vocabulary item used as the standard name for this stance-phase group I pathway.
# see above
- [Asif 2012 — On the Improvement of Multi-Legged Locomotion over Difficult Terrains Using a Balance Stab](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/reciUGAHPamTLce0k) · [doi](https://doi.org/10.5772/7789)  
  Pure control-engineering paper - no sensory-afferent content; Ben may want to decide whether it belongs in SADb.
  # delete if not done already
- [Beer 1998 — Biorobotic approaches to the study of motor systems](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recW2QvwSd3rFpEda) · [doi](https://doi.org/10.1016/s0959-4388(98)80121-9)  
  Abstract is a three-sentence editorial summary; no species, afferents, or pathways — minimal curation possible.
- [Blackburn 2004 — Sex comparison of extensibility, passive, and active stiffness of the knee flexors](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recxUQ8JuDHyBXIk0) · [doi](https://doi.org/10.1016/j.clinbiomech.2003.09.003)  
  Marginal for a sensory-afferent database: biomechanical stiffness comparison with no afferent or reflex pathway content in the abstract.
  # I think we deleted it already
- [Dean 2013 — Proprioceptive Feedback and Preferred Patterns of Human Movement](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recUWzWrumdQrLybs) · [doi](https://doi.org/10.1097/jes.0b013e3182724bb0)  
  Abstract is an 'In Brief' hypothesis statement with no data — too thin for afferent typing or feedback pathway tags.
  # delete
- [Dimitrijević 1998 — Evidence for a Spinal Central Pattern Generator in Humansa](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/reccXh6N1BuE4gy39) · [doi](https://doi.org/10.1111/j.1749-6632.1998.tb09062.x)  
  Stimulation site is 'posterior structures' of the lumbar cord — afferent or dorsal-root involvement implied but not named in the abstract; see Hachmann 2021 for the large-diameter afferent interpretation.
- [Edgley and Jankowska 1987 — An interneuronal relay for group I and II muscle afferents in the midlumbar segments of th](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recf2nnsCOlfEvGxV) · [doi](https://doi.org/10.1113/jphysiol.1987.sp016676)  
  Group I afferents decomposed to Ia+Ib in the afferents field (standard definition); group I disynaptic actions could not be assigned to a specific vocabulary item.
- [Ekeberg and Pearson 2005 — Computer simulation of stepping in the hind legs of the cat: an examination of mechanisms ](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recYC3Vr0fYcePhr4) · [doi](https://doi.org/10.1152/jn.00065.2005)  
  Abstract names load (ankle extensor force) and hip-angle signals without naming receptor types; 'Ib stance to swing' is tagged because the load-sensitive stance-termination pathway is the paper's core mechanism, but the receptor attribution is by convention.
  # sure 'Ib stance to swing' fits but we can do 'ankle Ib feedback to hip during stance' it sounds like
- [Gao 2019 — Dendritic Neuron Model With Effective Learning Algorithms for Classification, Approximatio](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recvwIzskqiUp9ikx) · [doi](https://doi.org/10.1109/tnnls.2018.2846646)  
  Off-topic for a sensory afferent database: pure machine-learning neuron-model paper with no afferent or locomotor content — Ben may want to reclassify or drop.
  # delete
- [Gerstner 2009 — How Good Are Neuron Models?](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recKCncNuVAa7rZ4g) · [doi](https://doi.org/10.1126/science.1181936)  
  Abstract is a two-sentence teaser with no methods or results; only a minimal note can be grounded, and commentary-versus-review classification is uncertain (is_review left false).
  # delete
- [Gomar 2014 — Digital Multiplierless Implementation of Biological Adaptive-Exponential Neuron Model](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recwAaYQRw2Psco5n) · [doi](https://doi.org/10.1109/tcsi.2013.2286030)  
  Off-topic for a sensory afferent database: digital hardware neuron implementation with no afferent or locomotor content — Ben may want to reclassify.
  # delete
- [Gossard 2011 — Chapter 2--the spinal generation of phases and cycle duration.](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recfIKtOOvoPaRCFc) · [doi](https://doi.org/10.1016/b978-0-444-53825-3.00007-3)  
  Ankle dorsiflexion prolongs extension during fictive locomotion, but the afferent type is not specified and no vocabulary item covers stance prolongation by extensor stretch.
  # maybe I'm just tired but you can say 'Ankle dorsiflexion prolongs extension during fictive locomotion,'
- [Gubina 1974 — On the Dynamic Stability of Biped Locomotion](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recZVayeFS2NnRjwz) · [doi](https://doi.org/10.1109/tbme.1974.324294)  
  Pure control-theory paper (state feedback laws, no afferent content); feedback vocabulary not applicable here.
  # delete
- [Guevremont 2007 — Physiologically based controller for generating overground locomotion using functional ele](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recDbjiLWM58hFrG9) · [doi](https://doi.org/10.1152/jn.01177.2006)  
  Sensory-driven flexion/extension phase transitions were implemented with external sensors (force plates, accelerometers) — no afferent type identified, so no Feedback vocabulary entry applies.
  # I'll make a note to look at this again. It sounds like you can get Ib, II, and maybe from Ia from what you describe.
- [Hachmann 2021 — Epidural spinal cord stimulation as an intervention for motor recovery after motor complet](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recZblEttLs7yl33F) · [doi](https://doi.org/10.1152/jn.00020.2021)  
  Afferents described only as 'large diameter dorsal root proprioceptive' with no exact match in the afferents vocabulary (Ia/Ib/II not dissociated); feedback set via the 'large diameter spinal afferent stimulation' item.
- [Haque 2019 — Characterization of WT1-expressing interneurons and investigation of their role in locomot](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recbdmPohBK1Qlvi7) · [doi](https://doi.org/10.7939/r3-jt04-kz83)  
  University thesis (repository DOI), not a journal article; neonatal mouse data, species taken from the abstract text.
- [Hart and Giszter 2010 — A Neural Basis for Motor Primitives in the Spinal Cord](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/reclFhqtYHPAnaJ7k) · [doi](https://doi.org/10.1523/jneurosci.5894-08.2010)  
  Afferents list Cutaneous only; the 28 'proprioceptive neurons' in the abstract are not typed (Ia/II/Ib), so no proprioceptive entry was assigned.
- [Hellekes 2012 — Control of reflex reversal in stick insect walking: effects of intersegmental signals, cha](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rech1JxuXg5zKkApI) · [doi](https://doi.org/10.1152/jn.00718.2011)  
  Core finding is the task- and intersegmentally gated reversal between chordotonal excitation and inhibition; no dedicated vocabulary item for reflex reversal itself.
- [Hodgkin 1952 — A quantitative description of membrane current and its application to conduction and excit](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recYcaqhK8ItCNuvf) · [doi](https://doi.org/10.1113/jphysiol.1952.sp004764)  
  Abstract field contains publisher boilerplate only (no scientific content); all fields grounded in title alone, and is_model inferred from the title. Suggest re-fetching the real abstract.
- [Hunt 2017 — Development and Training of a Neural Controller for Hind Leg Walking in a Dog Robot](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recM0Spre5cjj30My) · [doi](https://doi.org/10.3389/fnbot.2017.00018)  
  Robotics study with no animal experiments; animals entries record the modeled biology only (dog hind legs, mammalian locomotion models).
  ## dog, mammal animal categories. You should find the rules if you haven't already. This is my advisor's work. Get it right.

- [Hägglund 2013 — Optogenetic dissection reveals multiple rhythmogenic modules underlying locomotion](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recYd8bwl6gcq1KwR) · [doi](https://doi.org/10.1073/pnas.1304365110)  
  OpenAlex topic tag says Zebrafish but the abstract concerns the mammalian CPG; animal recorded as Mammals. Ben may want the species confirmed from the full text.
  ## zebrafish and Mammals. Did you get rules from this paper?

- [Izhikevich 2006 — Dynamical Systems in Neuroscience: The Geometry of Excitability and Bursting](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recTyf20KciVGaJVG) · [doi](https://doi.org/10.7551/mitpress/2526.001.0001)  
  Monograph/textbook, not primary research or a review essay; curated as methods background. No animal or afferent content in the abstract.
## delete all monographs and textbooks

- [Jankowska and Edgley 2010 — Functional subdivision of feline spinal interneurons in reflex pathways from group Ib and ](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recyaLy0ToIQpfrnw) · [doi](https://doi.org/10.1111/j.1460-9568.2010.07354.x)  
  Core claim - shared group I/II premotor integration rather than separate Ib and II populations - has no feedback vocabulary entry.
## we'll come back to it. but it sounds like "type I/II premotor integration"? The title has Ib, what other rules have you got from this paper?

- [Jessell 2000 — Neuronal specification in the spinal cord: inductive signals and transcriptional codes](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recmtKi7BqBJSH9my) · [doi](https://doi.org/10.1038/35049541)  
  Developmental review with no sensory-afferent content in the abstract; Ben may want to review its inclusion in the afferent database.
  ## delete

- [Karayannidou 2009 — Maintenance of Lateral Stability During Standing and Walking in the Cat](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recFd4pICpFKMQZ4h) · [doi](https://doi.org/10.1152/jn.90934.2008)  
  Lateral-perturbation corrective responses are a postural feedback pathway with no Feedback vocabulary entry.
  ## "lateral postural reflex"

- [Knikou 2010 — Neural control of locomotion and training-induced plasticity after spinal and cerebral les](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rechb1DliCuuKLIy7) · [doi](https://doi.org/10.1016/j.clinph.2010.01.039)  
  Abstract cites 'reduced animal preparations' for afferent regulation of the CPG without naming species, so animals restricted to Human.
  ## read the full text

- [Knuesel 2011 — Effects of muscle dynamics and proprioceptive feedback on the kinematics and CPG activity ](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recLvq2dSfUodmjTa) · [doi](https://doi.org/10.1186/1471-2202-12-s1-p158)  
  Abstract is two figure captions (conference abstract, about 50 words) - too thin to support afferent types, feedback vocabulary, or an animal_study field.
  ## new rule for SADb on airtable: If it's a conference abstract and there is no paper or poster accessible, the paper should probably just be deleted.

- [Koditschek 2004 — Mechanical aspects of legged locomotion control](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/reccLn6XYsnyAJmyX) · [doi](https://doi.org/10.1016/j.asd.2004.06.003)  
  'Biomechanically mediated preflexive feedback' inferred from the muscular-skeletal-neural integration framing; the abstract does not use the word preflex.
  ## Is the term in quotes one you made or was already there?

- [Kooij 2000 — An adaptive model of sensory integration in a dynamic environment applied to human stance ](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rec4qLdYPmyjgG02Z) · [doi](https://doi.org/10.1007/s004220000196)  
  Balance sensory-integration model (vestibular, visual, somatosensory weighting) lies outside the current afferent/feedback vocabulary; vestibular and visual inputs have no SADb tags.
  ## Great. sound like you can create at least 3 feedback tags to this article.

- [Laliberté 2019 — Propriospinal Neurons: Essential Elements of Locomotor Control in the Intact and Possibly ](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recBq2I1f3v1wg9i1) · [doi](https://doi.org/10.3389/fncel.2019.00512)  
  The abstract references unspecified 'animal models of SCI' without naming species, so only Human is listed; add the relevant species after checking the full text.
- [Mazzaro 2006 — Afferent-mediated modulation of the soleus muscle activity during the stance phase of huma](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recY2yVKiHmqH2TuE) · [doi](https://doi.org/10.1007/s00221-006-0451-5)  
  Ia contribution is stated as 'at least partially' mediating the increments with no pathway detail given, so no Ia feedback term tagged; group II result tagged 'type II excitatory'.
  ## ok. See if position, not just velocity, feedback is mentioned and use that as a stand in for Ia feedback.

- [Merlet 2021 — Inhibition and Facilitation of the Spinal Locomotor Central Pattern Generator and Reflex C](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recz45EOs5F1ttdpq) · [doi](https://doi.org/10.3389/fnins.2021.720542)  
  Lumbar-region (inhibitory) versus perineal (facilitatory) somatosensory modulation of the CPG and reflex gain has no feedback vocabulary entry.
- [Naris 2020 — A neuromechanical model exploring the role of the common inhibitor motor neuron in insect ](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recDBS8FWkbL8ogi1) · [doi](https://doi.org/10.1007/s00422-019-00811-y)  
  Core mechanism is efferent (common inhibitor MN plus slow-fiber dynamics), so no afferent Feedback vocabulary entry applies.
- [Nelson and Quinn 1998 — Posture control of a cockroach-like robot](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recifDtY9KNMQFxmg) · [doi](https://doi.org/10.1109/robot.1998.676348)  
  Robot-engineering paper - no afferent content; Ben may want to decide whether it belongs in SADb.
  ## Yeah. Quinn is my boss's boss, sorta. There's some Ia or insect equivalent feedback in there.

- [Orlovsky 1999 — Neuronal Control of LocomotionFrom Mollusc to Man](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recH0svuO344o1Stg) · [doi](https://doi.org/10.1093/acprof:oso/9780198524052.001.0001)  
  Monograph, so is_review kept false per the guide (reviews exclude monographs); abstract too general to ground animals beyond Human (leech is named but not in the animal list) or any mechanism field.
  ## delete then

- [Pearson and Rossignol 1991 — Fictive motor patterns in chronic spinal cats](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rec4xSUdAGn6omVtB) · [doi](https://doi.org/10.1152/jn.1991.66.6.1874)  
  Tagged 'Fictive locomotion without sensory feedback' for the paralyzed chronic-spinal fictive preparation; rhythms were evoked by cutaneous (perineal/paw) stimulation, and the leg-position afferent class modulating burst timing is not named.
I really don't like how you (and the papers, tbf), use cutaneous. Stimulus or TENS/EMS seems more fitting. To me, I also think that sometimes cutaneous means normal contact, shear contact, or something else. Don't we have "deafferented" already? So would this be two types of feedback tags, like "deafferented" and "Stimulus to PF layer".

- [Perreault 1999 — Depression of muscle and cutaneous afferent-evoked monosynaptic field potentials during fi](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recSh9tW4eoP1c5tq) · [doi](https://doi.org/10.1111/j.1469-7793.1999.00691.x)  
  Abstract specifies 'group I' without an Ia/Ib split; afferents tagged Ia+Ib by convention. The demonstrated depression is generalized locomotor-state gating of group I/II/cutaneous afferent transmission; closest vocabulary entry is 'Ia presynaptic inhibition' — consider a broader locomotor-state transmission-depression term.
- [Pratt 2012 — Capturability-based analysis and control of legged locomotion, Part 2: Application to M2V2](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recQrYCoH7uVZ4zvP) · [doi](https://doi.org/10.1177/0278364912452762)  
  Robotics control paper (capturability framework) with no neural or afferent content; kept as a balance-framework reference. Scope decision for Ben.
- [Pratt and Jordan 1987 — Ia inhibitory interneurons and Renshaw cells as contributors to the spinal mechanisms of f](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recBowk2bNdwFfR9o) · [doi](https://doi.org/10.1152/jn.1987.57.1.56)  
  The queue abstract is truncated mid-sentence at its end ('...in which the motoneuron was most'); curation deliberately excludes the cut-off clause.
- [Procházka and Prochazka 2011 — Proprioceptive Feedback and Movement Regulation](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recXastPRIyZzzaN1) · [doi](https://doi.org/10.1002/cphy.cp120103)  
  Feedback entries inferred from the abstract's table-of-contents section titles ('Proportional Control (Stretch Reflexes)' -> Ia monosynaptic excitation; 'Positive Force Feedback?' -> Ib excitatory). Invertebrate proprioceptors are covered but no specific invertebrate taxon is named in the abstract, so none tagged.
- [Quinn 2011 — Novel Locomotion via Biological Inspiration](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/reccqqnP7PstMuQ0j) · [doi](https://doi.org/10.1117/12.886413)  
  Earthworm peristalsis cited as soft-robot inspiration; no Earthworm entry in the animal vocabulary.
- [Roden-Reynolds 2015 — Hip proprioceptive feedback influences the control of mediolateral stability during human ](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recfdwgo7yvRq9Lt7) · [doi](https://doi.org/10.1152/jn.00551.2015)  
  Mediolateral foot-placement effects of hip proprioception have no corresponding feedback vocabulary item.
- [Rossignol 2011 — Neural Control of Stereotypic Limb Movements](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recgzIOzQxO3lmhT2) · [doi](https://doi.org/10.1002/cphy.cp120105)  
  Abstract is a section list only (no substantive text); species not stated, so animals left sparse and pathway fields unpopulated.
- [Rossignol and Dubuc 1994 — Spinal pattern generation](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/reclVERvbshs0fFlf) · [doi](https://doi.org/10.1016/0959-4388(94)90139-2)  
  Chick-embryo and tadpole-embryo preparations have no list entry; covered under 'Vertebrates' alongside Rat and Lamprey.
- [Sakaguchi 2006 — Instability of synchronized motion in nonlocally coupled neural oscillators](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recOlalJjpzPiMBKO) · [doi](https://doi.org/10.1103/physreve.73.031907)  
  Pure nonlinear-dynamics paper (chimera states); no animal or afferent content. Likely included as oscillator-coordination background; scope decision for Ben.
- [Schmitt 2000 — Mechanical models for insect locomotion: dynamics and stability in the horizontal plane - ](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recHnESCpMQBAkl1J) · [doi](https://doi.org/10.1007/s004220000180)  
  Tagged 'Biomechanically mediated preflexive feedback' because the abstract demonstrates passive-mechanics stability matching cockroach gait and forces, though it never uses the word 'feedback'; confirm.
- [Shamaei 2013 — Estimation of Quasi-Stiffness of the Human Knee in the Stance Phase of Walking](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rec4sOUJC8ZCeSYgW) · [doi](https://doi.org/10.1371/journal.pone.0059993)  
  Tagged 'Biomechanically mediated preflexive feedback' for the stance-phase knee quasi-stiffness (linear moment-angle relation), although the abstract does not use preflex framing - drop the tag if that is too liberal.
- [Steuer and Guertin 2019 — Central pattern generators in the brainstem and spinal cord: an overview of basic principl](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recV3anD1eqvb5qrV) · [doi](https://doi.org/10.1515/revneuro-2017-0102)  
  Abstract mentions in vitro preparations from invertebrate species and primates; neither 'Invertebrates' nor 'Primates' is in the allowed animals list, so only Vertebrates tagged.
- [Szczecinski 2017 — Design process and tools for dynamic neuromechanical models and robot controllers](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/reclFa8gRAcl3H377) · [doi](https://doi.org/10.1007/s00422-017-0711-4)  
  is_model marked true for the design-process/controller contribution; Ben may prefer classifying it as methods/tools rather than a model.
- [Tazerart 2008 — The Persistent Sodium Current Generates Pacemaker Activities in the Central Pattern Genera](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recBqaMMPYjKT0toi) · [doi](https://doi.org/10.1523/jneurosci.1437-08.2008)  
  Species is given only as 'neonatal rodents', listed here under 'Mammals'; narrow to Rat or Mice from the full text if the exact species matters.
- [Ueda 2015 — Multistate network model for the pathfinding problem with a self-recovery property](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recj6eKOmKGpYISSR) · [doi](https://doi.org/10.1016/j.neunet.2014.08.008)  
  Appears off-topic for a sensory afferent database (abstract graph pathfinding; no locomotion, reflex, or afferent content) - Ben should decide whether to keep.
- [Ullner 2016 — Self-Sustained Irregular Activity in an Ensemble of Neural Oscillators](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recmjUmtKNDBVnNir) · [doi](https://doi.org/10.1103/physrevx.6.011015)  
  Abstract is a single summary sentence, too thin for detailed curation; notes and is_model rest on the title plus that sentence.
- [Vogelstein 2008 — A Silicon Central Pattern Generator Controls Locomotion in Vivo](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recwZPgwSdttusQDb) · [doi](https://doi.org/10.1109/tbcas.2008.2001867)  
  Abstract does not state the animal species or afferent types (only 'sensory feedback recorded from the animal's legs'); no feedback vocabulary match for closed-loop neuromorphic sensory input.
- [Zehr and Stein 1999 — What functions do reflexes serve during human locomotion](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recgQ4H7oES3zg5BT) · [doi](https://doi.org/10.1016/s0301-0082(98)00081-1)  
  Stretch- and load-receptor afferents discussed functionally but not typed in the abstract (load could be Ib or plantar cutaneous); bird flight also mentioned, with no vocabulary animal.
- [van Arkel 2015 — The capsular ligaments provide more hip rotational restraint than the acetabular labrum an](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recApiY0mIaDT9YNn) · [doi](https://doi.org/10.1302/0301-620x.97b4.34638)  
  Orthopaedic cadaver-biomechanics study, not a sensory-afferent paper; relevant to SADb only as passive hip rotational-restraint data for biomechanical models, so Ben should decide on retention.

