# Batch 3 drafts (queue rows 11-20) — 2026-09-12, laptop session

Write order: resolve 10 Paper record ids (DOI search) -> create 3 Review Papers + 1 Models
-> update the 10 papers (Notes/Animals/Feedback) -> verify by reopening each -> log.

Feedback vocab ids in use:
- Ib stance to swing = recVAxf4i8eprMZ5X
- Ia or II stance to swing = recfRtpyUZxlj2uPR
- Ib disynaptic inhibition = recsaUPV6EKpf08tO
- trochanteral campaniform sensilla load→MN magnitude = recZCTtDpMMOAHRxd

New tag records: Review "Whelan 1996", "Tresch and Jarc 2009", "Frigon and Rossignol 2006";
Models "Rybak 2013".

---

1. LWWSHR63 — Sivertsen et al. 2016, Pontine reticulospinal projections in the neonatal mouse (J Comp Neurol, 10.1002/cne.23904; PMID 26400815)
Notes: "Retrograde tracing in the neonatal mouse shows the pontine reticulospinal projection consists of segregated ipsi- and contralaterally projecting populations (3:1 ipsilateral) that are predominantly excitatory (GAD-negative), with distinct brainstem and spinal axon trajectories; the larger size and greater number of the uncrossed population accounts for its more reliable transmission to spinal motoneurons. This is the anatomical substrate for the powerful descending reticulospinal drive assumed in spinal locomotor models."
Animals: Mice
Feedback: none (descending tract anatomy)
Grounding: PubMed abstract (Zotero local search failed on slashed DOI; item exists in personal lib)

2. ACANM79H — Whelan 1996, Control of locomotion in the decerebrate cat (Prog Neurobiol 49:481-515, 10.1016/0301-0082(96)00028-7; PMID 8895997)
Review Papers record: "Whelan 1996" (paper link)
Notes: "In the decerebrate cat, locomotion is initiated by the mesencephalic locomotor region acting through the medial medullary reticular formation and the ventrolateral funiculus, and the phase transitions are set by identified afferents: group I Golgi tendon organ afferents prolong stance, while length- and velocity-sensitive afferents from extensor muscles (hip) signal leg extension to permit swing — unloaded and extended is the swing-permitting state."
Animals: Cat
Feedback: Ib stance to swing; Ia or II stance to swing (mirrors Pearson 1995 pilot pairing; both mechanisms explicit in abstract)
Grounding: PubMed abstract PMID 8895997

3. 8AUBC3V7 — Kaliyamoorthy et al. 2005, Force Sensors in Hexapod Locomotion (Int J Robot Res, 10.1177/0278364905055381)
Notes: "Finite-element modeling of the cockroach leg shows that loads are sensed by leg force sensors located close to the body, where specific force vectors (body load versus propulsion) can be discriminated, and that this information is used in positive load feedback to regulate walking — design principles directly transferable to force-controlled legged robots."
Animals: Insects; Cockroach
Feedback: trochanteral campaniform sensilla load→MN magnitude (recZCTtDpMMOAHRxd) — JUDGMENT CALL: abstract says "positive load feedback" via insect leg force sensors (Zill-group CS work); flag if Ben wants own-work reviews tagged differently
Grounding: Zotero personal abstract (876 chars)

4. YSB9TKJX — Tresch and Jarc 2009, The case for and against muscle synergies (Curr Opin Neurobiol, 10.1016/j.conb.2009.09.002; PMID 19828310)
Review Papers record: "Tresch and Jarc 2009" (paper link)
Notes: "Tresch and Jarc 2009 weigh the evidence for and against muscle synergies as a hypothesized output level of the CNS, concluding that the field's task is to distinguish flexible combinations of muscle groups from alternative control variables (individual muscles, units, kinematics) rather than to accept or reject modularity wholesale."
Animals: (blank — review spans preparations not enumerable from abstract; flag in log)
Feedback: none
Grounding: Zotero personal abstract + PubMed

5. JU4WBGTJ — Ross and Nichols 2009, Heterogenic feedback between hindlimb extensors in the spontaneously locomoting premammillary cat (J Neurophysiol 101:184-197, 10.1152/jn.90338.2008; PMID 19005003)
Notes: "During spontaneous treadmill stepping in the premammillary decerebrate cat, force-dependent heterogenic inhibition between hindlimb extensors persists (quadriceps onto gastrocnemius, gastrocnemius onto plantaris/FHL) but distal-onto-proximal inhibition is weaker than during the crossed-extension reflex, yielding a proximal-to-distal gradient of Ib inhibition that supports interjoint coordination and limb stability."
Animals: Cat
Feedback: Ib disynaptic inhibition (recsaUPV6EKpf08tO)
Grounding: PubMed abstract PMID 19005003

6. QU973HAA — Frigon and Rossignol 2006, Functional plasticity following spinal cord lesions (Prog Brain Res 157:231-260, 10.1016/s0079-6123(06)57016-5; PMID 17167915)
Review Papers record: "Frigon and Rossignol 2006" (paper link)
Notes: "After spinal cord injury, reflex pathways caudal to the lesion are initially depressed by low motoneuron excitability and then recover — sometimes to exaggeration (spasticity) — and in spinal cats step training normalizes transmission in simple reflex pathways, suggesting that the modified afferent inflow must itself be normalized for a stable locomotor rhythm to be re-expressed."
Animals: Cat
Feedback: none
Grounding: PubMed abstract PMID 17167915

7. YJVT4HJU — Horchler et al. 2004, Highly mobile robots that run and jump (no DOI)
Notes: "Mini-WHEGS robots abstract cockroach locomotion principles into 9-cm four-wheel-leg vehicles using an alternating diagonal gait that run over 10 body lengths per second, climb obstacles taller than their leg length, and add a jumping mechanism for larger obstacles — an example of insect-derived mechanical intelligence transferable to bipedal foot placement hardware."
Animals: (blank — robot paper; note Zotero record carries a PATENT-style abstract "The present invention relates to..."; flag record-type question for Ben)
Feedback: none
Grounding: Zotero personal abstract only (970 chars) — patent-style text; no better text found without DOI

8. DQLYCBH2 — Rybak et al. 2013, Modelling genetic reorganization in the mouse spinal cord affecting left-right coordination during locomotion (J Physiol 591:5491-5508, 10.1113/jphysiol.2013.261115; PMID 24081162)
Models record: "Rybak 2013" (paper link; paper Models 2 auto-links)
Notes: "A computational model of the left-right commissural circuitry (inhibitory and excitatory commissural interneuron populations plus an EphA4-positive subpopulation) reproduces the gait phenotypes of axon-guidance mutations — EphA4 knockout converts to synchronized hopping via crossed excitation, Netrin-1 knockout by loss of contralateral inhibition, DCC knockout by loss of both — and shows amplified inhibition restores alternation in EphA4 and DCC knockouts but not Netrin-1 knockouts."
Animals: Mice
Feedback: none (CPG/genetic, not sensory)
Grounding: PubMed abstract PMID 24081162

9. LNCVWZ8E — Pearson 2000, Neural adaptation in the generation of rhythmic behavior (Annu Rev Physiol 62:723-753, 10.1146/annurev.physiol.62.1.723; PMID 10845109)
Review Papers record: "Pearson 2000" (paper link)
Notes: "Central pattern generators are extremely flexible: neuromodulators, central commands and afferent signals reshape the cellular and synaptic properties of the neurons and the coupling between populations, and afferent feedback from limb proprioceptors drives the long-term adaptation that keeps motor output matched to changed body mechanics and persistent performance errors."
Animals: (blank from abstract — review spans cat/insect examples not enumerable from abstract; flag)
Feedback: none
Grounding: PubMed abstract PMID 10845109

10. WA46PM2Y — Duysens et al. 2013, The flexion synergy, mother of all synergies and father of new models of gait (Front Comput Neurosci 7:14, 10.3389/fncom.2013.00014; PMID 23494365)
Notes: "The flexion synergy of Sherrington's flexor reflex, the withdrawal-reflex modules of Schouenborg, and the neonatal locomotor modules of Dominici largely overlap, and the end-of-stance facilitation of the flexion synergy — broadly afferent-facilitated but load-suppressed — points to a flexor burst generator as the core of an asymmetric CPG model whose afferent gating has already been implemented in walking bipedal robots."
Animals: Human
Feedback: Ib stance to swing (recVAxf4i8eprMZ5X) — JUDGMENT CALL: abstract states flexor burst activation is "suppression by load afferent input", i.e., load receptors gate the stance-to-swing transition; flag if Ben prefers a dedicated record
Grounding: PubMed abstract PMID 23494365
