#!/usr/bin/env python
# Build Airtable payloads for batch 3 (queue rows 11-20) from the reviewed drafts.
# Field ids (from README): Notes=fld3gPiUIKn26N6ji, Animals=fld1x2BXLKIdA2dCw,
# Feedback=fldK2H6RfaSdLLLM7. Record ids resolved at write time by DOI filter.
import json, os

FB = {
    "ib_sts_swing": "recVAxf4i8eprMZ5X",
    "ia2_sts_swing": "recfRtpyUZxlj2uPR",
    "ib_dis_inh": "recsaUPV6EKpf08tO",
    "cs_load_mn": "recZCTtDpMMOAHRxd",
}

ROWS = [
    dict(doi="10.1002/cne.23904", key="LWWSHR63",
         notes="Retrograde tracing in the neonatal mouse shows the pontine reticulospinal projection consists of segregated ipsi- and contralaterally projecting populations (3:1 ipsilateral) that are predominantly excitatory (GAD-negative) with distinct brainstem and spinal axon trajectories; the larger size and greater number of the uncrossed population accounts for its more reliable transmission to spinal motoneurons. This is the anatomical substrate for the powerful descending reticulospinal drive assumed in spinal locomotor models.",
         animals=["Mice"], fb=[]),
    dict(doi="10.1016/0301-0082(96)00028-7", key="ACANM79H",
         notes="In the decerebrate cat, locomotion is initiated by the mesencephalic locomotor region acting through the medial medullary reticular formation and the ventrolateral funiculus, and the phase transitions are set by identified afferents: group I Golgi tendon organ afferents prolong stance, while length- and velocity-sensitive afferents from extensor muscles signal leg extension to permit swing - unloaded and extended is the swing-permitting state.",
         animals=["Cat"], fb=[FB["ib_sts_swing"], FB["ia2_sts_swing"]]),
    dict(doi="10.1177/0278364905055381", key="8AUBC3V7",
         notes="Finite-element modeling of the cockroach leg shows that loads are sensed by leg force sensors located close to the body, where specific force vectors (body load versus propulsion) can be discriminated, and that this information is used in positive load feedback to regulate walking - design principles directly transferable to force-controlled legged robots.",
         animals=["Insects", "Cockroach"], fb=[FB["cs_load_mn"]]),
    dict(doi="10.1016/j.conb.2009.09.002", key="YSB9TKJX",
         notes="Tresch and Jarc 2009 weigh the evidence for and against muscle synergies as a hypothesized output level of the CNS, concluding that the field's task is to distinguish flexible combinations of muscle groups from alternative control variables (individual muscles, units, kinematics) rather than to accept or reject modularity wholesale.",
         animals=[], fb=[]),
    dict(doi="10.1152/jn.90338.2008", key="JU4WBGTJ",
         notes="During spontaneous treadmill stepping in the premammillary decerebrate cat, force-dependent heterogenic inhibition between hindlimb extensors persists (quadriceps onto gastrocnemius, gastrocnemius onto plantaris/FHL) but distal-onto-proximal inhibition is weaker than during the crossed-extension reflex, yielding a proximal-to-distal gradient of Ib inhibition that supports interjoint coordination and limb stability.",
         animals=["Cat"], fb=[FB["ib_dis_inh"]]),
    dict(doi="10.1016/s0079-6123(06)57016-5", key="QU973HAA",
         notes="After spinal cord injury, reflex pathways caudal to the lesion are initially depressed by low motoneuron excitability and then recover - sometimes to exaggeration (spasticity) - and in spinal cats step training normalizes transmission in simple reflex pathways, suggesting that the modified afferent inflow must itself be normalized for a stable locomotor rhythm to be re-expressed.",
         animals=["Cat"], fb=[]),
    dict(doi="", key="YJVT4HJU", title="Highly mobile robots that run and jump",
         notes="Mini-WHEGS robots abstract cockroach locomotion principles into 9-cm four-wheel-leg vehicles using an alternating diagonal gait that run over 10 body lengths per second, climb obstacles taller than their leg length, and add a jumping mechanism for larger obstacles - an example of insect-derived mechanical intelligence transferable to bipedal foot placement hardware.",
         animals=[], fb=[]),
    dict(doi="10.1113/jphysiol.2013.261115", key="DQLYCBH2",
         notes="A computational model of the left-right commissural circuitry (inhibitory and excitatory commissural interneuron populations plus an EphA4-positive subpopulation) reproduces the gait phenotypes of axon-guidance mutations - EphA4 knockout converts to synchronized hopping via crossed excitation, Netrin-1 knockout by loss of contralateral inhibition, DCC knockout by loss of both - and shows amplified inhibition restores alternation in EphA4 and DCC knockouts but not Netrin-1 knockouts.",
         animals=["Mice"], fb=[]),
    dict(doi="10.1146/annurev.physiol.62.1.723", key="LNCVWZ8E",
         notes="Central pattern generators are extremely flexible: neuromodulators, central commands and afferent signals reshape the cellular and synaptic properties of the neurons and the coupling between populations, and afferent feedback from limb proprioceptors drives the long-term adaptation that keeps motor output matched to changed body mechanics and persistent performance errors.",
         animals=[], fb=[]),
    dict(doi="10.3389/fncom.2013.00014", key="WA46PM2Y",
         notes="The flexion synergy of Sherrington's flexor reflex, the withdrawal-reflex modules of Schouenborg, and the neonatal locomotor modules of Dominici largely overlap, and the end-of-stance facilitation of the flexion synergy - broadly afferent-facilitated but load-suppressed - points to a flexor burst generator as the core of an asymmetric CPG model whose afferent gating has already been implemented in walking bipedal robots.",
         animals=["Human"], fb=[FB["ib_sts_swing"]]),
]

# New tag records to create first (then link from the papers)
NEW_REVIEWS = [
    dict(name="Whelan 1996", paper_doi="10.1016/0301-0082(96)00028-7"),
    dict(name="Tresch and Jarc 2009", paper_doi="10.1016/j.conb.2009.09.002"),
    dict(name="Frigon and Rossignol 2006", paper_doi="10.1016/s0079-6123(06)57016-5"),
]
NEW_MODELS = [
    dict(name="Rybak 2013", paper_doi="10.1113/jphysiol.2013.261115"),
]

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "batch3")
os.makedirs(out, exist_ok=True)
with open(os.path.join(out, "payloads.json"), "w", encoding="utf-8") as f:
    json.dump(dict(rows=ROWS, new_reviews=NEW_REVIEWS, new_models=NEW_MODELS), f, indent=1, ensure_ascii=False)
print("wrote", os.path.join(out, "payloads.json"))
