# Round-2 SADb: database-first opening, DOI qualification, methods compression.
import os

p = os.path.join(
    r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation"
    r"\Notes\advisor_revision_20261022\proposals\50-database.tex"
)
s = open(p, encoding="utf-8").read()
fails = []

def rep(old, new, tag):
    global s
    n = s.count(old)
    if n != 1:
        fails.append((tag, n))
        return
    s = s.replace(old, new)

# 1. database-first opening
rep(
    "The synthetic nervous system of Chapter~\\ref{ch:sns} is assembled from physiological parts: Ia, Ib, and II afferents, the reflex pathways they drive, and the spinal circuits those pathways modulate. Building that controller requires numbers. For each sensory channel, the designer has to know what the receptor encodes, over what ranges of length, velocity, and force it operates and with what receptor dynamics, and what the afferent does at the spinal end: which interneuron populations it contacts, whether the action is excitatory or inhibitory, whether the target is ipsilateral, contralateral, or both, and in which phase of the step cycle the pathway is open. Element by element, that specification was established experimentally, one pathway and one preparation at a time, across a century of work in cat, mouse, rat, and human preparations, and much of it remains unknown or contested outside the pathway-rich core. A controller designer is not the only reader stranded by how those results are distributed over an enormous primary literature: an animal researcher planning the next experiment wants to know which preparations established a pathway property, and how it compares across species, and a roboticist wants to know which receptor encoding a sensor should imitate. Those readers ask different questions of the same literature, and the literature answers none of them without being gathered, structured, and searchable first.",
    "The Sensory Afferent Database (SADb) gathers into one searchable resource the physiological specification from which the synthetic nervous system of Chapter~\\ref{ch:sns} is assembled: what each class of Ia, Ib, and II afferent encodes, over what ranges of length, velocity, and force it operates, which interneuron populations it contacts, whether the action is excitatory or inhibitory, ipsilateral, contralateral, or both, and in which phase of the step cycle the pathway is open. That specification was established experimentally, one pathway and one preparation at a time, across a century of work in cat, mouse, rat, and human preparations, and much of it remains unknown or contested outside the pathway-rich core, scattered over an enormous primary literature. A controller designer, an animal researcher planning the next experiment, and a roboticist choosing what a sensor should imitate ask different questions of that same literature, and it answers none of them without being gathered, structured, and searchable first.",
    "S1-opening",
)

rep(
    "The Sensory Afferent Database (SADb) is the structured afferent resource built to answer those questions. I built it within the AARL's NeuroNex Communication, Coordination, and Control in Nervous Systems (C3NS) effort and presented it at Neuroscience 2023 \\citep{bolen_sensory_2023}, and later presented the AI-assisted triage study that expanded it at Neuroscience 2024 \\citep{bolen_using_2024}. Its content is organized as five linked tables (Papers, Feedback, Models, Review Papers, and Recorders), with the DOI as the unique identifier of every record.",
    "I built it within the AARL's NeuroNex Communication, Coordination, and Control in Nervous Systems (C3NS) effort and presented it at Neuroscience 2023 \\citep{bolen_sensory_2023}, and later presented the AI-assisted triage study that expanded it at Neuroscience 2024 \\citep{bolen_using_2024}. Its content is organized as five linked tables (Papers, Feedback, Models, Review Papers, and Recorders), with the DOI as the unique identifier wherever a record carries one.",
    "S2-intro-merge",
)

# 2. assembly compression
rep(
    "The corpus was assembled from three sources that did not agree: my personal Zotero collection (538 items), the AARL group library's Biology subtree (911 items across 27 collections), and the Airtable base (99 papers). Matching record by record on DOI and normalized title showed that the personal collection was essentially a subset of the group library and that the Airtable base was a small curated subset rather than a mirror; the gap was imported in verified batches, with every written record re-read and checked after writing. The corpus then grew through citation mining in the style of Research Rabbit, newly saved references from the laboratory's Zotero libraries, and a bulk pass that imported 372 papers auto-classified from the laboratory's staged PDF collection, bringing the audited corpus to 943 papers, the full set against which every later decision in this chapter was made.",
    "The corpus was assembled from three sources that did not agree: my personal Zotero collection (538 items), the AARL group library's Biology subtree (911 items across 27 collections), and the Airtable base (99 papers). Matching record by record on DOI and normalized title showed that the personal collection was essentially a subset of the group library and that the Airtable base was a small curated subset rather than a mirror, and the gap was imported in verified batches, each written record re-read and checked after writing. The corpus then grew through citation mining, newly saved references, and a bulk pass that imported 372 papers auto-classified from the laboratory's staged PDF collection, bringing the audited corpus to 943 papers, the full set against which every later decision in this chapter was made.",
    "S3-assembly",
)

# 3. machine/human division compression
rep(
    "The curation was machine-assisted under a division of labor in which the machine drafts and the human rules; the comparative evaluation of the AI tools that motivated this division is reported in Section~\\ref{sec:sadb_triage} \\citep{bolen_using_2024}. The machine contributed scale: batch drafting from grounded text, 374 papers curated in fifteen batches, automated vocabulary, abstract-paste, and duplicate gates that rejected bad records before anything was written, and per-record verification against the grounding text after the write. I held every decision that required judgment or carried risk: the curation vocabulary and its amendments, the 164 flagged judgment calls, and every single deletion in the prune. A PDF-first rule recovered grounding for 47 of 55 flagged papers from the laboratory's own library, and a university-authenticated browser session recovered abstracts for 24 of the 70 papers that had no findable text; the rest stayed bare by rule.",
    "The curation was machine-assisted under a division of labor in which the machine drafts and the human rules; the comparative evaluation of the AI tools that motivated this division is reported in Section~\\ref{sec:sadb_triage} \\citep{bolen_using_2024}. The machine contributed scale: batch drafting from grounded text, automated vocabulary, abstract-paste, and duplicate gates that rejected bad records before anything was written, and per-record verification against the grounding text after the write. I held every decision that required judgment or carried risk: the curation vocabulary and its amendments, the 164 flagged judgment calls, and every single deletion in the prune. Papers for which no text could be found stayed bare by rule rather than guessed, with grounding recovered where the laboratory's own library or a university-authenticated session could provide it.",
    "S4-division",
)

# 4. prune compression
rep(
    "Pruning followed, in three stages, each approved by me before execution, and together the stages removed 436 records: 382 papers plus 54 orphaned review, model, and pathway records that lost their anchors to the paper cuts. First, 95 records left as off-topic or orphaned, the satellite topics named above plus duplicates and isolated records stranded by an earlier strict flagging pass. Second, a cut list of 290 papers assembled from explicit rules for duplicates, monographs, and conference abstracts, and from a weighted rubric ranking the remainder by citation count, age, rating, network isolation, and record thinness, while 351 protected records (foundational authors, feedback-linked papers, highly cited work, and pre-1940 classics) were excluded from consideration. Third, a final round removed four papers and 47 review records whose links had been emptied by the earlier cuts. The prune was deliberate on two grounds: the working base had swollen past the laboratory's free-plan record cap, and a curated resource that answers a design question in ten records beats a haystack of a thousand. No record was destroyed: every deletion was snapshotted, and the packaged explorer carries the 382 cut papers under an archived flag, so the audited corpus remains reconstructible in full.",
    "Pruning then removed 436 records, 382 papers plus 54 orphaned review, model, and pathway records, in three stages: off-topic satellites, duplicates, and orphans first; a rule-built cut list of 290 papers next, ranked by a weighted rubric over citation count, age, rating, network isolation, and record thinness, with 351 protected records (foundational authors, feedback-linked papers, highly cited work, and pre-1940 classics) excluded from consideration; and a final sweep of records whose links the earlier cuts had emptied. The prune was deliberate on two grounds: the working base had swollen past the laboratory's free-plan record cap, and a curated resource that answers a design question in ten records beats a haystack of a thousand. No record was destroyed; every deletion was snapshotted, and the packaged explorer carries the 382 cut papers under an archived flag, so the audited corpus remains reconstructible in full.",
    "S5-prune",
)

with open(p, "w", encoding="utf-8", newline="") as f:
    f.write(s)

if fails:
    print("FAILED:", fails)
else:
    print("all SADb round-2 replacements applied")
print("words now:", len(s.split()))
