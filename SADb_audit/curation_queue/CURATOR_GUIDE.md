# CURATOR GUIDE — SADb curation layer v2 (2026-09-28)

You are curating Sensory Afferent Database papers. You receive a batch file
(`curation_queue/queue_NN.json`, a JSON array of papers with `id`, `title`,
`primary` author, `year`, `doi`, `abstract`, `topic`). You produce ONE output
file `curation_out/batch_NN.json` — a JSON array, SAME ORDER, SAME `id`s.

## Grounding rule (absolute) — PDF FIRST (Ben, 2026-10-01)
Ground on the paper's FULL TEXT whenever a PDF exists — Ben's rule: "You need
to look at the actual pdf if it is attached." Source order:
1. `curation_queue/pdf_grounding/<recordId>.txt` (extracted from the local
   Zotero copy via `pdf_grounding.py` — run it for your batch's papers first)
2. The paper's `abstract` from the batch file (fallback when no PDF exists)
3. Scholar snippet in `export/grounding_scholar.json` (second fallback)
NEVER write fields from model memory. If none of the three exist, leave the
field empty — never guess. Note: OpenAlex abstracts are occasionally attached
to the WRONG record (e.g. Carlson-Kuhta 'Forms II' carried the 'Forms III'
abstract) — when a PDF is present it ALWAYS outranks the abstract.

## Fields to produce per paper

- **notes** (string, required, 20–1500 chars): 1–3 sentences of distilled
  insight usable verbatim in a dissertation with `\citep{}` — mechanism
  statements, what the paper established, and where useful how it feeds
  models. NOT an abstract paraphrase: no background filler, no "we
  investigated", no verbatim abstract sentences. Gold style examples:
  - "Loading of ankle extensors during stance inhibits flexor burst
    generation through Ib afferent pathways, providing a load-dependent
    gate on the swing-to-stance transition." (mechanism + significance)
  - "Established that deafferented chronic spinal cats retain fictive
    locomotion, localizing rhythm generation to spinal networks
    independent of sensory input."

- **animals** (array, exact names from this list ONLY — all preparations the
  paper actually discusses, empty if none/unclear):
  Cat, Dog, Frog, Hexapod, Human, Insects, **Invertebrates** (Ben 2026-10-02),
  Lamprey, Mammals, Mice, Rat, Salamander, Stick Insect, Vertebrates, Zebrafish,
  Arthropods, Cockroach, Turtle. ("Mice", never "Mouse".)

- **afferents** (array, exact names ONLY, empty if the abstract does not
  specify afferent types): "Type 1 (legacy)" (old papers using the type 1/2
  numbering), "Ia", "Ib", "II", "III/IV", "Mechanosensory", "Cutaneous",
  "Heat", "Nociceptive", "Flexor reflex afferents". Muscle spindle primary
  = Ia; spindle secondary = II; Golgi tendon organ = Ib; group III/IV =
  metabolic/fatigue-nociceptive fine afferents; hair plates, campaniform
  sensilla, chordotonal organs = Mechanosensory.

- **feedback** (array of EXACT names from the vocabulary below; stay
  conservative — only when the abstract explicitly demonstrates/describes
  the pathway; empty otherwise; NEVER invent names):
  "Type 1 swing to stance", "type 1 stance to swing", "Ia stance to swing",
  "Ia swing to stance", "Ia or II stance to swing", "Ia or II swing to stance",
  "Ib stance to swing", "Ib swing to stance", "Ib disynaptic excitation",
  "Ib disynaptic inhibition", "Ib inhibition", "Ib excitatory",
  "Ib contralateral inhibition", "Ia presynaptic inhibition",
  "Ia monosynaptic excitation", "Ia monosynaptic", "Ia disynaptic inhibition",
  "Ia reciprocal inhibition", "Ia inhibitory", "II inhibitory",
  "type II excitatory", "Cutaneous stance modification", "Cutaneous flexor excitation",
  "Fictive locomotion without sensory feedback", "Total afferent inhibition",
  "group III/IV fatigue", "Mechanosensory monosynaptic excitation",
  "Biomechanically mediated preflexive feedback",
  "trochanteral hair plate multi-synaptic excitation",
  "trochanteral hair plate multi-synaptic inhibition",
  "trochanteral campaniform sensilla load signals to adjust MN magnitude",
  "Chordotonal organ multi-synaptic excitation",
  "Chordotonal organ multi-synaptic inhibition",
  "Ipsilateral excitatory - swimming edge cell", "Cross inhibitory - swimming edge cell",
  "large diameter spinal afferent stimulation".
  If the paper demonstrates a pathway that clearly fits NONE of these, do NOT
  force it: leave feedback empty and put a one-line description in `flags`,
  **composed in Ben's naming grammar (2026-10-02)**:
  `afferent + (inhibit|excite) + (contralateral|ipsilateral|agonist|antagonist) + joint + phase`
  — e.g. "ankle Ib feedback to hip during stance" (his Ekeberg-Pearson
  example), "lateral postural reflex" (his Karayannidou naming). Standalone
  phase terms ("stance", "swing", "stance to swing", "swing to stance") are
  combinable primitives. REVERSE direction (phase/state MODULATING afferent
  gain, as in Akazawa 1982) is its own dimension: "phase → afferent gain",
  not a forward pathway.

- **animal_study** (string, ≤1000 chars): how the paper's information could
  seed an animal study, or how its hypothesis could be tested ON an animal
  (e.g. which perturbation/lesion/recording would test the claim). "" if
  genuinely N/A (e.g. pure simulation papers with no testable animal claim).

- **robot_sim** (string, ≤1000 chars): which neural pathway, synapse type,
  or loss-of-function result from the paper's animal work would let a
  simulated or robotic model test its hypotheses — i.e. what to implement or
  ablate in a model to reproduce/predict the paper's finding. "" if N/A.

- **is_model** (boolean) + **model_name** (e.g. "Taga 1995" style
  "Surname Year"): true ONLY when the paper's CORE contribution is a
  computational/conceptual model (not merely using one).

- **is_review** (boolean) + **review_name** ("Surname Year"): true for
  literature reviews / review essays (not primary research, not monographs).

- **flags** (string, optional): anything Ben should decide — proposed new
  Feedback vocabulary, ambiguity, abstract too thin for a specific field.

## Output schema (exact keys, same order as input)
```json
[{"id": "recXXXXXXXXXXXXXX", "notes": "...", "animals": ["Cat"],
  "afferents": ["Ib"], "feedback": ["Ib stance to swing"],
  "animal_study": "...", "robot_sim": "...",
  "is_model": false, "model_name": "", "is_review": false, "review_name": "",
  "flags": ""}]
```

## Hard rules
- Output must be VALID JSON (no trailing commas, no comments).
- Never drop or reorder papers; never change ids.
- One output file only; no Airtable/web/network access.
