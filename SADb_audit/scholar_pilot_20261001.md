# Scholar-via-browser pilot (2026-10-01) — grounding WITHOUT API calls

Ben's insight: the in-app browser holds his logged-in Google Scholar (PSU
library) + Airtable sessions, so article lookups can run through the browser
instead of the rate-limited public API. Pilot: 3 no-text papers, exact-phrase
searches at human pace (~2-3 s between queries), snippets read from the
DOM. All three returned abstract-grade text that OpenAlex + Europe PMC
never had.

| Record | Paper | Scholar result |
|---|---|---|
| [rec1fq8mn2oqPeLiQ](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rec1fq8mn2oqPeLiQ) | Armstrong 1986, Supraspinal contributions (Prog Neurobiol) | exact match; snippet = abstract opening (Muybridge/Marey/Sherrington/Graham Brown lineage) |
| [rec2h1B7gvrjFk4I4](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rec2h1B7gvrjFk4I4) | Orlovsky 1972, vestibulospinal neurons (Brain Res) | exact match; snippet carries real results text (stance-phase peak discharge, cerebellum-intact modulation) |
| [rec2BZRpatIkJvjp0](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rec2BZRpatIkJvjp0) | Sedlackova 2020, SNS optomotor response (Springer) | exact match; snippet = abstract (MantisBot, wide-field vision, dynamical neural model) |

## Workflow for the remaining ~107 no-text papers
1. For each `curation_queue/notext.json` entry: exact-phrase Scholar search
   (2-3 s pauses, batches of ~10 per session — Scholar rate-limits bots).
2. Snippet (or the linked PDF/landing page where PSU access opens it) becomes
   the grounding text in `export/grounding.json`.
3. Curate + write once the API budget recovers (or via the Airtable UI tab
   for small batches — also browser, also free of API calls).

Snippets are stored below for the three pilots (verify against the paper when
curating — Scholar snippets can truncate mid-sentence).

## Armstrong 1986
Interest in analyzing the mechanisms which control locomotion in mammals (and
for that matter in other vertebrates and in invertebrates) has never been more
intense than at present. A promising start was made in the late nineteenth and
early twentieth centuries by such luminaries as Muybridge, Marey and
Phillipson who analysed movements and Sherrington and Graham Brown who studied
the underlying neural mechanisms (fo… [truncated]

## Orlovsky 1972
The activity of vestibulospinal neurons giving axons to the lumbosacral
spinal cord was recorded during locomotion (walking and running on the
treadmill) in mesencephalic and thalamic cats. The overall activity of most
neurons increases to a considerable degree during locomotion, and periodic
alternations of this activity in relation to the locomotor cycle (modulation)
were observed in cats with intact cerebellum. The peak discharge usually
occurs at the beginning of the stance phase of the ipsi… [truncated]

## Sedlackova 2020
We seek to increase the sophistication of our insect-like hexapod robot
MantisBot's visual system. We assembled and tested a benchtop robotic testbed
with which to test our dynamical neural model of the insect visual system.
Here we specifically model wide-field vision and the optomotor response. The
system is composed of a Raspberry Pi with a camera outfitted with a 360°
lens. The camera sits on a motorized turntable, which represents the
"robot". Above the turntable sits another motorized syst… [truncated]
