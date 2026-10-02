# Prune review — ROUND 2 (Ben's 2026-10-01 patterns) — nothing executed

Patterns from your message, applied to the 855 live papers (row numbers
in the grid don't map to records, so these are content matches on
title + curation note; each entry links to its Airtable record).

## CS-style neural networks (not biologically inspired) — 3 papers

Prediction/classification/optimization NN usage with no biological framing — your rows-736-759 / 'computer-science neural networks' pattern.

- [Hilts 2000 — Emulating Balance Control Observed in Human Test Subjects with a Neural Network](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recV3YQtMTF5hreIm) · [doi](https://doi.org/10.1007/978-3-319-95972-6_21)
- [Kim 2021 — Multifrequency Hebbian plasticity in coupled neural oscillators](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recXW34i5WwjznlgT) · [doi](https://doi.org/10.1007/s00422-020-00854-6)
- [Vogels 2005 — NEURAL NETWORK DYNAMICS](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recpiClFNWFTFnOCA) · [doi](https://doi.org/10.1146/annurev.neuro.28.061604.135637)
## delete


## Traditional control (no neural basis) — 5 papers

PID/MPC/robust/adaptive controller papers — your rows 844/854/827/843 pattern. NOTE: Park 2003 and Wang 2020 are human postural-control EXPERIMENTS (biological subjects, EMG/force data) — auto-matched on 'control' wording but likely IN scope; Hyun 2014 is a classical hierarchical controller using proprioceptive signals — the robotics-boundary case you described. Giesseler and Wang 2013 are pure aerospace control — clear prunes.

- [Giesseler 2012 — Model Predictive Control for Gust Load Alleviation](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recNxBNmZRvcI2Sfr) · [doi](https://doi.org/10.3182/20120823-5-nl-3013.00049)
- [Hyun 2014 — High speed trot-running: Implementation of a hierarchical controller using proprioceptiv](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recHO6cXTO7rowMc5) · [doi](https://doi.org/10.1177/0278364914532150)
- [Park 2003 — Postural feedback responses scale with biomechanical constraints in human standing](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rechvvMAgwOmRzUTK) · [doi](https://doi.org/10.1007/s00221-003-1674-3)
- [Wang 2013 — Attitude and Altitude Controller Design for Quad-Rotor Type MAVs](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recJ7cgrCX7pokhZO) · [doi](https://doi.org/10.1155/2013/587098)
- [Wang 2020 — Standing Balance Experiment with Long Duration Random Pulses Perturbation](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recKffjcACq4E7Jfe) · [doi](https://doi.org/10.5281/zenodo.3631958)
## delete Giesseler and Wang only


## General neuron models (no biology) — 1 papers

Neuron-modeling math without a biological/preparation context.

- [Fairhurst 2010 — Observers for Canonic Models of Neural Oscillators](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recx6cU1RkrYsx1TX) · [doi](https://doi.org/10.1051/mmnp/20105206)
## delete

## Duplicate candidates — 12 pairs (your rows-275/276 hint)

Same primary author + similar title (or same year). AUTO-VERDICT: a bioRxiv DOI (10.1101/...) on one side = true preprint/published duplicate (prune the preprint); titles differing by part numbering (I/II), hindlimb/forelimb, or data/model = companion papers, KEEP BOTH.

- sim 0.97 · same year 2024 — companion papers (different studies — keep both) 
  - [Mari 2024 — Changes in intra- and interlimb reflexes from hindlimb cutaneous afferents after](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rec8MVmUOx1oDk02j) (10.1113/jp286151)
  - [Mari 2024 — Changes in intra- and interlimb reflexes from forelimb cutaneous afferents after](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recB7n3thv2EEgswz) (10.1113/jp286808)
  ## keep both. Use the full title name. Don't be lazy.

- sim 0.92 · same year 2000 — review: companion or duplicate?
  - [Schmitt 2000 — Mechanical models for insect locomotion: dynamics and stability in the horizonta](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recBBTuyfoRBF4qKO) (10.1007/s004220000181)
  - [Schmitt 2000 — Mechanical models for insect locomotion: dynamics and stability in the horizonta](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recHnESCpMQBAkl1J) (10.1007/s004220000180)
## keep both. Use full article name to avaoid confusion.

- sim 0.91 — **TRUE DUPLICATE (preprint vs published — prune the 10.1101 side)**
  - [Merlet 2021 — Cutaneous inputs from perineal region facilitate spinal locomotor activity and m](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recIpNrSMGJDUWbbu) (10.1002/jnr.24791)
  - [Merlet 2020 — Cutaneous inputs from perineal region facilitates and modulates spinal locomotor](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recZxoFgqdQh2aUJO) (10.1101/2020.07.29.226530)
## delete preprint, obviously

- sim 0.9 · same year 1998 — review: companion or duplicate?
  - [Procházka and Gorassini 1998 — Ensemble firing of muscle afferents recorded during normal locomotion in cats](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recXn5YzGAlqnE7K9) (10.1111/j.1469-7793.1998.293bu.x)
  - [Procházka and Gorassini 1998 — Models of ensemble firing of muscle spindle afferents recorded during normal loc](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recayW4BudURUfmSb) (10.1111/j.1469-7793.1998.277bu.x)
## keep both

- sim 0.86 · same year 2006 — review: companion or duplicate?
  - [Burkitt 2006 — A review of the integrate-and-fire neuron model: II. Inhomogeneous synaptic inpu](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recErzmL73q2PEhtl) (10.1007/s00422-006-0082-8)
  - [Burkitt 2006 — A Review of the Integrate-and-fire Neuron Model: I. Homogeneous Synaptic Input](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recgPmjGMrRBLOuA5) (10.1007/s00422-006-0068-6)
#keep both

- sim 0.78 · same year 1980 — review: companion or duplicate?
  - [Forssberg 1980 — The locomotion of the low spinal cat. II. Interlimb coordination](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recBgnEcuiN6cSJBL) (10.1111/j.1748-1716.1980.tb06534.x)
  - [Forssberg 1980 — The locomotion of the low spinal cat. I. Coordination within a hindlimb.](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/reczvPaUJrVXapwLB) (10.1111/j.1748-1716.1980.tb06533.x)
  ## keep both

- sim 0.75 · same year 2004 — review: companion or duplicate?
  - [Donelan and Pearson 2004 — Contribution of Force Feedback to Ankle Extensor Activity in Decerebrate Walking](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rec1YbgpV0dpRyTwr) (10.1152/jn.00325.2004)
  - [Donelan and Pearson 2004 — Contribution of sensory feedback to ongoing ankle extensor activity during the s](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recUDhhSVfIJ5cJg9) (10.1139/y04-043)


- sim 0.75 · same year 2006 — review: companion or duplicate?
  - [Rybak 2006 — Modelling spinal circuitry involved in locomotor pattern generation: insights fr](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recV7qGPJWiTDrPsM) (10.1113/jphysiol.2006.118711)
  - [Rybak 2006 — Modelling spinal circuitry involved in locomotor pattern generation: insights fr](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recy0oLAVRRLclNXl) (10.1113/jphysiol.2006.118703)
## Keep both. Look how important they are ffs. One is feed-forward, one is feedback. They both probably have a lot of citations.

- sim 0.74 — review: companion or duplicate? 
  - [Perreault 1993 — Activity of medullary reticulospinal neurons during fictive locomotion.](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recjsiiIz1XKKqefz) (10.1152/jn.1993.69.6.2232)
  - [Perreault 1994 — Microstimulation of the medullary reticular formation during fictive locomotion.](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recmVTrWXTVfSOWnt) (10.1152/jn.1994.71.1.229)
## I can't tell if they're unique. I assume so.

- sim 0.7 · same year 2005 — review: companion or duplicate?
  - [Quevedo 2005 — Intracellular analysis of reflex pathways underlying the stumbling corrective re](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/reclD5keu6AOXNPQL) (10.1152/jn.00176.2005)
  - [Quevedo 2005 — Stumbling corrective reaction during fictive locomotion in the cat.](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recoEJr11AdVO4jyw) (10.1152/jn.00175.2005)
## I can't tell since you don't have the papers in here

- sim 0.51 · same year 1991 — review: companion or duplicate?
  - [Pratt 1991 — Functionally complex muscles of the cat hindlimb. IV. Intramuscular distribution](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/rec5sjYrRbJ9t4e20) (10.1007/bf00229407)
  - [Pratt 1991 — Functionally complex muscles of the cat hindlimb. I. Patterns of activation acro](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recKW2KMXXHAHCpKi) (10.1007/bf00229406)
## keep. They're both different

- sim 0.5 · same year 2023 — review: companion or duplicate?
  - [Ramalingasetty 2023 — On All Fours: A 3D Framework to Study Closed-loop Control of Quadrupedal Mouse L](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recPnKgpWB1oLiyWJ) (10.1162/isal_a_00788)
  - [Ramalingasetty 2023 — An integrated neuromechanical model of the mouse to study neural control of loco](https://airtable.com/appMQTnobUNRytIp7/tblnnMrZszhboU4uD/recwy6OuthAArEsQ9) (10.18910/92326)
  ## I can't tell since you don't have the papers in here

