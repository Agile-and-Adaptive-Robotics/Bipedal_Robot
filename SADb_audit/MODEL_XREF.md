# MODEL_XREF — SADb rules/studies/models for neuromechanical modeling

**Standing rule (Ben, 2026-10-02):** any session building, tuning, or
debugging a neuromechanical model in **AnimatLab**, **SNS Simulink /
SNS_Simscape**, or the **MuJoCo-SNS toolbox** — on any machine — reads
this file FIRST and grounds pathway wiring/gains in the cited studies.
Regenerate after every curation batch:
`myo python SADb_audit/build_model_xref.py`. Deeper browsing:
`SADb_audit/knowledge_base/INDEX.md` and the SADb Explorer app
(`SADb_audit/app/sadb_app.html`, 943 papers incl. archived).

Corpus snapshot: 943 papers (382 app-only/archived).

## 1. Feedback-pathway RULES (the wiring vocabulary)

Each rule below is demonstrated by the linked papers (author year,
OpenAlex citations; [A] = archived from Airtable, still in the app).
When you implement a pathway, cite the demonstrating paper in comments.

- **Ib stance to swing** — 17 paper(s): Kiehn 2016, Duysens and Pearson 1980, Conway 1987, Pearson 1995
- **Cutaneous flexor excitation** — 16 paper(s): Sherrington 1910, Zehr and Stein 1999, Forssberg 1979, Forssberg 1977
- **Cutaneous stance modification** — 13 paper(s): Rossignol 2006, Forssberg 1979, Forssberg 1977, Duysens and Pearson 1976
- **Biomechanically mediated preflexive feedback** — 13 paper(s): Schmitt 2000, Koditschek 2004, Chiel 2009, Yakovenko 2004
- **Ib excitatory** — 12 paper(s): Pearson and Collins 1993, Procházka and Prochazka 2011, Prochazka 1997, Procházka 1997
- **Fictive locomotion without sensory feedback** — 12 paper(s): Cazalets 1992, Grillner and Zangger 1975, Guertin 2013, Pearson and Rossignol 1991
- **Ia or II stance to swing** — 10 paper(s): Kiehn 2016, Pearson 1995, Akay 2014, Kriellaars 1994
- **Ia monosynaptic excitation** — 7 paper(s): Capaday et al. 1986, Zehr and Stein 1999, Eccles and Lundberg 1958, Procházka and Prochazka 2011
- **Ia reciprocal inhibition** — 6 paper(s): Kiehn 2016, Jankowska et al. 1967, Rybak 2006, Pratt and Jordan 1987
- **Ib disynaptic excitation** — 6 paper(s): McCrea 1995, Jankowska 1981, Pearson 1998, Nichols 2018
- **Chordotonal organ multi-synaptic excitation** — 5 paper(s): Ekeberg 2004, Hess and Büschges 1999, Hellekes 2012, Ayali 2015
- **Chordotonal organ multi-synaptic inhibition** — 5 paper(s): Ekeberg 2004, Hess and Büschges 1999, Hellekes 2012, Ayali 2015
- **Ib disynaptic inhibition** — 5 paper(s): McCrea 1995, Jankowska 1981, Nichols 2018, Stephens and Yang 1996
- **type II excitatory** — 5 paper(s): Edgley and Jankowska 1987, Donelan and Pearson 2004, Mazzaro 2006, Hatz 2012
- **Ia or II swing to stance** — 4 paper(s): Akay 2014, Lam and Pearson 2002, Lam and Pearson 2001, McVea 2005
- **Ib inhibition** — 4 paper(s): Duysens 2000, Pearson and Collins 1993, Procházka 1997, Nichols and Ross 2009
- **Ib swing to stance** — 3 paper(s): Akay 2014, Ivashko 2003, Domínguez-Rodríguez 2020
- **Ia monosynaptic** — 3 paper(s): Akazawa 1982, Nichols 2018, Ivashko 2003
- **trochanteral campaniform sensilla load signals to adjust MN magnitude** — 3 paper(s): Akay 2004, Ayali 2015, Kaliyamoorthy 2005
- **Ia presynaptic inhibition** — 3 paper(s): Jankowska 1967, Gosgnach 2000, Perreault 1999
- **Ia stance to swing** — 2 paper(s): Verschueren 2003, Ivashko 2003
- **large diameter spinal afferent stimulation** — 2 paper(s): Minassian 2007, Hachmann 2021
- **Ia swing to stance** — 1 paper(s): Ivashko 2003
- **Mechanosensory monosynaptic excitation** — 1 paper(s): Ivashko 2003
- **Ia disynaptic inhibition** — 1 paper(s): Ivashko 2003
- **Somatosensory feedback to postural control** — 1 paper(s): Kooij 2000
- **Vestibular feedback to postural control** — 1 paper(s): Kooij 2000
- **Visual feedback to postural control** — 1 paper(s): Kooij 2000
- **trochanteral hair plate multi-synaptic excitation** — 1 paper(s): Ayali 2015
- **trochanteral hair plate multi-synaptic inhibition** — 1 paper(s): Ayali 2015
- **Ia inhibitory** — 1 paper(s): Hultborn 1971
- **group III/IV fatigue** — 1 paper(s): Amann 2020
- **Total afferent inhibition** — 1 paper(s): Grillner and Zangger 1979
- **Ia disynaptic excitation** — 1 paper(s): Angel 1996
- **I disynaptic excitation (mixed group I)** — 1 paper(s): Angel 2005
- **Type 1 swing to stance** — 1 paper(s): Guertin 1995
- **type 1 stance to swing** — 1 paper(s): Guertin 1995
- **II inhibitory** — 1 paper(s): Edgley and Jankowska 1987
- **Ib contralateral inhibition** — 1 paper(s): Hochman 2013

## 2. Key MODELS (papers that ARE models — 'Is the model paper' links)

NOTE: this field has known multi-link pollution — treat names as
leads, verify against the paper before citing.

- **Schumacher 2025** — Emergence of natural and robust bipedal walking by learning from biologically pl (2025) · 10 cites
- **Shevtsova 2025** — Reorganization of spinal neural connectivity following recovery after thoracic s (2025) · 0 cites
- **Wang 1995** — Emergent synchrony in locally coupled neural oscillators (1995) · 186 cites
- **Kiemel 2008** — Identification of the Plant for Upright Stance in Humans: Multiple Movement Patt (2008) · 112 cites
- **Shriki 2003** — Rate Models for Conductance-Based Cortical Neuronal Networks (2003) · 224 cites
- **Zhu 2024** — A spinal circuit model with asymmetric cervical-lumbar layout controls backward  (2024) · 2 cites
- **Young 2019** — Analyzing Moment Arm Profiles in a Full-Muscle Rat Hindlimb Model (2019) · 16 cites
- **Ekeberg 2004** — Dynamic simulation of insect walking. (2004) · 150 cites
- **Ivashko 2003** — Modeling the spinal cord neural circuitry controlling cat hindlimb movement duri (2003) · 57 cites
- **Owaki 2013** — Simple robot suggests physical interlimb communication is essential for quadrupe (2013) · 189 cites
- **Mo 2025** — A multi-Layered neural control framework: Combining central pattern generators u (2025) · 0 cites
- **Deng 2022** — Biomechanical and Sensory Feedback Regularize the Behavior of Different Locomoto (2022) · 8 cites

## 3. Key STUDIES (top-cited; cite these for mechanisms)

- **LeCun 2015** (84638 cites · [A]) — Deep learning  
  3Department of Computer Science and Operations Research Université de Montréal, Pavillon André-Aisenstadt, PO Box 6128 Centre-Ville STN Montréal, Quebec H3C 3J7.
- **Hodgkin 1952** (23440 cites) — A quantitative description of membrane current and its application to conduction and   
  Title-level grounding only, because the supplied abstract is publisher boilerplate with no scientific content: this is the foundational quantitative description.
- **Hornik 1989** (21727 cites · [A]) — Multilayer feedforward networks are universal approximators  
  Deep Learning: “Multilayer feedforward networks are universal approximators” (Hornik, 1989) 2..
- **Shahriari 2016** (6221 cites) — Taking the Human Out of the Loop: A Review of Bayesian Optimization  
  Review of taking human out loop:.
- **Izhikevich 2003** (4933 cites) — Simple model of spiking neurons  
  Computational model study.
- **Izhikevich 2006** (3878 cites · [A]) — Dynamical Systems in Neuroscience: The Geometry of Excitability and Bursting  
  Standard reference connecting electrophysiology to nonlinear dynamical systems theory: neuronal information processing depends on dynamical as well as electroph.
- **Wilson 1972** (3802 cites) — Excitatory and Inhibitory Interactions in Localized Populations of Model Neurons  
  The Wilson-Cowan coupled nonlinear equations for localized excitatory and inhibitory neuron populations exhibit simple and multiple hysteresis, multiple stable .
- **Horak 2006** (2749 cites) — Postural orientation and equilibrium: what do we need to know about neural control of  
  Frames postural control as two interacting goals - orientation (active alignment of trunk and head to gravity, support surface, and visual surround) and equilib.
- **Tanaka 2019** (2147 cites · [A]) — Recent advances in physical reservoir computing: A review  
  Review of recent advances physical reservoir computing: review.
- **Jessell 2000** (2121 cites · [A]) — Neuronal specification in the spinal cord: inductive signals and transcriptional code  
  Reviews the developmental logic of spinal neuronal specification: inductive signals and transcriptional codes assign class identity to progenitors at defined po.
- **Ijspeert 2008** (1879 cites) — 2008 Special Issue: Central pattern generators for locomotion control in animals and   
  The canonical review of locomotor CPGs spanning neurobiology and robotics: CPG circuits turn simple low-dimensional inputs into coordinated high-dimensional rhy.
- **Ling 2016** (1655 cites · [A]) — Reynolds averaged turbulence modelling using deep neural networks with embedded invar  
  Computational model study.
- **Duraisamy 2019** (1517 cites · [A]) — Turbulence Modeling in the Age of Data  
  Review of turbulence modeling age data.
- **Stevens 2006** (1514 cites) — The costs of fatal and non-fatal falls among older adults  
  Review of costs fatal non-fatal falls among.
- **Grillner 2011** (1376 cites) — Control of Locomotion in Bipeds, Tetrapods, and Fish  
  Comprehensive synthesis of locomotor control across bipeds, tetrapods, and fish: spinal unit-burst-generator organization, interlimb coordination, load- and hip.
- **Maass 1997** (1325 cites) — Networks of spiking neurons: The third generation of neural network models  
  Computational model study.
- **Sherrington 1910** (1325 cites) — Flexion-reflex of the limb, crossed extension-reflex, and reflex stepping and standin  
  Defined the flexion-reflex as a protective type-reflex of the whole limb and its companion crossed extension-reflex in the spinal cat and dog, and showed that r.
- **Marder and Bucher 2001** (1275 cites) — Central pattern generators and the control of rhythmic movements  
  Review of central pattern generators: rhythmic motor output can be produced with no sensory or descending timing cues, and neuromodulation reconfigures the same.
- **Burkitt 2006** (1244 cites) — A Review of the Integrate-and-fire Neuron Model: I. Homogeneous Synaptic Input  
  Review of review integrate-and-fire neuron model: i..
- **Brown 1911** (1221 cites) — The intrinsic factors in the act of progression in the mammal  
  Rhythmic stepping persists in the low-spinal cat after complete deafferentation of the hindlimb muscles, and Graham Brown concluded the cycle is generated by pa.
- **Rossignol 2006** (1115 cites) — Dynamic Sensorimotor Interactions in Locomotion  
  Definitive synthesis of locomotion as a dynamic interaction between a genetically determined spinal CPG and sensory feedback: extensor proprioceptive input adju.
- **Ivanenko 2004** (1102 cites) — Five basic muscle activation patterns account for muscle activity during human locomo  
  We recorded from12–16 ipsilateral leg and trunk muscles
using both surface and intramuscular recording and determined the average, normalized EMG
of each record.
- **Grillner 1985** (1085 cites) — Neurobiological Bases of Rhythmic Motor Acts in Vertebrates  
  Classic formulation of the principles underlying nervous control of innate rhythmic motor acts in vertebrates, establishing the in vitro lamprey spinal cord pre.
- **Taga 1991** (1036 cites) — Self-organized control of bipedal locomotion by neural oscillators in unpredictable e  
  Proposed that stable and flexible locomotion emerges as a global limit cycle produced by global entrainment between coupled neural oscillators, the musculoskele.
- **Capaday et al. 1986** (978 cites) — Amplitude modulation of the soleus H-reflex in the human during walking and standing  
  The soleus H-reflex is strongly modulated over the step cycle during walking versus standing; central gating of the monosynaptic Ia pathway is phase- and task-d.
- **Takakusaki 2017** (968 cites) — Functional Neuroanatomy for Posture and Gait Control  
  Copyright © 2017 The Korean Movement Disorder Society 1 Functional Neuroanatomy for Posture and Gait Control Kaoru Takakusaki The Research Center for Brain Func.
- **Grillner 2006** (940 cites) — Biological Pattern Generation: The Cellular and Computational Logic of Networks in Mo  
  Review articulating the 'circuit doctrine' for motor pattern generation: the mode of operation of central motor-program (CPG) networks, the neural mechanisms by.
- **Jankowska 1992** (933 cites) — INTERNEURONAL RELAY IN SPINAL PATHWAYS FROM PROPRIOCEPTORS  
  review pathway and sensory.
- **Herman 1976** (920 cites) — Neural control of locomotion  
  Landmark edited volume (Advances in Behavioral Biology, v.
- **Kiehn 2006** (906 cites) — LOCOMOTOR CIRCUITS IN THE MAMMALIAN SPINAL CORD  
  Landmark review of locomotor circuits in the mammalian spinal cord: rhythmogenesis, pattern formation, and the genetic dissection of interneuron classes..

## 4. Afferent cheat sheet (encoding conventions used across the stacks)

- **Ia** (spindle primary): length + velocity. MuJoCo-SNS encoder:
  `(L−Lmid)/Lhalf` and `L̇/0.6` → afferent current (spinal/ build_network).
- **Ib** (Golgi tendon organ): force. Encoder: `F/Fmax`.
- **II** (spindle secondary): static length; often stance-gated to avoid
  a global co-contraction floor (the ungated-i0 lesson in DESIGN.md).
- **group III/IV**: metabolic/fatigue-nociceptive — 'group III/IV fatigue' rule.
- **Cutaneous / mechanosensory**: heel/toe contact sensors → HEEL_IN/TOE_IN
  ports; insect hair plates/campaniform/chordotonal = Mechanosensory.
- **Phase → afferent gain** (Akazawa 1982 direction): reflex gain is
  locomotor-phase dependent — the RG/PF state gates every afferent pathway.

## 5. Platform cross-reference (where rules live in each toolchain)

| Rule family | MuJoCo-SNS (spinal/) | SNS Simulink (SNS_Simscape) | AnimatLab (Neuromechanical_Models) |
|---|---|---|---|
| Rhythm generation (RG half-centers, Brown 1911 lineage) | per-leg RG; NaP/tau-h traps in AGENTS | SNS_Library NonSpikingNeuron; SNS_SpinalNetwork | Biped_2xCPG_wSubs RG (contact-driven; W2L) |
| Ia reciprocal inhibition | params ia_in / IaIN population | KneeReflexDemo reciprocal Ia pair | Deng-style Ia-IN chains in aproj |
| Ib autogenic + stance-gated reversal (IBEXC) | full_rules Ib branch; ib_rge | — | W2L Ib autogenic excitation |
| Ib/Ia disynaptic excitation & inhibition (Angel/Jankowska) | full_rules group-I INs | — | Ia relay INs in W2L |
| Heel/toe contact reset at PF layer | heel_pf_layer, toe_df_inh (Ben's rules JSON) | — | contact neurons → RG (W2L reference) |
| Renshaw recurrent inhibition | --renshaw G | — | Renshaw chains in aproj |
| Stance-gated phase reset / PRESET transients | PRESET_E/F high-pass onset | — | — |
| II afferents stance-gated | ii gating in build_network | — | — |
| Knee convention traps | flexion-negative knee (audit_signs) | — | — |

Platform facts above are summaries of AGENTS.md's toolchain sections —
read the full sections there before implementation work.

## 6. Non-negotiables when porting rules across platforms

- Tau-h semantics differ (sns_toolbox tau_h(V) quenches NaP half-centers;
  fixed tau_h works — spinal/_tau_h_check.py, deng_cpg_ode.py).
- Ia/Ib/II afferents must be pure-signal (resting tone → constant
  reciprocal inhibition — the i0 lesson).
- Summation-order chaos flips marginal limit cycles between backends
  (basin gate: basin_gate.py).
- SNS_Library vs sns_toolbox use OPPOSITE synapse-saturation conventions.

