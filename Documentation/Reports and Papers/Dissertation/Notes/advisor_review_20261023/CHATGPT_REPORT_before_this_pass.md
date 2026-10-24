# Advisor review and revision handoff

Updated from the live **Bolen-Dissertation** Overleaf Review → Overview, final capture **22 October 2026, 12:22:09 Pacific**. This was a review of advisor contents and their anchors, not another full-document quality review. The advisor is actively reviewing; this is a snapshot, not a claim that no later comments exist.

## Current inventory

**73 open advisor threads (76 advisor messages) and 51 pending advisor tracked edits.** Ben's replies are explicitly separated below. Resolved comments and Ben-only appendix/main entries are excluded from the advisor totals. The previous 36-thread/8-edit inventory is superseded and preserved in Notes/advisor_review_20261022/CHATGPT_REPORT_before_this_pass.md. Existing SNS-01..15, DB-01..26, INTRO-01, BG-01, and BG-02 identifiers are retained.

| File | Advisor threads | Pending advisor edits |
|---|---:|---:|
| 10-introduction.tex | 14 | 20 |
| 15-background.tex | 24 | 25 |
| 50-database.tex | 24 | 2 |
| 60-sns.tex | 11 | 4 |

Line numbers and context refer to the fresh pre-upload source snapshot in tmp/advisor_review_20261022/overleaf_before/chapters, except the later live background snapshot saved as Notes/advisor_review_20261022/background_final.tex. They do not refer to the differing local ProofFinal chapters or PDF line numbers. A long paragraph occupies one source line. Positions are the rendered review panel's character offsets. Exact highlighted text was captured from the editor where available; a context window is explicitly labeled when no highlight was returned. Thread IDs are retained in the JSON evidence. Source offsets can drift after concurrent edits. Tracked edits remain pending; recording them does not accept them.

## Ben's current controlling instructions

These supersede conflicting older campaign plans and the earlier ZCode prompt:

1. Reduce the dissertation, especially SADb and SNS. Quality and a clear argument take precedence over more text. Do not bury a cut by adding an equally long appendix.
2. SNS does **not** require Methods/Results/Discussion headings. Lead with Ben's PhD-time combination of the Deng CPG and Li walker in AnimatLab: stepping-like movement in air, but no walking achieved by hand tuning. Use that limitation to motivate the knowledge/capability gaps in Ben's replies below. Follow with brief results of new software tools and concrete verification showing what works and what gap is closed. Future work belongs mainly in Chapter 7.
3. Preserve the advisor's active review. Supply narrowly targeted changes wherever possible; prefer a few paragraphs or one section to replacement of entire chapters. Log moves/deletions, retain labels and citations, and do not overwrite Overleaf from stale local chapters.
4. Bring Chapter 4 into line with verified recent torque-paper work from GitHub, distinguishing new measurements, uncorrected estimates, completed reductions, and predicted results.
5. Preserve Chapter 3's isometric-force figures and MATLAB figures. Reports and prompts stay local; **never upload CHATGPT_REPORT or handoff/administrative files to Overleaf.**

## Advisor comments and tracked changes

## 10-introduction.tex

### INTRO-01 · line 8 · comment · 14 October, 10:59 am

**Highlighted source:**

```latex
Artificial
```

**Advisor — 14 October, 10:59 am:**

> This paragraph does a great job identifying one of the gaps of knowledge you worked to fill, but there are a few others you can worked on too.


### INTRO-02 · line 12 · tracked edit · 22 October, 9:38 am

**Context near the review anchor (not an exact highlight):**

```latex
rque data, so even the nominal torque envelope of the robot leg rested on unverified calculations.

There are additional gaps concerning the development of more detailed neural control models of locomotion and balance. No simulation pipeline existed for closing the loop between a musculoskeletal robot model of this kin
```

**Pending advisor edit:**

```text
Changed: The remaining
to There are additional
```


### INTRO-03 · line 12 · tracked edit · 22 October, 9:38 am

**Context near the review anchor (not an exact highlight):**

```latex
rque envelope of the robot leg rested on unverified calculations.

There are additional gaps concerning the development of more detailed neural control models of locomotion and balance. No simulation pipeline existed for closing the loop between a musculoskeletal robot model of this kind and a biologically grounded neu
```

**Pending advisor edit:**

```text
Added: ing
```


### INTRO-04 · line 12 · tracked edit · 22 October, 9:36 am

**Context near the review anchor (not an exact highlight):**

```latex
elope of the robot leg rested on unverified calculations.

There are additional gaps concerning the development of more detailed neural control models of locomotion and balance. No simulation pipeline existed for closing the loop between a musculoskeletal robot model of this kind and a biologically grounded neural cont
```

**Pending advisor edit:**

```text
Changed: controller
to development of more detailed neural control models of locomotion and balance
```


### INTRO-05 · line 12 · tracked edit · 22 October, 9:36 am

**Context near the review anchor (not an exact highlight):**

```latex
 existed for closing the loop between a musculoskeletal robot model of this kind and a biologically grounded neural controller. The neuromechanical modeling literature provides the ingredients, including two-level spinal central pattern generator architectures \citep{rybak_modelling_2006, mccrea_organization_2008}, neu
```

**Pending advisor edit:**

```text
Changed: realistic
to grounded
```


### INTRO-06 · line 12 · tracked edit · 22 October, 9:37 am

**Context near the review anchor (not an exact highlight):**

```latex
}, neuromechanical simulation frameworks \citep{markin_afferent_2010}, and synthetic nervous system design tools \citep{nourse_sns_2023}, but it does not provide their assembly for a BPA-driven bipedal robot with human-like muscle routing. 

Wiring realistic feedback into such a controller exposes a fourth gap, this on
```

**Pending advisor edit:**

```text
Changed: tooling
to design tools
```


### INTRO-07 · line 12 · comment · 22 October, 9:40 am

**Highlighted source:**

```latex
but it does not provide their assembly for a BPA-driven bipedal robot with human-like muscle routing.
```

**Advisor — 22 October, 9:40 am:**

> expand


### INTRO-08 · line 12 · comment · 22 October, 9:37 am

**Highlighted source:**

```latex
their
```

**Advisor — 22 October, 9:37 am:**

> undefined


### INTRO-09 · line 12 · tracked edit · 22 October, 9:41 am

**Context near the review anchor (not an exact highlight):**

```latex
t it does not provide their assembly for a BPA-driven bipedal robot with human-like muscle routing. 

Wiring realistic feedback into such a controller exposes a fourth gap, this one in the literature rather than in methods: the physiological properties of the group Ia, group Ib, and group II afferents that should infor
```

**Pending advisor edit:**

```text
Added: 
```


### INTRO-10 · line 14 · comment · 22 October, 9:42 am

**Context near the review anchor (not an exact highlight):**

```latex
it does not provide their assembly for a BPA-driven bipedal robot with human-like muscle routing. 

Wiring realistic feedback into such a controller exposes a fourth gap, this one in the literature rather than in methods: the physiological properties of the group Ia, group Ib, and group II afferents that should inform 
```

**Advisor — 22 October, 9:42 am:**

> Awkward


### INTRO-11 · line 14 · comment · 22 October, 9:41 am

**Highlighted source:**

```latex
hundreds if not thousands
```

**Advisor — 22 October, 9:41 am:**

> You should have a good idea of how many now that you have gathered them.


### INTRO-12 · line 16 · tracked edit · 22 October, 9:43 am

**Context near the review anchor (not an exact highlight):**

```latex
er or an animal researcher could draw them.

Four gaps thus stood between the existing robot design methods and a robot that can walk under neural control: no actuator force model valid at arbitrary resting length, no quantified correction for the torque lost through brackets, tendons, and wrapping, no closed-loop simu
```

**Pending advisor edit:**

```text
Added: methods
```


### INTRO-13 · line 16 · comment · 22 October, 9:44 am

**Highlighted source:**

```latex
no closed-loop simulation pipeline joining the robot's musculoskeletal model to a biologically realistic neural controller, and no consolidated resource for the afferent data that such a controller requires.
```

**Advisor — 22 October, 9:44 am:**

> use mirror language here and for your chapter order. Your chapters should be in the order you present the gaps/you should present the gaps in the order your chapters are.


### INTRO-14 · line 22 · comment · 22 October, 9:49 am

**Highlighted source:**

```latex
In parallel, a biomimetic four-bar knee linkage with a migrating instantaneous center of rotation was designed and validated by Steele \citep{steele_development_2017, steele_biomimetic_2018, steele_biomimetic_2018-1}, and a canine-inspired quadruped platform, Muscle Mutt, formerly called DoggyDeux, was developed with antagonistic BPA pairs and synthetic nervous system control \citep{scharzenberger_design_2019}. That platform extended an earlier AARL demonstration of central pattern generator control of a pneumatically actuated robot \citep{hunt_development_2017}.
```

**Advisor — 22 October, 9:49 am:**

> It is not clear how this work connects to the Biomimetic Humanoid Project and how your work builds on it.


### INTRO-15 · line 26 · comment · 22 October, 9:49 am

**Highlighted source:**

```latex
first two gaps
```

**Advisor — 22 October, 9:49 am:**

> undef

**Advisor — 22 October, 9:50 am:**

> You should name them specifically


### INTRO-16 · line 26 · comment · 22 October, 9:49 am

**Highlighted source:**

```latex
fourth
```

**Advisor — 22 October, 9:49 am:**

> undef


### INTRO-17 · line 26 · comment · 22 October, 9:50 am

**Highlighted source:**

```latex
third
```

**Advisor — 22 October, 9:50 am:**

> undef


### INTRO-18 · line 41 · tracked edit · 22 October, 9:54 am

**Context near the review anchor (not an exact highlight):**

```latex
are adapted from the author's publication of this work \citep{bolen_2026}.

Chapter~\ref{ch:torque} improves the actuator model at the joint. It describes the robot leg model, the pinned-joint and biomimetic knee test stands, and the equipment used to collect isometric torque data; formulates the hybrid torque calculat
```

**Pending advisor edit:**

```text
Changed: carries
to improves
```


### INTRO-19 · line 41 · tracked edit · 22 October, 9:54 am

**Context near the review anchor (not an exact highlight):**

```latex
's publication of this work \citep{bolen_2026}.

Chapter~\ref{ch:torque} improves the actuator model at the joint. It describes the robot leg model, the pinned-joint and biomimetic knee test stands, and the equipment used to collect isometric torque data; formulates the hybrid torque calculation method and the improved
```

**Pending advisor edit:**

```text
Deleted: s
```


### INTRO-20 · line 41 · tracked edit · 22 October, 9:54 am

**Context near the review anchor (not an exact highlight):**

```latex
s publication of this work \citep{bolen_2026}.

Chapter~\ref{ch:torque} improves the actuator model at the joint. It describes the robot leg model, the pinned-joint and biomimetic knee test stands, and the equipment used to collect isometric torque data; formulates the hybrid torque calculation method and the improved 
```

**Pending advisor edit:**

```text
Changed: to
to at
```


### INTRO-21 · line 41 · comment · 22 October, 9:55 am

**Highlighted source:**

```latex
Gait2392
```

**Advisor — 22 October, 9:55 am:**

> Too specific


### INTRO-22 · line 45 · tracked edit · 22 October, 9:57 am

**Context near the review anchor (not an exact highlight):**

```latex
on living subjects. It opens with the preliminary AnimatLab bipedal walker results, whose hand-tuned, two-layer CPG produced suspended stepping, with neither weight-bearing walking nor any role for afferent feedback in weight-bearing motion demonstrated; that outcome set the verification-first approach the developed to
```

**Pending advisor edit:**

```text
Added: ,
```


### INTRO-23 · line 45 · comment · 22 October, 9:58 am

**Highlighted source:**

```latex
nor any role for afferent feedback in weight-bearing motion demonstrated
```

**Advisor — 22 October, 9:58 am:**

> huh?


### INTRO-24 · line 45 · tracked edit · 22 October, 9:59 am

**Context near the review anchor (not an exact highlight):**

```latex
eedback in weight-bearing motion demonstrated; that outcome set the verification-first approach the developed tools use. The tools themselves are the OpenSim-to-MuJoCo conversion of the robot's musculoskeletal model with its synthetic BPA knee, the Simulink SNS library built to prototype reflex topologies, the spiking-
```

**Pending advisor edit:**

```text
Added: developed
```


### INTRO-25 · line 45 · tracked edit · 22 October, 9:59 am

**Context near the review anchor (not an exact highlight):**

```latex
t-bearing motion demonstrated; that outcome set the verification-first approach the developed tools use. The tools themselves are the OpenSim-to-MuJoCo conversion of the robot's musculoskeletal model with its synthetic BPA knee, the Simulink SNS library built to prototype reflex topologies, the spiking-neuron mirrors, 
```

**Pending advisor edit:**

```text
Changed: follow
to use
```


### INTRO-26 · line 45 · comment · 22 October, 10:00 am

**Highlighted source:**

```latex
the spiking-neuron mirrors
```

**Advisor — 22 October, 10:00 am:**

> huh?


### INTRO-27 · line 45 · tracked edit · 22 October, 10:02 am

**Context near the review anchor (not an exact highlight):**

```latex
, and the cross-platform verification protocol that ties the implementations together. Each tool is validated before the complete model is assembled. The stimulus and deletion experiments complete the validation process, and the chapter closes by discussing a five-stage training curriculum, the spiking-neuron mirrors, 
```

**Pending advisor edit:**

```text
Changed: demonstrated on the fundamental exercises that establish what it does and does not show: cross-engine reproduction of the contact-driven rhythm core, headless re-runs of the archived AnimatLab reference models, the BPA-actuated walker in MuJoCo, a gait-validation sweep against a locomotion reference library, and a comparison of the production walker variants
to validated before the complete model is assembled
```


### INTRO-28 · line 45 · tracked edit · 22 October, 10:02 am

**Context near the review anchor (not an exact highlight):**

```latex
alidated before the complete model is assembled. The stimulus and deletion experiments complete the validation process, and the chapter closes by discussing a five-stage training curriculum, the spiking-neuron mirrors, and the new tools the program leaves other researchers.

Chapter~\ref{ch:overall} discusses the work 
```

**Pending advisor edit:**

```text
Changed: demonstrations
to validation process
```


### INTRO-29 · line 45 · tracked edit · 22 October, 10:03 am

**Context near the review anchor (not an exact highlight):**

```latex
ulus and deletion experiments complete the validation process, and the chapter closes by discussing a five-stage training curriculum, the spiking-neuron mirrors, and the new tools the program leaves other researchers.

Chapter~\ref{ch:overall} discusses the work as a whole. It shows how the four chapters work together 
```

**Pending advisor edit:**

```text
Changed: the
to a
```


### INTRO-30 · line 45 · comment · 22 October, 10:03 am

**Highlighted source:**

```latex
spiking-neuron mirrors
```

**Advisor — 22 October, 10:03 am:**

> again, huh?


### INTRO-31 · line 47 · tracked edit · 22 October, 10:04 am

**Context near the review anchor (not an exact highlight):**

```latex
s the program leaves other researchers.

Chapter~\ref{ch:overall} discusses the work as a whole. It shows how the four chapters work together to develop improved biomimetic humanoid robots. Its also scopes the future work these results enable for other researchers, from BPA characterization and routing expansion throug
```

**Pending advisor edit:**

```text
Changed: opens with the outcome, that improved first-pass robot design is now achievable, and attributes it: the design factors that remove the rigid-body idealization, the improved BPA characterization, the tiered attachment and routing optimization, the sensory afferent database, and the new neuromechanical simulation tools
to shows how the four chapters work together to develop improved biomimetic humanoid robots
```


### INTRO-32 · line 47 · tracked edit · 22 October, 10:08 am

**Context near the review anchor (not an exact highlight):**

```latex
e. It shows how the four chapters work together to develop improved biomimetic humanoid robots. Its also scopes the future work these results enable for other researchers, from BPA characterization and routing expansion through database and simulation tools to real-time neural control of the physical robot. Chapter~\re
```

**Pending advisor edit:**

```text
Changed: remaining subsections
to also
```


### INTRO-33 · line 47 · tracked edit · 22 October, 10:08 am

**Context near the review anchor (not an exact highlight):**

```latex
s how the four chapters work together to develop improved biomimetic humanoid robots. Its also scopes the future work these results enable for other researchers, from BPA characterization and routing expansion through database and simulation tools to real-time neural control of the physical robot. Chapter~\ref{ch:concl
```

**Pending advisor edit:**

```text
Added: s
```


### INTRO-34 · line 47 · tracked edit · 22 October, 10:08 am

**Context near the review anchor (not an exact highlight):**

```latex
 other researchers, from BPA characterization and routing expansion through database and simulation tools to real-time neural control of the physical robot. Chapter~\ref{ch:conclusion} consolidates the principal contributions of the dissertation. Figure~\ref{fig:roadmap} summarizes the flow.

\begin{figure}[htbp]
\cent
```

**Pending advisor edit:**

```text
Changed: tooling
to tools
```


## 15-background.tex

### BG-03 · line 36 · comment · 22 October, 10:42 am

**Highlighted source:**

```latex
\begin{figure}[htbp]
```

**Advisor — 22 October, 10:42 am:**

> Motor center is a little off


### BG-01 · line 48 · tracked edit · 14 October, 1:31 pm

**Context near the review anchor (not an exact highlight):**

```latex
th relationship. Predicting their force as a function of length and pressure is therefore the first improvement of the torque tools developed in this dissertation, and the state of that prediction is the subject of the next subsection.

Figure~\ref{fig:festo_bpa_cad} shows the commercial actuator assembly used in the l
```

**Pending advisor edit:**

```text
Changed: requirement
to improvement
```


### BG-04 · line 100 · comment · 22 October, 10:45 am

**Highlighted source:**

```latex
\begin{figure}[htbp]
```

**Advisor — 22 October, 10:45 am:**

> This figure is confusing and probably not needed,


### BG-05 · line 111 · comment · 22 October, 10:53 am

**Highlighted source:**

```latex
Gait2392
```

**Advisor — 22 October, 10:53 am:**

> I assume the 23 corresponds to the DoF. What about the 92?


### BG-06 · line 113 · comment · 22 October, 10:49 am

**Highlighted source:**

```latex
the full-muscle rat hindlimb model of \citet{young_analyzing_2019} makes that adjacent-joint dependence a first-class result of its attachment-point technique.
```

**Advisor — 22 October, 10:49 am:**

> awkward


### BG-07 · line 113 · comment · 22 October, 10:49 am

**Highlighted source:**

```latex
first-class
```

**Advisor — 22 October, 10:49 am:**

> primary?


### BG-31 · line 113 · tracked edit · 22 October, 10:49 am

**Context near the review anchor (not an exact highlight):**

```latex
-point technique. OpenSim evaluates moment arms as muscle-length sensitivities of this kind, so the moment arms used throughout this dissertation are properties of the Gait2392 geometry rather than measured inputs \citep{delp_opensim_2007, seth_opensim_2018}.

The human knee is the joint where these definitions are str
```

**Pending advisor edit:**

```text
Added: moment
```


### BG-09 · line 113 · comment · 22 October, 10:50 am

**Highlighted source:**

```latex
used throughout this dissertation
```

**Advisor — 22 October, 10:50 am:**

> not on your collected data though, right? Be more specific.


### BG-10 · line 115 · comment · 22 October, 10:50 am

**Highlighted source:**

```latex
stressed
```

**Advisor — 22 October, 10:50 am:**

> huh?

**Advisor — 22 October, 10:50 am:**

> definitions are stressed?


### BG-11 · line 115 · comment · 22 October, 10:51 am

**Highlighted source:**

```latex
 Its instantaneous center of rotation migrates posteriorly with flexion, a behavior captured by four-bar linkage descriptions of the joint \citep{romero_four-bar_2020, steele_biomimetic_2018}, and \citet{romero_four-bar_2020} treat locating that moving center as a kinematic problem in its own right, solved with natural coordinates so that the answer stays cheap to evaluate when an optimization must recompute it many times
```

**Advisor — 22 October, 10:51 am:**

> run on. awkward


### BG-12 · line 115 · tracked edit · 22 October, 10:52 am

**Context near the review anchor (not an exact highlight):**

```latex
femoral joint give different moment arms for identical muscle paths, so the modeling choice must be determined carefully. The ISB joint coordinate system recommendations were written against exactly this definitional freedom, to keep reports of joint motion comparable across laboratories \citep{wu_isb_2002}. Throughout
```

**Pending advisor edit:**

```text
Changed: matched to the question asked of it
to determined carefully
```


### BG-13 · line 115 · comment · 22 October, 10:52 am

**Highlighted source:**

```latex
written against
```

**Advisor — 22 October, 10:52 am:**

> huh?


### BG-14 · line 115 · comment · 22 October, 10:53 am

**Highlighted source:**

```latex
judged
```

**Advisor — 22 October, 10:53 am:**

> I'm now a little lost...which method does 2392 use?


### BG-15 · line 119 · comment · 22 October, 10:55 am

**Highlighted source:**

```latex
It is a one-degree-of-freedom four-bar linkage whose crossed links reproduce the posterior migration of the human knee's instantaneous center of rotation (ICR) during flexion, so the rolling-and-sliding behavior of the human joint emerges from the link geometry itself rather than from cams or gear profiles; the design was calculated to reproduce an average sliding-to-rolling ratio of 0.079 over a \ang{151} range of motion \citep{steele_biomimetic_2018-1}.
```

**Advisor — 22 October, 10:55 am:**

> run on(s)


### BG-16 · line 119 · comment · 22 October, 10:56 am

**Highlighted source:**

```latex
which is exactly the behavior a torque prediction tool must capture
```

**Advisor — 22 October, 10:56 am:**

> non sequitor


### BG-17 · line 131 · tracked edit · 22 October, 10:58 am

**Context near the review anchor (not an exact highlight):**

```latex
mination_2019}. Hand placement of muscle origins and insertions was a painstaking and iterative task. This resulted in the creation of a placement optimization algorithm \citep{morrow_optimization_2020}, discussed in Chapter~\ref{ch:torque}. Figure~\ref{fig:aarl_opensim_comparison} shows the laboratory's bipedal skelet
```

**Pending advisor edit:**

```text
Changed: , and the
to . This resulted in the creation of a
```


### BG-18 · line 131 · tracked edit · 22 October, 10:58 am

**Context near the review anchor (not an exact highlight):**

```latex
ions was a painstaking and iterative task. This resulted in the creation of a placement optimization algorithm \citep{morrow_optimization_2020}, discussed in Chapter~\ref{ch:torque}. Figure~\ref{fig:aarl_opensim_comparison} shows the laboratory's bipedal skeleton beside the OpenSim muscle-routing reference developed fo
```

**Pending advisor edit:**

```text
Changed: that replaced it is dissertation work
to algorithm
```


### BG-19 · line 131 · comment · 22 October, 11:49 am

**Highlighted source:**

```latex
igure~\ref{fig:aarl_opensim_comparison} shows the laboratory's bipedal skeleton beside the OpenSim muscle-routing reference developed for it
```

**Advisor — 22 October, 11:49 am:**

> non seq


### BG-20 · line 131 · tracked edit · 22 October, 10:59 am

**Context near the review anchor (not an exact highlight):**

```latex
_optimization_2020}, discussed in Chapter~\ref{ch:torque}. Figure~\ref{fig:aarl_opensim_comparison} shows the laboratory's bipedal skeleton beside the OpenSim muscle-routing reference developed for it.

\begin{figure}[htbp]
\centering
\begin{minipage}[c]{0.67\textwidth}
\centering
\textbf{(A)}\par\smallskip
\includegra
```

**Pending advisor edit:**

```text
Changed: places
to shows
```


### BG-21 · line 152 · tracked edit · 22 October, 11:50 am

**Context near the review anchor (not an exact highlight):**

```latex
n_determination_2019, morrow_optimization_2020}.}
\label{fig:aarl_opensim_comparison}
\end{figure}

An overview of the mechanics of that master's study are described below. Moment arms were computed geometrically from attachment points in the body frame, as the cross product of a joint-axis unit vector with the vector 
```

**Pending advisor edit:**

```text
Changed: T
to An overview of t
```


### BG-22 · line 152 · tracked edit · 22 October, 11:50 am

**Context near the review anchor (not an exact highlight):**

```latex
{fig:aarl_opensim_comparison}
\end{figure}

An overview of the mechanics of that master's study are described below. Moment arms were computed geometrically from attachment points in the body frame, as the cross product of a joint-axis unit vector with the vector from the joint center to the muscle's line of action, as
```

**Pending advisor edit:**

```text
Changed: worth recording because the dissertation's torque work reuses them
to described below
```


### BG-23 · line 152 · tracked edit · 22 October, 11:51 am

**Context near the review anchor (not an exact highlight):**

```latex
t of a joint-axis unit vector with the vector from the joint center to the muscle's line of action, as described in Section~\ref{sec:bg_biomechanics}. PAM routings were then placed by hand in OpenSim beside the Gait2392 paths they imitate, moved to respect Festo BPA contraction limits and to keep actuators from interfe
```

**Pending advisor edit:**

```text
Changed: the same construction
to as
```


### BG-24 · line 152 · comment · 22 October, 11:56 am

**Highlighted source:**

```latex
PAM routings were then placed by hand in OpenSim beside the Gait2392 paths they imitate, moved to respect Festo BPA contraction limits and to keep actuators from interfering with one another, and the result was saved as the gait2392\_robotbody model that the later chapters of this dissertation continue to refine \citep{bolen_bipedal_2019}. 
```

**Advisor — 22 October, 11:56 am:**

> Delete? 


### BG-25 · line 152 · comment · 22 October, 11:51 am

**Highlighted source:**

```latex
routings
```

**Advisor — 22 October, 11:51 am:**

> origin and insertion points?


### BG-26 · line 152 · comment · 22 October, 11:51 am

**Highlighted source:**

```latex
 moved to respect Festo BPA contraction limits
```

**Advisor — 22 October, 11:51 am:**

> huh?


### BG-27 · line 152 · comment · 22 October, 11:52 am

**Highlighted source:**

```latex
and the result was saved
```

**Advisor — 22 October, 11:52 am:**

> I don't need to be told that you saved a file.


### BG-28 · line 152 · comment · 22 October, 11:56 am

**Highlighted source:**

```latex
PAM routings were then placed by hand in OpenSim beside the Gait2392 paths they imitate, moved to respect Festo BPA contraction limits and to keep actuators from interfering with one another, and the result was saved as the gait2392\_robotbody model that the later chapters of this dissertation continue to refine \citep{bolen_bipedal_2019}. 
```

**Advisor — 22 October, 11:56 am:**

> I feel like this is jumping all over the place


### BG-29 · line 152 · comment · 22 October, 11:59 am

**Highlighted source:**

```latex
followed
```

**Advisor — 22 October, 11:59 am:**

> is calculated?


### BG-30 · line 152 · comment · 22 October, 11:57 am

**Highlighted source:**

```latex
The biceps femoris long head shows what the comparison exposed:
```

**Advisor — 22 October, 11:57 am:**

> what comparison was made? You don't say


### BG-32 · line 182 · tracked edit · 22 October, 12:05 pm

**Context near the review anchor (not an exact highlight):**

```latex
tor controller developed in Chapter~\ref{ch:sns} is built directly on this two-level organization.


\section{Sensory Feedback in Locomotion}\label{sec:bg_sensory}

Sensory feedback is not a passive correction layer but an integral part of locomotor pattern generation. Reviews of afferent control of the CPG \citep{van_
```

**Pending advisor edit:**

```text
Deleted: \begin{sidewaysfigure}[p] \centering \includegraphics[width=0.97\textheight,height=0.82\textwidth,keepaspectratio]{circuit_literature.pdf} % ALT TEXT: A two-panel synthesis shows a bilateral spinal locomotor network % and a detailed knee-control motif with rhythm-generator, pattern-formation, % motoneuron, commissural, and proprioceptive pathways. \caption[Locomotor circuit architecture synthesized from the literature]{Locomotor circuit architecture synthesized from the literature. The bilateral organization combines rhythm-generation, pattern-formation, motoneuron, commissural, and proprioceptive pathways used across two-level CPG and neuromechanical models \citep{rybak_modelling_2006, rybak_organization_2015, shevtsova_modeling_2016, deng_neuromechanical_2019}. The detailed knee-control motif makes the reciprocal antagonist and sensory pathways explicit.} \label{fig:locomotor_circuit_synthesis} \end{sidewaysfigure}
```


### BG-33 · line 185 · tracked edit · 22 October, 12:04 pm

**Context near the review anchor (not an exact highlight):**

```latex
stems to locomotion are the subject of ongoing study \citep{akay_relative_2021, akay_sensory_2020}. Tonic central and sensory stimuli have been shown to facilitate involuntary air-stepping in humans \citep{selionov_tonic_2009}, demonstrating that simple graded stimulation of motoneuronal and rhythmic circuits can elici
```

**Pending advisor edit:**

```text
Changed: Of particular relevance to the tests planned in Chapter~\ref{ch:futurework}, t
to T
```


### BG-34 · line 185 · tracked edit · 22 October, 12:07 pm

**Context near the review anchor (not an exact highlight):**

```latex
ve and contact feedback channels designed in Chapter~\ref{ch:sns} follow this afferent organization.

\begin{figure}[htbp]
\centering
\includegraphics[width=0.55\textwidth]{rybak_cpg_circuit.pdf}
% ALT TEXT: A three-tier schematic shows the rhythm-generator half-centers
% above a pattern-formation network of interneuro
```

**Pending advisor edit:**

```text
Deleted:  Three published circuits, reproduced in Figures~\ref{fig:deng_connectome}, \ref{fig:afferent_connectomes}, and~\ref{fig:shinohara_connectome}, place these afferent channels in their circuit context. The rat hindlimb controller of \citet{deng_neuromechanical_2019} is organized exactly as Section~\ref{sec:bg_neural} describes: a rhythm-generator half-center and a pattern-formation layer sit above a sensory--motor network in which extensor and flexor muscle receptors feed back through interneurons and Renshaw cells onto the same motoneurons the CPG drives. The human-scale model of \citet{dirusso_investigating_2023} carries the idea to a full lower limb: each agonist--antagonist muscle pair carries Ia, II, and Ib reflex populations with Renshaw cells, and pattern-formation nodes from the locomotor controller terminate on the same motoneuron pools. The cat interlimb model of \citet{shinohara_mechanisms_2025} adds the left--right half of the problem, coupling the two rhythm generators through V3 and commissural interneuron populations while force, length, and velocity feedback returns from every muscle. The schematics differ in animal, scale, and emphasis, but they wire the same afferent classes of this section into the circuit at the level each hypothesis concerns, whether that is the rhythm generator, the pattern-formation layer, or the motor pool. \begin{figure}[htbp] \centering \includegraphics[width=1.0\textwidth]{deng_cpg_connectome.pdf} % ALT TEXT: A four-panel circuit schematic: single-limb control hierarchy, % rhythm-generator configuration, pattern-formation network, and a % sensory-motor network with extensor and flexor receptors, interneurons, % Renshaw cells, motoneurons, and muscles. \caption[Deng two-layer CPG and sensory--motor network]{Synthetic nervous system of the rat hindlimb model: \textbf{(A)} single limb control hierarchy, \textbf{(B)} rhythm-generator configuration, \textbf{(C)} pattern-formation network, and \textbf{(D)} sensory--motor network. Reproduced from \citet[Fig.~4]{deng_neuromechanical_2019} under CC BY 4.0.} \label{fig:deng_connectome} \end{figure} \begin{figure}[htbp] \centering \includegraphics[width=0.62\textwidth]{dirusso_spinal_network.pdf} % ALT TEXT: A spinal network between a muscle and its antagonist with pattern % nodes, interneurons, sensory neurons, motoneurons, and two muscles. \caption[Di Russo spinal network with afferent feedback]{Spinal network between a muscle and its antagonist in the human locomotion model of \citet{dirusso_investigating_2023}: Ia, II, and Ib reflex populations with Renshaw cells, and pattern input from the locomotor controller onto the same motoneuron pools. Reproduced from \citet[Fig.~5]{dirusso_investigating_2023} under CC BY 4.0.} \label{fig:afferent_connectomes} \end{figure} \begin{figure}[htbp] \centering \includegraphics[width=0.62\textwidth]{shinohara_connectome.pdf} % ALT TEXT: A bilateral neuromusculoskeletal cat model: rhythm-generator and % pattern-formation layers, commissural and V3 interneurons, motoneuron % pools, and a cat hindlimb on a treadmill with labeled muscles. \caption[Shinohara bilateral neuromusculoskeletal model]{Neuromusculoskeletal cat model of \citet{shinohara_mechanisms_2025}: bilateral rhythm-generator and pattern-formation layers coupled by V3 and commissural interneurons, motoneuron pools for seven muscles per hindlimb, and force, length, and velocity afferent feedback. Reproduced from \citet[Fig.~1]{shinohara_mechanisms_2025} under CC BY 4.0.} \label{fig:shinohara_connectome} \end{figure} The two-level organization these circuits share descends from the modeling program of the Rybak group, and its foundational schematics are reproduced in Figures~\ref{fig:rybak_cpg_circuit} and~\ref{fig:shevtsova_spinal_circuit}. In the concept model of \citet{rybak_modelling_2006}, afferent feedback ascends to both the rhythm generator and the pattern-formation network, so a sensory perturbation can reset the rhythm at the generator or delete motoneuron activity at the pattern level without necessarily doing both, precisely the deletion phenomenology that motivated the two-level split. The extended model of \citet{shevtsova_modeling_2016} carries the same organization to the two-joint muscles, embedding dedicated pattern-formation populations for the bifunctional posterior biceps--semitendinosus and rectus femoris pools, and that biarticular convergence problem is the one a robot muscle routed across two joints presents to the controllers developed in this dissertation.
```


### BG-35 · line 210 · tracked edit · 22 October, 12:09 pm

**Context near the review anchor (not an exact highlight):**

```latex


Computational neuroscience expresses hypotheses about locomotor circuitry as executable models. A computational model is a system of coupled equations in which elements of the schematic correspond to identified neural populations or pathways. The founding hypothesis is Brown's half-center: alternating flexor and exte
```

**Pending advisor edit:**

```text
Changed: locomotor
to computational
```


### BG-36 · line 210 · tracked edit · 22 October, 12:08 pm

**Context near the review anchor (not an exact highlight):**

```latex
otheses about locomotor circuitry as executable models. A computational model is a system of coupled equations in which elements of the schematic correspond to identified neural populations or pathways. The founding hypothesis is Brown's half-center: alternating flexor and extensor activity arises from two reciprocally
```

**Pending advisor edit:**

```text
Deleted:  differential
```


### BG-37 · line 210 · tracked edit · 22 October, 12:09 pm

**Context near the review anchor (not an exact highlight):**

```latex
tor circuitry as executable models. A computational model is a system of coupled equations in which elements of the schematic correspond to identified neural populations or pathways. The founding hypothesis is Brown's half-center: alternating flexor and extensor activity arises from two reciprocally inhibitory central 
```

**Pending advisor edit:**

```text
Deleted: every
```


### BG-38 · line 210 · tracked edit · 22 October, 12:09 pm

**Context near the review anchor (not an exact highlight):**

```latex
cuitry as executable models. A computational model is a system of coupled equations in which elements of the schematic correspond to identified neural populations or pathways. The founding hypothesis is Brown's half-center: alternating flexor and extensor activity arises from two reciprocally inhibitory central populat
```

**Pending advisor edit:**

```text
Added: s
```


### BG-39 · line 210 · tracked edit · 22 October, 12:09 pm

**Context near the review anchor (not an exact highlight):**

```latex
A computational model is a system of coupled equations in which elements of the schematic correspond to identified neural populations or pathways. The founding hypothesis is Brown's half-center: alternating flexor and extensor activity arises from two reciprocally inhibitory central populations with opposing bias \cite
```

**Pending advisor edit:**

```text
Deleted: s
```


### BG-40 · line 210 · tracked edit · 22 October, 12:09 pm

**Context near the review anchor (not an exact highlight):**

```latex
omputational model is a system of coupled equations in which elements of the schematic correspond to identified neural populations or pathways. The founding hypothesis is Brown's half-center: alternating flexor and extensor activity arises from two reciprocally inhibitory central populations with opposing bias \citep{b
```

**Pending advisor edit:**

```text
Deleted:  an
```


### BG-41 · line 210 · tracked edit · 22 October, 12:09 pm

**Context near the review anchor (not an exact highlight):**

```latex
m of coupled equations in which elements of the schematic correspond to identified neural populations or pathways. The founding hypothesis is Brown's half-center: alternating flexor and extensor activity arises from two reciprocally inhibitory central populations with opposing bias \citep{brown_intrinsic_1911}. Modern 
```

**Pending advisor edit:**

```text
Added: s
```


### BG-42 · line 210 · tracked edit · 22 October, 12:09 pm

**Context near the review anchor (not an exact highlight):**

```latex
 equations in which elements of the schematic correspond to identified neural populations or pathways. The founding hypothesis is Brown's half-center: alternating flexor and extensor activity arises from two reciprocally inhibitory central populations with opposing bias \citep{brown_intrinsic_1911}. Modern instantiatio
```

**Pending advisor edit:**

```text
Changed: ; the wiring diagram is not an illustration of the model, it is the model
to s
```


### BG-43 · line 210 · comment · 22 October, 12:10 pm

**Highlighted source:**

```latex
The founding hypothesis is Brown's half-center: alternating flexor and extensor activity arises from two reciprocally inhibitory central populations with opposing bias \citep{brown_intrinsic_1911}.
```

**Advisor — 22 October, 12:10 pm:**

> non seq - this is specifically about CPGs


### BG-44 · line 210 · comment · 22 October, 12:10 pm

**Highlighted source:**

```latex
Modern instantiations of the half-center endow each neuron with persistent sodium channels, whose slow inactivation terminates each burst so that reciprocal inhibition switches the populations without any phasic input \citep{deng_neuromechanical_2019}. Rhythmogenesis, in this tradition, is a property of the intrinsic membrane dynamics, not of a clock elsewhere in the network.
```

**Advisor — 22 October, 12:10 pm:**

> More stuff specific to CPGs, not Comp Neuro


### BG-45 · line 212 · comment · 22 October, 12:11 pm

**Highlighted source:**

```latex
The two-level organization described in Section~\ref{sec:bg_neural} is the field's current working architecture, and its models account for specific electrophysiological observations: the fictive-locomotion deletions reproduced by \citet{rybak_modelling_2006}, the genetically identified interneuron classes mapped onto the model by \citet{rybak_organization_2015}, and the organization of circuits controlling two-joint muscles reviewed by \citet{shevtsova_modeling_2016}. Coupling such circuits to musculoskeletal models shows how afferent feedback shapes CPG output \citep{markin_afferent_2010} and how biomechanics and sensory feedback regularize the behavior of different CPG architectures \citep{deng_biomechanical_2022}. Physiologically plausible models of human locomotion use the same machinery to dissect the relative contributions of reflexes and CPGs, with locomotor primitives, Ia, II, and Ib pathways, and Renshaw cells represented explicitly \citep{dirusso_investigating_2023}, and recent models address interlimb coordination under abrupt perturbations such as the sudden loss of ground support \citep{shinohara_mechanisms_2025}.
```

**Advisor — 22 October, 12:11 pm:**

> Seems like it should be in the CPG section


### BG-46 · line 214 · tracked edit · 22 October, 12:12 pm

**Context near the review anchor (not an exact highlight):**

```latex
 abrupt perturbations such as the sudden loss of ground support \citep{shinohara_mechanisms_2025}.

One way of interpreting computational neuroscience diagrams is from a controls tradition in which the block diagram is the controller itself. Each block is an identifiable neural population with conductance-based dynamic
```

**Pending advisor edit:**

```text
Changed: For a mechanical engineer, the productive way to read this literature is as
to One way of interpreting computational neuroscience diagrams is from
```


### BG-47 · line 214 · tracked edit · 22 October, 12:12 pm

**Context near the review anchor (not an exact highlight):**

```latex
uroscience diagrams is from a controls tradition in which the block diagram is the controller itself. Each block is an identifiable neural population with conductance-based dynamics, and each connection is a synaptic conductance with a physiologically interpretable weight. The controller developed in Chapter~\ref{ch:sn
```

**Pending advisor edit:**

```text
Changed: :
to .
```


### BG-48 · line 214 · tracked edit · 22 October, 12:12 pm

**Context near the review anchor (not an exact highlight):**

```latex
oscience diagrams is from a controls tradition in which the block diagram is the controller itself. Each block is an identifiable neural population with conductance-based dynamics, and each connection is a synaptic conductance with a physiologically interpretable weight. The controller developed in Chapter~\ref{ch:sns}
```

**Pending advisor edit:**

```text
Changed: e
to E
```


### BG-49 · line 214 · tracked edit · 22 October, 12:13 pm

**Context near the review anchor (not an exact highlight):**

```latex
dynamics, and each connection is a synaptic conductance with a physiologically interpretable weight. The controller developed in Chapter~\ref{ch:sns} comes from this tradition: its two-layer organization follows the rat hindlimb model of \citet{deng_neuromechanical_2019}, and its verification tests are drawn from the d
```

**Pending advisor edit:**

```text
Deleted:  The connection schematics published with these models are therefore worth reading as carefully as their equations.
```


### BG-50 · line 214 · comment · 22 October, 12:13 pm

**Context near the review anchor (not an exact highlight):**

```latex
dynamics, and each connection is a synaptic conductance with a physiologically interpretable weight. The controller developed in Chapter~\ref{ch:sns} comes from this tradition: its two-layer organization follows the rat hindlimb model of \citet{deng_neuromechanical_2019}, and its verification tests are drawn from the d
```

**Advisor — 22 October, 12:13 pm:**

> non seq


### BG-02 · line 260 · tracked edit · 14 October, 2:05 pm

**Context near the review anchor (not an exact highlight):**

```latex
cientific case are the same case: the robot is the instrument that makes these experiments possible.

```

**Pending advisor edit:**

```text
Deleted: these tools were
```


## 50-database.tex

### DB-01 · line 6 · comment · 21 October, 11:13 am

**Highlighted source:**

```latex
and design needs parameters, not only principles
```

**Advisor — 21 October, 11:13 am:**

> huh?


### DB-02 · line 6 · comment · 21 October, 11:13 am

**Highlighted source:**

```latex
For each sensory channel, the designer has to know what the receptor encodes, over what ranges of length, velocity, and force it operates and with what receptor dynamics,
```

**Advisor — 21 October, 11:13 am:**

> not sure what this has to do with the database


### DB-03 · line 6 · comment · 21 October, 11:14 am

**Highlighted source:**

```latex
whether the target is ipsilateral or contralateral
```

**Advisor — 21 October, 11:14 am:**

> there are more options


### DB-04 · line 6 · comment · 21 October, 11:14 am

**Highlighted source:**

```latex
Every element of that specification was established experimentally
```

**Advisor — 21 October, 11:14 am:**

> There are still lots of unknowns


### DB-05 · line 6 · comment · 21 October, 11:16 am

**Highlighted source:**

```latex
literature
```

**Advisor — 21 October, 11:16 am:**

> After this sentence, define the gap in knowledge and what you did to fill it.


### DB-06 · line 6 · tracked edit · 21 October, 11:14 am

**Context near the review anchor (not an exact highlight):**

```latex
se, rat, and human preparations, and the results are distributed over an enormous primary literature. A controller designer is not the only reader stranded by that scatter. An animal researcher planning the next experiment wants to know which preparations established a pathway property, and how the property compares ac
```

**Pending advisor edit:**

```text
Deleted:  rather than collected anywhere
```


### DB-07 · line 6 · comment · 21 October, 11:15 am

**Highlighted source:**

```latex
A controller designer is not the only reader stranded by that scatter. An animal researcher planning the next experiment wants to know which preparations established a pathway property, and how the property compares across species, before committing a preparation to it. A roboticist wants to know which receptor encoding a sensor should imitate. Those readers ask different questions of the same literature, and the literature answers none of them without being gathered, structured, and searchable first.
```

**Advisor — 21 October, 11:15 am:**

> This should probably be the beginning of the paragraph


### DB-08 · line 8 · comment · 21 October, 11:19 am

**Highlighted source:**

```latex
Reviews
```

**Advisor — 21 October, 11:19 am:**

> The tone of this paragraph is as an opening paragraph, and yet it is the 2nd. Probably doesn't need to be here. Discussion of the reviews as a starting point probably goes in your 'methods' section of this chapter


### DB-09 · line 10 · comment · 21 October, 11:20 am

**Highlighted source:**

```latex
Neuroscience 2023 \citep{bolen_sensory_2023}
```

**Advisor — 21 October, 11:20 am:**

> I thought it was presented 2x, once without the AI help, and the other with.


### DB-10 · line 10 · comment · 21 October, 11:20 am

**Highlighted source:**

```latex
distilled mechanism notes
```

**Advisor — 21 October, 11:20 am:**

> huh?


### DB-11 · line 10 · comment · 21 October, 11:21 am

**Highlighted source:**

```latex
first-class record
```

**Advisor — 21 October, 11:21 am:**

> huh?


### DB-12 · line 10 · comment · 21 October, 11:21 am

**Highlighted source:**

```latex
pathway records
```

**Advisor — 21 October, 11:21 am:**

> huh?


### DB-13 · line 10 · comment · 21 October, 11:25 am

**Highlighted source:**

```latex
The database lives in a five-table Airtable base (Papers, Feedback, Models, Review Papers, and Recorders)
```

**Advisor — 21 October, 11:25 am:**

> I thought you migrated it to a different system


### DB-14 · line 10 · comment · 21 October, 11:25 am

**Highlighted source:**

```latex
canonical
```

**Advisor — 21 October, 11:25 am:**

> database?


### DB-15 · line 12 · comment · 21 October, 11:26 am

**Highlighted source:**

```latex
This chapter and the controller chapter are therefore built as one program.
```

**Advisor — 21 October, 11:26 am:**

> huh?


### DB-16 · line 12 · comment · 21 October, 11:28 am

**Highlighted source:**

```latex
The database informs the sensory channels of Chapter~\ref{ch:sns}: the afferent classes the controller monitors, the signal encoding each channel carries, and the connection rules its circuit implements all trace to records held here, and Section~\ref{sec:sadb_use} documents that mapping.
```

**Advisor — 21 October, 11:28 am:**

> This is a discussion point for the end of the chapter.


### DB-17 · line 13 · tracked edit · 21 October, 11:29 am

**Context near the review anchor (not an exact highlight):**

```latex
 follow describe how the corpus was assembled, audited, pruned, and packaged into its present form.

Methods

The audit and curation effort reported here began from a three-way reconciliation, because the database's three homes did not agree. The author's personal Zotero collection held 538 items, the AARL group librar
```

**Pending advisor edit:**

```text
Added:  Methods
```


### DB-18 · line 14 · comment · 21 October, 11:33 am

**Highlighted source:**

```latex
Methods
```

**Advisor — 21 October, 11:33 am:**

> You don't have to follow this exactly, but you should generally put this in the Chochrane form. https://www.cochrane.org/sites/default/files/uploads/PDFs/guide_to_the_contents_of_a_cochrane_methodology_protocol_and_review.pdf


### DB-19 · line 16 · comment · 21 October, 11:29 am

**Highlighted source:**

```latex
audit
```

**Advisor — 21 October, 11:29 am:**

> ?


### DB-20 · line 16 · comment · 21 October, 11:37 am

**Highlighted source:**

```latex
The author's personal Zotero collection held 538 items, the AARL group library's Biology subtree held 911 items across 27 collections, and the Airtable base held 99 papers. Matching record by record on DOI and normalized title showed that the personal collection was essentially a subset of the group library (all 483 personal DOIs present) and that the Airtable base was a small curated subset rather than a mirror: 95 of its 99 papers resolved to records in the group library's Biology subtree, while 443 personal items had never been entered at all. The gap was imported in verified batches, with every written record re-read and checked after writing, bringing the Papers table to 532 records. The corpus then kept growing while curation caught up: citation mining in the style of Research Rabbit contributed a Crossref-verified digest of recent papers, newly saved references flowed in from the laboratory's Zotero libraries, and a bulk pass imported 372 papers auto-classified from the laboratory's staged collection of PDFs.
```

**Advisor — 21 October, 11:37 am:**

> This is super confusing. We don't need all these details. The airtable base is the database...it can't 'start' with anything if you 'start from the beginning'. Also, starting from the beginning does not mean you tell a chronological story of everything that happened for the database curation. Start with the finished product, then tell us what we need to know about how it could be re-created. 


### DB-21 · line 18 · comment · 21 October, 11:37 am

**Highlighted source:**

```latex
note
```

**Advisor — 21 October, 11:37 am:**

> what are notes?


### DB-22 · line 18 · comment · 21 October, 11:39 am

**Highlighted source:**

```latex
A note is written from the paper's own text, from the abstract at minimum and from the full text when available, never from the curator's recollection of the paper;
```

**Advisor — 21 October, 11:39 am:**

> You are giving details before giving us the big picture. What does that database look like? What is in it? How is it organized?


### DB-23 · line 22 · comment · 21 October, 11:40 am

**Highlighted source:**

```latex
approved
```

**Advisor — 21 October, 11:40 am:**

> approved by who?


### DB-24 · line 22 · comment · 21 October, 11:40 am

**Highlighted source:**

```latex
left as off-topic or orphaned
```

**Advisor — 21 October, 11:40 am:**

> how was this decided?

**Advisor — 21 October, 11:41 am:**

> how were they identified?


### DB-25 · line 24 · comment · 21 October, 11:41 am

**Highlighted source:**

```latex
the laboratory
```

**Advisor — 21 October, 11:41 am:**

> don't you mean you? Take credit in your dissertation.


### DB-26 · line 24 · comment · 21 October, 11:42 am

**Highlighted source:**

```latex
curation tool, not a reader
```

**Advisor — 21 October, 11:42 am:**

> huh?


## 60-sns.tex

### SNS-01 · line 6 · comment · 21 October, 11:48 am

**Highlighted source:**

```latex
Two definitions cover most of the neuroscience vocabulary used below. First, the locomotor controller is a two-level central pattern generator: a rhythm generator sets the timing of the step cycle, and a pattern-formation layer shapes that timing into phase-appropriate activation of antagonist muscle groups \citep{rybak_modelling_2006}. Each leg's rhythm generator is a pair of mutually inhibiting antagonist populations, one extensor and one flexor, so that only one is strongly active at a time. Second, in the non-spiking (graded) formulation used for most of this work, each neuron is electrically a nonlinear RC circuit: synaptic current charges a membrane capacitance through a leak conductance, and the resulting membrane voltage, tens of millivolts, is both the neuron's state and its output \citep{nourse_sns_2023}. Graded chemical synapses couple these units: a depolarized presynaptic cell opens a conductance in the postsynaptic membrane whose sign makes the connection excitatory or inhibitory. The sensory channels that close the loop are the muscle afferents, Ia signaling muscle length and velocity, Ib signaling muscle force, and group II signaling length, together with heel and toe contact sensors.
```

**Advisor — 21 October, 11:48 am:**

> So many details and I don't understand what the goal is. What is the gap in knowledge? What did you do to fill it?

**Ben's reply (author context, not an advisor instruction or independent verification) — 21 October, 2:09 pm:**

> Gaps: Gait2392_simbody with CPG in Mujoco-SNS Toolbox
> - Gait2392_robotbody (BPAs and BPA routing) with CPG
> - Block-programming on a fast physics simulator (Animatlab no longer supported)
> - Robot solidworks assembly with BPAs from Matlab and a neural controller
> - Li walker with Deng neural controller and foot contact feedback.


### SNS-02 · line 8 · comment · 21 October, 12:45 pm

**Highlighted source:**

```latex
the reflex-based and learned alternatives \citep{geyer_musclereflex_2010, song_deep_2021} are treated there.
```

**Advisor — 21 October, 12:45 pm:**

> huh?


### SNS-03 · line 10 · comment · 21 October, 12:47 pm

**Highlighted source:**

```latex
tooling
```

**Advisor — 21 October, 12:47 pm:**

> tooling tends to be a manufacturing term. Software Tool Development is more like what you are looking for here.


### SNS-04 · line 10 · tracked edit · 21 October, 12:48 pm

**Context near the review anchor (not an exact highlight):**

```latex
hapters~\ref{ch:methods}--\ref{ch:torque}.

This chapter concerns tooling: its object is a verified virtual neuromechanical laboratory for the BPA-actuated robot. In this laboratory a controller can be built in any of several simulation engines, checked against its own re-implementations in the other engines and agains
```

**Pending advisor edit:**

```text
Added: virtual
```


### SNS-05 · line 10 · tracked edit · 21 October, 12:48 pm

**Context near the review anchor (not an exact highlight):**

```latex
tooling: its object is a verified virtual neuromechanical laboratory for the BPA-actuated robot. In this laboratory a controller can be built in any of several simulation engines, checked against its own re-implementations in the other engines and against the published behavior of the models it descends from, and then 
```

**Pending advisor edit:**

```text
Changed: that
to this
```


### SNS-06 · line 10 · comment · 21 October, 12:49 pm

**Highlighted source:**

```latex
built in any of several simulation engines
```

**Advisor — 21 October, 12:49 pm:**

> why would we need several?

**Ben's reply (author context, not an advisor instruction or independent verification) — 21 October, 2:04 pm:**

> Well, Mujoco doesn't have as good of a physics engine as HYFYDY, for example, and is limited to rigid tendons. Python also is harder to use because it is code without a GUI. Simscape multibody allows for importing parts from solidworks. Some sort of block-programming language like Simulink is a natural choice for building our style of connectome. Mujoco-SNS Toolbox allows for real-time control and computer in the loop commands.


### SNS-07 · line 10 · comment · 21 October, 12:53 pm

**Highlighted source:**

```latex
answer
```

**Advisor — 21 October, 12:53 pm:**

> This is all super awkward. You need to identify the gap in knowledge and how you fill that gap in knowledge. Much of what you wrote here should be saved for discussion, where you talk about the implications of the work.

**Ben's reply (author context, not an advisor instruction or independent verification) — 21 October, 2:00 pm:**

> I thought that was for research hypothesis not for "here are some tools we made".


### SNS-08 · line 14 · tracked edit · 21 October, 3:34 pm

**Context near the review anchor (not an exact highlight):**

```latex
e closing with the new tools and the program they enable for Chapter~\ref{ch:futurework}.

\section{The AnimatLab Bipedal Walker}\label{sec:prior_walker}

As preliminary work of this dissertation, built before the simulation toolchain of this chapter existed, the author and Morrow constructed a neuromechanical model of
```

**Pending advisor edit:**

```text
Deleted: Preliminary Results:
```


### SNS-09 · line 20 · comment · 21 October, 12:59 pm

**Highlighted source:**

```latex
verification-first approach
```

**Advisor — 21 October, 12:59 pm:**

> I don't understand the connection.


### SNS-10 · line 24 · comment · 21 October, 1:02 pm

**Highlighted source:**

```latex
build order
```

**Advisor — 21 October, 1:02 pm:**

> huh?


### SNS-11 · line 26 · comment · 21 October, 1:00 pm

**Highlighted source:**

```latex
\begin{figure}[htbp]
```

**Advisor — 21 October, 1:00 pm:**

> Overlapping text in figure


### SNS-12 · line 35 · tracked edit · 21 October, 1:04 pm

**Context near the review anchor (not an exact highlight):**

```latex
peline}
\end{figure}

\subsection{Existing Tools and Data Resources}\label{sec:sns_existing_tools}

OpenSim supplies the musculoskeletal modeling side, the robot's source model, and the moment-arm and torque standards the conversions are audited against \citep{delp_opensim_2007, seth_opensim_2018}. MuJoCo supplies the 
```

**Pending advisor edit:**

```text
Deleted: The laboratory stands on tools this program did not build.
```


### SNS-13 · line 37 · comment · 21 October, 1:06 pm

**Highlighted source:**

```latex
The SADb of Chapter~\ref{ch:database} is the data resource among the existing tools \citep{bolen_sensory_2023}. Its segments are the five-table Airtable base (Papers, Feedback, Models, Review Papers, and Recorders), the laboratory's Zotero libraries that feed it, citation mining in the Research Rabbit style that expanded it, and the laboratory's own packaged explorer application that distributes it. The SADb framed this chapter's benchmarks and validation measurements: its pathway records fix what each afferent class encodes and what each connection rule does, so the sensory channels, the reflex topologies, and the acceptance phenomena against which the controllers are validated (the deafferentation and stimulation results cited in Section~\ref{sec:walker_variant_results}) are looked up rather than re-derived from the primary literature.
```

**Advisor — 21 October, 1:06 pm:**

> This whole paragraph seems unnecessary. 


### SNS-14 · line 39 · comment · 21 October, 1:08 pm

**Highlighted source:**

```latex
\subsection{New Tools}\label{sec:sns_new_tools}
```

**Advisor — 21 October, 1:08 pm:**

> Wait, that is all you are saying about the other tools? It is hard for the reader to understand the limitations and gaps without a deeper description of them.


### SNS-15 · line 41 · comment · 21 October, 1:11 pm

**Highlighted source:**

```latex
Simulink blocks
```

**Advisor — 21 October, 1:11 pm:**

> Did you compare and discuss against the libraries myself and Fletcher made here?

**Ben's reply (author context, not an advisor instruction or independent verification) — 21 October, 2:00 pm:**

> No. You've never encouraged us to use Simulink.


## What the advisor's new notes require

The introduction has newly annotated undefined references to numbered gaps, unexplained spiking-neuron mirrors, an overly specific roadmap, and a weak connection to the Biomimetic Humanoid Project. Name the gaps directly and present them in chapter order. Explain Ben's work as a continuation of the humanoid project, rather than inserting adjacent laboratory projects without a clear connection.

The final background batch also identifies CPG-specific material misplaced in the broader computational-neuroscience section (BG-31 and BG-43–45). Move or condense it rather than repeating it in both sections. The new background notes target the motor drawing, the graphical abstract's usefulness, Gait2392's name, the moment-arm discussion, joint-center modeling, run-on sentences, and the chronology of the master's model. State what is compared and how it supports the present work. Replace abstractions such as “first-class,” “definitions are stressed,” and “written against” with concrete descriptions. Do not simply add explanatory paragraphs after every “huh?”; replace or remove the unclear statement and its unnecessary detail.

SADb should start with the finished resource: what it contains, its organization, and who can use it. Then give only the methods needed to reproduce the curation, adapting the advisor's Cochrane suggestion without pretending it was a registered systematic review. Verify the two conference presentations, AI-assisted workflow, decision authority, and actual inclusion criteria from evidence. Ben's role should be explicit. Drop the import chronology and operational record counts unless needed to interpret the final corpus.

For SNS, the old suggestion to compress AnimatLab to one generic paragraph is superseded by Ben's requested causal narrative. The account must identify what Ben combined, what the air-stepping result shows, what failed under hand tuning, and why the new tools were needed. A compact capability/evidence table can replace repeated inventories. Verification should identify the expected behavior, the test or comparison, the observed outcome, and its limit. Do not expand this into a catalogue of all experiments.

## GitHub torque update checked in this pass

Fetch succeeded. The active branch is KneeTestSetup_BenBo_stw. During the pass HEAD and origin both reached **65b05aae161d6a4e87f561a38e57dbae90ec6f95** (22 October, 12:00 Pacific), committed by another session/user; this agent did not commit or pull over local changes. Re-fetch at the start of ZCode work because other machines remain active.

- **fd1c13f4**, 22 October 02:22: updates Knee_Torque_Test/ProofFinal/Main.tex:540–542 with the six-point extensor cam-route sweep, explicitly pending reaction-point geometry and adjoint torque reduction. Read as manuscript status, not independently verified final torque values. Chapter 4 should not continue to say no full angle sweep exists if the test record confirms this sweep.
- The **same manuscript still says “full sweep ... pending” at line 556 and “first ... points ... near full extension” at line 560**. Resolve this internal inconsistency against the actual test sheets before adapting it. Do not copy both claims.
- **7862e4e2**, 22 October 04:31: adds ExtTest_3_01..03 records, updates Results_table_20mm.xlsx, Knee_Extensor_20mm.m, and Ext20mm_compare_20261022.png, and adds the new BiPam_pulley2_20261022.m plus sweep33/escalation records. These are newer than the six-point paragraph; inspect their provenance before selecting the current validation status.
- The 21 October commits also update torque and margin surfaces and the transmission doctrine. PULLEY_SWEEP_DIRECTIVE_20261021.md is the current scope: nBPA and tackle parts each span 1–3, re-solve actuator sizes, account for moving-via load and Bowden compliance, and audit OpenSim path artifacts. The existence of CSVs/mats is not proof that every configuration is feasible or implemented in the simulation model.

The dissertation's 30-torque.tex also differs between live Overleaf and local ProofFinal. Neither was overwritten in this pass. Integrate in a small evidence-backed patch against fresh Overleaf text, not by replacing the chapter with the torque manuscript.

## Figure synchronization in this pass

Uploaded the 11 revised PDF assets to their existing figs/10-introduction, figs/15-background, figs/60-sns, and figs/70-overall folders. Kept editable PPTX, EPS, alt text, and generators local. Chapter 3 assets were not uploaded or edited.

Applied seven exact-fragment changes to figure references, sizing, captions, and the sentence introducing the functional neuron view. No whole chapter was uploaded. Subsequent layout and technical checks led to nine further exact-fragment changes to the two combined figures: shorter captions and 0.94-textwidth panels retain figure text above 10 pt, and “fix axial displacement” precisely identifies the transverse constraint. The separate patches are recorded in the JSON evidence. Advisor tracked edits and comment status were not accepted, rejected, or resolved.

No report, prompt, or administrative file was uploaded to Overleaf; the fresh pre-upload ZIP contains no chatgpt_report file. The source snapshot revealed six divergent chapters (introduction, background, torque, SADb, SNS, and Discussion). Local drafts and live advisor edits must be reconciled, not treated as interchangeable. A fresh source ZIP is evidence at its download time, not a permanent canonical version.

## ZCode execution contract

Use Notes/new_session_prompt_20261022.md. It replaces the stale broad task list with a focused reduction/revision pass, explicit prose and caption limits, the new SNS structure, small patch batches, and evidence requirements. The old prompt and report are preserved locally in Notes/advisor_review_20261022. No chapter reduction has been claimed completed in this pass: ZCode is assigned that heavy revision work.

## Final synchronization checks

- **Live Overleaf:** 295-page current preview; Errors 0, Warnings 4. Both new combined-figure size warnings were eliminated. The four remaining warning cards concern pre-existing floats and a database table outside this figure-sync scope. The advisor's live source was not replaced with the different local draft.
- **Local verification:** ProofFinal/build_advisor_figures_20261022_final/main.pdf, 298 pages; four MiKTeX passes exit 0, zero undefined references/citations, and zero missing files. This local draft differs from the live review copy, so its page count is not a shortening result and must not be presented as the live Overleaf page count.
- Final live SNS and Discussion sources match the fresh pre-upload source plus only the intended figure patches (UTF-8 BOM ignored for comparison). The earlier copied introduction matches its intended single caption patch. No general prose rewrite was applied on Overleaf.
- The final rendered pages for Figures 6.6 and 7.2 pass independent visual and technical review. Minimum editable figure text is 10.152 pt. Their lowest text lies at 672.51 and 671.01 pt on a 792 pt page, respectively. No clipping remains.
- Chapter 3's 19 protected files remain byte-identical. Reports, prompts, native editing sources, and evidence stayed local. No advisor comments were resolved, no tracked edits were dispositioned, and no commits or pushes were made by this agent.
- GitHub was fetched again at completion; HEAD and origin still match 65b05aae. The latest torque-manuscript source commit remains fd1c13f4. No remote change was merged over the dirty working tree.
- Overleaf source/PDF export attempts later in the pass were blocked by Chrome. Final verification therefore uses the visible live PDF, compile UI, and exact source-copy comparisons, plus the independent local build. A completed earlier PDF download was inspectable, but it is not claimed as the final live PDF.

### Per-directive disposition

| Ben's request | Completed action / evidence |
|---|---|
| Review advisor contents, not the whole dissertation | Review Overview capture, exact anchors, pending edits, and relevant reply threads; no new full quality review. |
| Update local CHATGPT_REPORT | Current 73-thread/51-edit inventory above; previous report preserved locally. |
| Insert changed figures | Eleven PDF assets uploaded to the existing chapter figure folders; exact-fragment integration patches and Overleaf proof image saved in Notes/advisor_review_20261022/. |
| Prepare ZCode for reduction and new SNS structure | Replaced Notes/new_session_prompt_20261022.md; explicit AnimatLab → gap → tools/results → verification sequence, limited future-work handoff, and no mandatory Methods/Results/Discussion headings. |
| Avoid more bloated prose and captions | Prompt requires measured cuts, brief captions, comment-specific repairs, and independent review; both combined-figure captions shortened in this pass. |
| Check latest torque work | Fetch and commit/source inspection recorded above, including contradictory sweep status and newer test/sweep records. No unverified numerical transplant. |
| Keep active advisor review intact | No chapter-file uploads, wholesale replacements, comment resolution, or tracked-edit acceptance/rejection; future prose changes remain targeted local proposals. |
| Never save reports on Overleaf | No report/prompt uploaded; only figure PDFs and targeted figure-related source fragments changed there. |
| Protect isometric-force/MATLAB figures | All 19 protected Chapter 3 files retain their hashes. |

Start the new ZCode session using the entire prompt file. It is a handoff for work still to do, not a claim that SADb/SNS reduction or torque revision has already been completed.


## Final Revision C sweep — 2026-10-23

The expanded live Overleaf Review Overview was checked against this ledger. The inventory remains **73 advisor comment threads and 51 pending advisor tracked edits**: Introduction 34 items, Background 49, SADb 26, SNS 15. Ben’s replies and the four Ben-only appendix/main items remain separate. No new open advisor items were found. No Overleaf source, comment resolution, or tracked-edit disposition was changed.

Every advisor ID now has a local disposition and source-line evidence in `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C/COMMENT_DISPOSITIONS.md`. The original expanded entries are retained in `review/final_sweep_20261023/recorded_comments.json`. These are editorial dispositions, not claims of advisor approval.

This pass makes the SNS knowledge gap match all five combinations in Ben’s reply, keeps the interconnected bilateral RG construction in the causal opening, adds an early SADb evidence-gap statement, and puts the bubble chart in the SADb main chapter. The current app payload was independently read: it now includes all 1,393 records, so the old paper-only app limitation is superseded. The separate monthly scraper remains under development.

The completed physical-test plots from the latest commits are inserted; no spreadsheet sanity-check torque values are used. Four full fit tables move to Appendix C, with a compact mean-error summary in Chapter 4. The cat caption is shortened, and Figures 2.8 and 2.9 share one page. The active bipulley run is not reported as completed. Original sources and the protected force chapter/figures remain unchanged by this pass.


## Revision C uploaded; Chapter 4 tables restored — 2026-10-23

Uploaded Revision C into the existing Bolen-Dissertation project at Ben's explicit authorization. Replaced individual chapter files and main/support files, preserving folders and review threads. Updated the roadmap, Background figures, torque figures including final physical-test exports, SADb screenshots, SNS figures, Discussion schematics, and Appendix B circuit figures. Existing unchanged force-chapter and reference figures remain. No reports or review artifacts uploaded.

Chapter 4 now contains separate complete pre-optimization, optimized pinned-flexor, and optimized pinned-extensor GoF tables (4.1–4.3), plus a two-row flexor/extensor correction table (4.4), with four parameter columns chi0–chi3 as Ben confirmed. The earlier cross-configuration fit remains in Appendix C. No GoF results were recalculated or changed. The extensor-identified wrapping coefficient is explicitly distinguished from the pinned-flexor fit.

Final Overleaf compilation: Errors 0, Warnings 0; inherited typesetting warnings remain. Overleaf PDF has 259 physical pages; local MiKTeX verification PDF has 257. Both show Tables 4.2 and 4.3 on printed page 96. The Overleaf table page was visually inspected. No source text or comments were added to the protected isometric-force chapter.

Advisor review: 40 straightforward comments resolved; 34 substantive threads remain open, each with a new response describing the revision. This does not claim advisor acceptance. Accepted all 13 remaining grouped tracked-edit markers after replacement; Overview shows no pending tracked edits. The upload regrouped the original 51 advisor tracked edits and Ben's 3 tracked deletions, so grouped-marker counts differ from the original inventory. review_actions.json retains thread text, submitted replies, resolutions, and accepted groups.

The source-ZIP download event could not be verified. A labeled pre-upload history checkpoint was created instead. User authorization covered uploads and routine advisor replies/resolutions. No commit or push was performed.

Evidence: `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C/review/overleaf_upload_20261023/`.


## Ben's figure/prose correction sweep — 2026-10-23

Working copy: `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C`. Originals in ProofFinal remain unchanged. Ben explicitly authorizes Chapter 3 layout changes in this pass. ZCode-owned items are requests to verify, not facts to invent. Progress is updated in this checklist and `review/user_sweep_20261023/checklist.json`.

- [x] U01 [DONE] Remove the quoted stiffness-pair nonuniqueness/external-overshoot-selection paragraph; do not substitute a defensive Pareto justification. Evidence/status: 30-torque.tex: quoted sentence removed; matching Appendix C selection narrative removed. Before snapshot preserves it.
- [x] U02 [DONE] Replace Fig. 2.3 CAD BPA with supplied inflated/deflated photo; preserve originals, give both camera photos descriptive filenames. Evidence/status: Copied original JPEGs without resampling to figs/15-background/festo_bpas_inflated_deflated.jpg and figs/30-torque/cam_extensor_hardware.jpg; original phone files untouched. Photo wired into Fig.2.3; camera EXIF rotation handled in TeX. Visual page check pending.
- [x] U03 [DONE] Fig. 2.4: uppercase panel letters. Evidence/status: steeleknee.pdf: uppercase A/B native PDF overlays; caption uppercase too.
- [x] U04 [DONE] Fig. 2.5: panel letters at upper left. Evidence/status: 15-background.tex: three top-aligned minipages, letters at upper left.
- [x] U05 [DONE] Reorganize Figs. 2.5–2.6: skeleton, robot PAM routes, human Gait2392 paths together; remove repeated models and awkward letters. Evidence/status: Fig.2.5 now skeleton / robot PAM routes / human Gait2392 paths; split original masters_leg_model.pdf into robot_pam_routes.pdf and human_gait2392_paths.pdf. Fig.2.6 replaced duplicate renderings with moment_arm_methods.pdf.
- [x] U06 [DONE] Resize fig:Fmax_fit; fit Figs. 3.3 and 3.4 within page text area without altering MATLAB plot data. Evidence/status: 20-muscles.tex: 3 figures use page floats, centered, bounded panel width/height. Source EPS/MATLAB plots unchanged. Full page verification pending. Full 5-pass build passed; figures3.2–3.4 rendered and visually checked within page bounds.
- [x] U07 [DONE] Center Table 3.3 and increase text size; preserve every value. Evidence/status: Table model_compare unpivoted into 12 readable rows (model per row), centered, small rather than tiny, no negative adjustwidth. Numeric conservation/page check pending. All 36 numeric results checked against original; centered readable table visually verified.
- [x] U08 [DONE] Add diagram connecting Eq. 4.2 dot-product moment-arm construction to cross-product torque method; reuse appropriate Fig. 2.6 geometry where useful. Evidence/status: figs/15-background/moment_arm_methods.pdf: cross-product torque and dot-product moment-arm projection; Fig.2.6 label preserved and cross-referenced near Eq.4.2. Fixed prior self-orthogonal definitions of p_proj versus p_f in nearby prose; equation unchanged. Diagram visually checked.
- [x] U09 [DONE] Fig. 4.4 panel letters bold. Evidence/status: bktFrame.pdf: native text A–D replaced with bold labels; CAD/image content preserved.
- [x] U10 [DONE] Locate and show bending schematic with Eq. 4.24. Evidence/status: No existing BPA-wrapping schematic was located in searched dissertation/paper assets (bracket-bending assets are a different effect). Added native vector bpa_wrapping_loss.pdf as panel B of existing Fig.4.5, referenced at Eq.4.24, preserving later figure numbering. Diagram states fitted loss is not arc-minus-chord length.
- [x] U11 [DONE] Shorten captions of Tables 4.1–4.3; move essential scoring/fold context to prose, preserve values. Evidence/status: Tables4.1–4.3 captions shortened; point-set restrictions and counts moved before Table4.1; fold membership already in rows; excluded flexor test retained in prose. Numbers unchanged.
- [x] U12 [DONE] Table 4.4 chi0 in meters, with correct unit conversion. Evidence/status: Table4.4 header meters; +3.94mm -> +0.00394m, -6.38mm -> -0.00638m. Other parameters unchanged.
- [x] U13 [DONE] Explain anatomical origin/insertion search bounds, biomimetic routing/RoM tradeoff, and subsequent force/part-failure tradeoff. Evidence/status: Added paragraph before route redesign: anatomical bounds, relaxation for RoM, subsequent lower-force hardware tradeoff; cam paragraph simplified. No new bounds or optimizer outcomes invented.
- [x] U14 [DONE 2026-10-23 ZCode] Fig. 4.11 panel C showing torque magnitude across RoM: ZCode-owned; do not duplicate concurrent plot work. Evidence/status: DONE. Panel (C) (Flx20mm_meetsHuman — measured torque vs corrected 620 kPa prediction vs human BFSH maximum, 15 holds −119..+5°) added to fig:FlxGrp in 30-torque.tex; margins +0.42 N·m (@−105°) to +17.0 N·m (@−14°) stated in the panel-B prose paragraph; caption numbers match the torque paper's fig:Flx20mm record and were re-verified by re-running gen_Flx20mm_3panel_20261023.m (15/15 meet-or-exceed). Asset copied byte-identically into figs/30-torque/. A/B panels and caption preserved; label survives once; build fits (float-too-large resolved at 0.72/0.30 widths). Report: Notes/advisor_revision_20261022/Revision_C/review/u14u21_20261023/REPORT.md.
- [x] U15 [DONE] Add supplied cam-knee hardware photograph to Fig. 4.12. Evidence/status: cam_extensor_hardware.jpg copied byte-identically and included below cam route panels in Fig.4.12. Needs full page visual check.
- [x] U16 [DONE] SADb gap: structured database rather than structured account. Evidence/status: 50-database.tex: structured account -> structured database.
- [x] U17 [DONE] SADb source libraries: Zotero libraries of AARL and Ben, plus papers already in Airtable. Evidence/status: 50-database.tex: Zotero libraries of the AARL and myself, and papers already curated in Airtable.
- [x] U18 [DONE] Read SfN24_Bolen.pdf; replace vague AI claims with evidence-based account of limited help and categorization difficulties. Evidence/status: Read SfN24_Bolen.pdf, extracted text and visually inspected poster. Paragraph names ResearchRabbit, Litmaps, txyz, SciSpace, MirrorThink; categorization: none good, txyz adequate, others insufficient. Limited benefits and substantial correction stated. Existing bolen_using_2024 citation retained.
- [x] U19 [DONE] SNS 6.1: initial Connor/Ben controller had ONE RG controlling both sides; air stepping at NIH2021; Ben added afferent and heel/toe feedback; distinguish later bilateral RG expansion. Evidence/status: 60-sns.tex: initial Connor/Ben model ONE RG both sides; air stepping at NIH2021; distinguishes Ben feedback/contact additions and later bilateral interconnected RG model.
- [x] U20 [DONE] State AnimatLab unsupported, based on old physics engine, and slow in this workflow. Evidence/status: 60-sns.tex: AnimatLab unsupported, older physics engine, slow in this workflow.
- [x] U21 [DONE 2026-10-23 ZCode] ZCode verification: Table 6.1 split verification/limits; identify actual authored models versus existing MuJoCo/SNS tools, external MuJoCo–Simulink bridge attribution and direction, controller ports/demos/curriculum scope, resolved welded joint, actual muscle integration status. Do not guess answers. Evidence/status: DONE (verification pass; three text corrections). (1) Bridge attributed (official MathWorks mujoco-simulink-blockset, MuJoCo 2.3.7, AARL keyframe patch + a/b/c tests — BRIDGE_REPORT.md) and direction stated (MuJoCo model embeds inside Simulink; ctrl in, sensors out) in 60-sns.tex §6.4.3. (2) "92 actuators 5 s closed-loop run" corrected to open-loop without afferent return (E2 e2_cpg_open_loop.m); closed-loop credit moved to the E1 reflex experiment. (3) Table 6.1 row-4 limit now states the imported lower-body mechanism's ankle joints remain rigid (digest_2026-09-30: AH_001 ankle DOF open; knee DOF present in the pivoted BA_003 import). Verified-no-edit-needed: authored-vs-existing attributions, W2L rhythm-core numbers (1.000 s, 10/10, r≥0.99996 — CROSS_PLATFORM_REPORT §1.3), Appendix D claims incl. anchor −161.57, ports/demos/curriculum scope, muscle-integration account. Welded-joint claim stays dead (zero hits). Supervisor four-gate GO. Report: Notes/advisor_revision_20261022/Revision_C/review/u14u21_20261023/REPORT.md.
- [x] U22 [DONE] Remove unsupported Hunt/Young neuromechanical Simulink-library comparison; flag prior advisor reply for correction. Evidence/status: Unsupported comparison removed from manuscript. Earlier Overleaf reply edited in place: it no longer claims a Hunt/Young library comparison and leaves comparison with the actual libraries pending verification.
- [x] U23 [DONE] Define W2L at first use, grounded in actual model naming. Evidence/status: Defined Walker Two-Layer (W2L) at first mention; grounded in spinal/w2l_cpg/README.md / Walker_2_Layer_CPG_BilateralRG model naming.
- [x] U24 [DONE] Gap includes tested neuromechanical models with Gait2392-scale actuator counts and degrees of freedom; scope novelty honestly. Evidence/status: 60-sns.tex gap now explicitly includes many degrees of freedom and 92 musculotendon actuators in Gait2392. Does not claim no such models exist anywhere.
- [x] U25 [DONE] Fresh full build, render changed pages, conservation/line-evidence audit, clean upload bundle, sync authorized finished changes to Overleaf. No report uploads. Evidence/status: Fresh 5-pass local build passed. All 36 Table3.3 results and all Table4.1–4.3 data bodies conserved; Chapter3 source plots byte-identical. Changed pages rendered and visually inspected. Uploaded six changed chapters and eight changed/new figure assets individually to Overleaf. Live compilation: Errors0, Warnings0; existing typesetting warnings remain. Live photo and Fig3.3 page inspected; plot and caption fit. Prior misleading reply corrected in place. Reports not uploaded. ZCode-owned U14/U21 remain outstanding.


### Handoff checkpoint: user correction sweep

All edits are in Revision_C, NOT ProofFinal. Nothing from this new sweep has been uploaded to Overleaf yet. Do not run older one-shot revision scripts or packaging scripts: they would overwrite current work or assert the old table/protection state. Current sources are UTF-8; read/write explicitly with UTF-8. Bibliography remains untouched. Source MATLAB/experimental plots remain untouched; Chapter3 layout and table reformatting were explicitly authorized by Ben in this request.

Before snapshot: `Revision_C/review/user_sweep_20261023/before/` (all chapters and figs). Checklist: same folder `checklist.json`; evidence includes poster text/render, asset renders, and original copies. New scripts in root tmp: `revision_c_user_edit.py` (ONE-SHOT, already applied), `revision_c_user_assets.py` (vector PDF authoring, do not rerun on modified assets without restoring snapshot), `revision_c_user_layout2.py` (already applied), `revision_c_user_progress.py` (this checkpoint). Generating vector diagram code is in the assets script. Existing full build helper: `tmp/revision_c_build.py`.

ZCode priorities: Figure4.11 panel C, Table6.1 factual verification and attribution/direction questions, final layout QA, and correcting the prior Overleaf Hunt/Young reply. Do not invent a benchmark or claim Hunt/Young shared a neuromechanical Simulink library. Ben's original technical questions remain in U21 above. Fig4.5 now has wrapping schematic panel B to avoid renumbering the downstream torque figures.


### User correction sweep: final uploaded status

The pre-upload checkpoint above is superseded by this status. All six revised chapter sources and eight figure assets are now on Overleaf. Compile panel shows Errors 0 and Warnings 0; inherited typesetting warnings remain. The live BPA photograph and resized Figure 3.3 were visually checked. The earlier Hunt/Young advisor reply was corrected in place. No new comments were resolved in this sweep; substantive threads remain open.

Completed: U01–U13, U15–U20, U22–U25. ZCode-owned: U14 (Figure4.11 panel C) and factual verification under U21 (Table6.1 attribution, MuJoCo–Simulink bridge direction, actual imported mechanism/actuator integration, and scope of demos/curriculum). U21 layout and the clearly false welded-knee statement were corrected; do not infer that the remaining technical questions have been answered. No performance benchmark was invented.

Latest local review PDF: `Revision_C/Revision_C_user_sweep_20261023.pdf`. Clean upload folder: `Revision_C/OVERLEAF_UPLOAD`. Latest ZIP: `Revision_C/Revision_C_user_sweep_Overleaf.zip`. Evidence, before snapshots, complete diff, checklist, numeric conservation, and screenshots: `Revision_C/review/user_sweep_20261023/`. Original phone photos and MATLAB plot files are preserved. Photo copies have descriptive names. No commit/push or ProofFinal edits.


## U14 + U21 ZCode pass — 2026-10-23

Both ZCode-owned items are done locally against the fresh Overleaf ZIP (downloaded 08:25, extracted to tmp/overleaf_live_20261023). U14: panel (C) added to Fig. 4.11 (measured torque magnitude vs human BFSH across the 15-hold RoM, paper numbers of record re-verified by re-running the generator). U21: bridge attribution/direction stated, the 92-actuator run corrected to open-loop, and the Table 6.1 imports row now carries the ankle-rigid limit; all other table claims verified against campaign records with no correction needed. Supervisor four-gate GO. Local verification build 256 pp, 0 errors, 0 undefined refs; renders of printed pages 101/122/129 checked. NOTHING uploaded; no comments resolved; no commits or pushes. Full record: Notes/advisor_revision_20261022/Revision_C/review/u14u21_20261023/REPORT.md.
