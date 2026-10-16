# Dissertation comments and tracked changes

Read from [Bolen-Dissertation on Overleaf](https://www.overleaf.com/project/644dc82bcd6a5481e3c3a8bb) on October 15, 2026.

## Scope and reference conventions

- 46 open comment threads, containing 50 comment messages and replies.
- 15 pending tracked-change entries across ten files.
- 32 resolved comment threads, preserved separately below with their replies and resolution records.
- Line numbers refer to the current live Overleaf LaTeX source, not PDF pages or the local dissertation copy. Quoted source preserves LaTeX commands. Reviewer wording is transcribed verbatim.
- 24 open threads have no surviving source highlight and sit at the current end of their file. Their current anchor line is reported for navigation, but their original highlighted text, section, and line cannot be recovered from the review interface. No original location is inferred from the message.
- Overleaf groups some older deletions into a single entry at the end of a file. The full displayed deletion text is included, but the original locations of individual fragments are unavailable. Text in these grouped entries may be concatenated.
- Resolved threads do not expose source line anchors in the archive. Saved quoted text is included where present; blank quotes and unavailable historical locations are explicitly identified.
- Only CHATGPT_REPORT.md was replaced. No Overleaf source, comments, tracked-change decisions, or other local files were modified.

## Open comments

### 10-introduction.tex

#### C01 — ajh26 — 14 October, 10:59 am

**Section:** Introduction → Motivation

**Source line:** 8.

**Highlighted text:**

```text
Artificial
```

**Comment — ajh26 — 14 October, 10:59 am:**

```text
This paragraph does a great job identifying one of the gaps of knowledge you worked to fill, but there are a few others you can worked on too.
```

#### C02 — ajh26 — 14 October, 11:44 am

**Section:** Introduction → Problem Statement

**Source line:** 22.

**Highlighted text:**

```text
Second
```

**Comment — ajh26 — 14 October, 11:44 am:**

```text
Here are the other ones. These can probably be moved up to the motivation section with more explanation, like the other 1, and then eliminate this section.
```

### 15-background.tex

#### C03 — ajh26 — 14 October, 1:36 pm

**Section:** Background → Human Lower-Limb Biomechanics and the OpenSim Benchmark

**Source line:** 103.

**Highlighted text:**

```text
Figure~\ref{fig:aarl_opensim_comparison} places the laboratory's bipedal skeleton beside the OpenSim muscle-routing reference developed for it.
```

**Comment — ajh26 — 14 October, 1:36 pm:**

```text
feels out of place
```

#### C04 — ajh26 — 14 October, 1:50 pm

**Section:** Background → Neural Control of Locomotion: CPGs and Sensory Feedback

**Source line:** 142.

**Highlighted text:**

```text
begin
```

**Comment — ajh26 — 14 October, 1:50 pm:**

```text
this figure isn't really readable.
```

#### C05 — ajh26 — 14 October, 1:47 pm

**Section:** Background → Synthetic Nervous Systems and the Sensory Afferent Database

**Source line:** 154.

**Highlighted text:**

```text
others
```

**Comment — ajh26 — 14 October, 1:47 pm:**

```text
like us?
```

#### C06 — ajh26 — 14 October, 1:48 pm

**Section:** Background → Synthetic Nervous Systems and the Sensory Afferent Database

**Source line:** 154.

**Highlighted text:**

```text
(AARL), the rat hindlimb modeling program that supplies this dissertation's controller architecture
```

**Comment — ajh26 — 14 October, 1:48 pm:**

```text
?
```

**Reply — ajh26 — 14 October, 1:49 pm:**

```text
this is weirdly formatted
```

#### C07 — ajh26 — 14 October, 2:03 pm

**Section:** Background → Synthetic Nervous Systems and the Sensory Afferent Database

**Source line:** 154.

**Highlighted text:**

```text
survey its past, present, and future, and physical testbeds in this lineage now extend from insect-scale models to a robotic cat hindlimb 
```

**Comment — ajh26 — 14 October, 2:03 pm:**

```text
awkward
```

#### C08 — ajh26 — 14 October, 2:06 pm

**Section:** Background → Simulation Tooling: From OpenSim to MuJoCo

**Source line:** 162.

**Highlighted text:**

```text
Neuromechanical
```

**Comment — ajh26 — 14 October, 2:06 pm:**

```text
Does this need to be its own subsection?
```

**Reply — ajh26 — 14 October, 2:06 pm:**

```text
The transition here is awkward
```

#### C09 — ajh26 — 14 October, 2:52 pm

**Section:** Background → Testing Neuromechanical Hypotheses on Robots: The Ethical Case

**Source line:** 185.

**Highlighted text:**

```text
possible
```

**Comment — ajh26 — 14 October, 2:52 pm:**

```text
Something that should probably be worked in though is that such models can only be possible because of animal experiments. This make the ethical framework much more murky in my opinion.
```

#### C10 — ajh26 — 10 September, 11:41 am

**Section:** Original section unavailable; detached/end-of-file anchor.

**Source line:** 187 (current end-of-file anchor only; original line unavailable).

**Highlighted text:**

Unavailable: no surviving highlighted source text is displayed.

**Comment — ajh26 — 10 September, 11:41 am:**

```text
My opinion is to feel not guilty about borrowing images from other papers for a dissertation background section (Just make sure to cite it!). That way it won't take a ton of time and is helpful in getting your readers up to speed without sending them directly to the papers.
```

#### C11 — ajh26 — 10 September, 11:46 am

**Section:** Original section unavailable; detached/end-of-file anchor.

**Source line:** 187 (current end-of-file anchor only; original line unavailable).

**Highlighted text:**

Unavailable: no surviving highlighted source text is displayed.

**Comment — ajh26 — 10 September, 11:46 am:**

```text
At the end of each section/paragraph, you can add a sentence directly connecting it to your PhD work. Example: The work in this dissertation will help develop biomimetic robots that more accurately captures animal joint torque profiles.
```

### 20-methods.tex

#### C12 — ajh26 — 14 October, 3:02 pm

**Section:** Materials and Methods → Balance Platform and Inverted-Pendulum Testbeds

**Source line:** 387.

**Highlighted text:**

```text
\section{Balance Platform and Inverted-Pendulum Testbeds}\label{sec:balance_testbeds}
```

**Comment — ajh26 — 14 October, 3:02 pm:**

```text
I agree with Leo now that this dissertation has gotten much bigger. This is a lot of methods to keep in your head without getting some results along the way. A chapter on just the muscles (Your paper); A chapter on the torque; a chapter for the ankle; a chapter for the SNS; A chapter for the database. Include methods, results, and a  discussion in each of the chapters. Then you can discuss the work as a whole in the final chapter.
```

#### C13 — ajh26 — 14 October, 3:14 pm

**Section:** Materials and Methods → Balance Platform and Inverted-Pendulum Testbeds

**Source line:** 389.

**Highlighted text:**

```text
The testbeds deliberately reduce the robot to a single actuated joint and a single sensory question: which measurement of body orientation, fed through which feedback law, holds the leg upright on a moving support? This reduction allowed vestibular and proprioceptive feedback strategies to be exercised in hardware and provides a baseline for the full-body comparisons defined in Chapter~\ref{ch:futurework}.
```

**Comment — ajh26 — 14 October, 3:14 pm:**

```text
I would avoid talking too much about the balance questions in your work here, and stick to the design aspects.
```

#### C14 — ajh26 — 14 October, 3:20 pm

**Section:** Materials and Methods → Balance Platform and Inverted-Pendulum Testbeds

**Source line:** 397.

**Highlighted text:**

```text
The final balance platform stands the inverted pendulum, its ankle driven by an antagonist BPA pair, on a motorized platform that oscillates the support surface beneath the foot. An IMU on the tibia serves as a pseudo-vestibular sensor, feeding the same proportional law. This is the configuration in which balance was demonstrated both on flat ground and under the oscillating support (Section~\ref{sec:results_balance}).
```

**Comment — ajh26 — 14 October, 3:20 pm:**

```text
I question whether the inverted pendulum work should be in the dissertation. Not much 'scientific' happened here. Mostly just some interesting design work. It also isn't clear how each design follows from each other, what the overall purpose is, and there is a significant lack of details in this section to really constitute a methods section.
```

#### C15 — ajh26 — 14 October, 4:09 pm

**Section:** Materials and Methods → Prior Neuromechanical Simulation: The AnimatLab Bipedal Walker

**Source line:** 406.

**Highlighted text:**

```text
\section{Prior Neuromechanical Simulation: The AnimatLab Bipedal Walker}\label{sec:prior_walker}
```

**Comment — ajh26 — 14 October, 4:09 pm:**

```text
This neuromechanical simulation stuff lacks 'direction'. You need a clear storyline of what is being built and why up front. Prior and preliminary work goes in the background. Everything that isn't background is to be treated as part of your dissertation work. 
```

#### C16 — ajh26 — 14 October, 4:10 pm

**Section:** Materials and Methods → Preliminary Neuromechanical Simulation Workflows → MuJoCo Conversion and Synthetic BPA Knee

**Source line:** 463.

**Highlighted text:**

```text
scripted model repair
```

**Comment — ajh26 — 14 October, 4:10 pm:**

```text
I don't know what this means
```

#### C17 — ajh26 — 14 October, 4:11 pm

**Section:** Materials and Methods → Preliminary Neuromechanical Simulation Workflows → MuJoCo Conversion and Synthetic BPA Knee

**Source line:** 463.

**Highlighted text:**

```text
documented software environment
```

**Comment — ajh26 — 14 October, 4:11 pm:**

```text
I don't know what this means
```

#### C18 — ajh26 — 14 October, 4:12 pm

**Section:** Materials and Methods → Preliminary Neuromechanical Simulation Workflows → MuJoCo Conversion and Synthetic BPA Knee

**Source line:** 463.

**Highlighted text:**

```text
ported
```

**Comment — ajh26 — 14 October, 4:12 pm:**

```text
what do you mean ported? Do you just mean you implemented the BPA model in Python?
```

#### C19 — ajh26 — 14 October, 4:12 pm

**Section:** Materials and Methods → Preliminary Neuromechanical Simulation Workflows → MuJoCo Conversion and Synthetic BPA Knee

**Source line:** 463.

**Highlighted text:**

```text
coupled to MuJoCo tendon kinematics
```

**Comment — ajh26 — 14 October, 4:12 pm:**

```text
I don't know what this means
```

#### C20 — ajh26 — 14 October, 4:16 pm

**Section:** Materials and Methods → Preliminary Neuromechanical Simulation Workflows → MuJoCo Conversion and Synthetic BPA Knee

**Source line:** 465.

**Highlighted text:**

```text
At the time of writing the conversion stage is complete and verified
```

**Comment — ajh26 — 14 October, 4:16 pm:**

```text
The wording of this is awkward. Of course you can only give results 'at the time of writing'. This should be true for everything in this dissertation...
```

#### C21 — ajh26 — 14 October, 4:13 pm

**Section:** Materials and Methods → Preliminary Neuromechanical Simulation Workflows → MuJoCo Conversion and Synthetic BPA Knee

**Source line:** 465.

**Highlighted text:**

```text
from simulation torques rather than static model tables, because coupled path points move
```

**Comment — ajh26 — 14 October, 4:13 pm:**

```text
I don't know what this means
```

#### C22 — ajh26 — 14 October, 4:18 pm

**Section:** Materials and Methods → Preliminary Neuromechanical Simulation Workflows → Simulink SNS Library and Reflex Topology

**Source line:** 483.

**Highlighted text:**

```text
settled on
```

**Comment — ajh26 — 14 October, 4:18 pm:**

```text
You need to go through this whole section and rewrite pretending you are writing for a journal.
```

#### C23 — ajh26 — 14 October, 4:22 pm

**Section:** Materials and Methods → Preliminary Neuromechanical Simulation Workflows → Cross-Platform Verification of the Contact-Driven Rhythm Core

**Source line:** 503.

**Highlighted text:**

```text
four persistent-sodium half-centers
```

**Comment — ajh26 — 14 October, 4:22 pm:**

```text
In general, when talking about your simulation work, you have not provided enough background for people to understand the computational neuroscience you are bringing up here. You either need to provide significant background in computational neuroscience and the models you are using (A LOT of work), or simplify your language to be more understandable by MEs, pointing people to the CompNeuro work to read more (less work).
```

#### C24 — ajh26 — 11 September, 9:14 am

**Section:** Original section unavailable; detached/end-of-file anchor.

**Source line:** 521 (current end-of-file anchor only; original line unavailable).

**Highlighted text:**

Unavailable: no surviving highlighted source text is displayed.

**Comment — ajh26 — 11 September, 9:14 am:**

```text
Stop conditions?
```

#### C25 — ajh26 — 11 September, 9:18 am

**Section:** Original section unavailable; detached/end-of-file anchor.

**Source line:** 521 (current end-of-file anchor only; original line unavailable).

**Highlighted text:**

Unavailable: no surviving highlighted source text is displayed.

**Comment — ajh26 — 11 September, 9:18 am:**

```text
I'm not sure I understand the point of this. Is this for the optimization later? Or are you modeling the compliance change in the attachment? You then next talk about optimized stiffness terms, but it isn't clear yet that those values were optimized.
```

#### C26 — ajh26 — 11 September, 9:21 am

**Section:** Original section unavailable; detached/end-of-file anchor.

**Source line:** 521 (current end-of-file anchor only; original line unavailable).

**Highlighted text:**

Unavailable: no surviving highlighted source text is displayed.

**Comment — ajh26 — 11 September, 9:21 am:**

```text
x and z are not independent. Make sure you discuss why. Perhaps set them up as independent, then discuss why they equal each other and make the new equation eliminated x3
```

#### C27 — ajh26 — 11 September, 9:22 am

**Section:** Original section unavailable; detached/end-of-file anchor.

**Source line:** 521 (current end-of-file anchor only; original line unavailable).

**Highlighted text:**

Unavailable: no surviving highlighted source text is displayed.

**Comment — ajh26 — 11 September, 9:22 am:**

```text
New section start here to talk about otpimization of parameters.
```

**Reply — ajh26 — 11 September, 9:22 am:**

```text
Make sure you do a good job justifying the use of optimization instead of analytical or theoretical methods.
```

#### C28 — ajh26 — 11 September, 9:27 am

**Section:** Original section unavailable; detached/end-of-file anchor.

**Source line:** 521 (current end-of-file anchor only; original line unavailable).

**Highlighted text:**

Unavailable: no surviving highlighted source text is displayed.

**Comment — ajh26 — 11 September, 9:27 am:**

```text
Some random thoughts I just had. Since muscles have similar force-length curves, do we not see a lot of large squishy creatures because of the loss in force? Also, not a lot of 'protruding' attachment points? Might be worth pointing out in the discussion after more rumination.
```

#### C29 — ajh26 — 11 September, 9:28 am

**Section:** Original section unavailable; detached/end-of-file anchor.

**Source line:** 521 (current end-of-file anchor only; original line unavailable).

**Highlighted text:**

Unavailable: no surviving highlighted source text is displayed.

**Comment — ajh26 — 11 September, 9:28 am:**

```text
This is confusing. Do you just mean because it wraps around and makes contact with the knee?
```

**Reply — ajh26 — 11 September, 9:29 am:**

```text
Also, are you starting a new subject here...new paragraph?
```

#### C30 — ajh26 — 11 September, 10:08 am

**Section:** Original section unavailable; detached/end-of-file anchor.

**Source line:** 521 (current end-of-file anchor only; original line unavailable).

**Highlighted text:**

Unavailable: no surviving highlighted source text is displayed.

**Comment — ajh26 — 11 September, 10:08 am:**

```text
This whole section about optimization is quite confusing and without much context.
```

#### C31 — ajh26 — 11 September, 10:10 am

**Section:** Original section unavailable; detached/end-of-file anchor.

**Source line:** 521 (current end-of-file anchor only; original line unavailable).

**Highlighted text:**

Unavailable: no surviving highlighted source text is displayed.

**Comment — ajh26 — 11 September, 10:10 am:**

```text
?
```

#### C32 — ajh26 — 11 September, 10:14 am

**Section:** Original section unavailable; detached/end-of-file anchor.

**Source line:** 521 (current end-of-file anchor only; original line unavailable).

**Highlighted text:**

Unavailable: no surviving highlighted source text is displayed.

**Comment — ajh26 — 11 September, 10:14 am:**

```text
?
```

#### C33 — ajh26 — 11 September, 10:23 am

**Section:** Original section unavailable; detached/end-of-file anchor.

**Source line:** 521 (current end-of-file anchor only; original line unavailable).

**Highlighted text:**

Unavailable: no surviving highlighted source text is displayed.

**Comment — ajh26 — 11 September, 10:23 am:**

```text
This is a much better way to introduce the optimization.
```

#### C34 — ajh26 — 11 September, 10:27 am

**Section:** Original section unavailable; detached/end-of-file anchor.

**Source line:** 521 (current end-of-file anchor only; original line unavailable).

**Highlighted text:**

Unavailable: no surviving highlighted source text is displayed.

**Comment — ajh26 — 11 September, 10:27 am:**

```text
Things in this section need to be explained more clearly and clearly justified.
```

#### C35 — ajh26 — 11 September, 10:28 am

**Section:** Original section unavailable; detached/end-of-file anchor.

**Source line:** 521 (current end-of-file anchor only; original line unavailable).

**Highlighted text:**

Unavailable: no surviving highlighted source text is displayed.

**Comment — ajh26 — 11 September, 10:28 am:**

```text
Should probably be in the background.
```

#### C36 — ajh26 — 11 September, 10:39 am

**Section:** Original section unavailable; detached/end-of-file anchor.

**Source line:** 521 (current end-of-file anchor only; original line unavailable).

**Highlighted text:**

Unavailable: no surviving highlighted source text is displayed.

**Comment — ajh26 — 11 September, 10:39 am:**

```text
These types of sentences can maybe be rephrased to stand alone better and be more like the structure of a background.
```

#### C37 — ajh26 — 11 September, 10:41 am

**Section:** Original section unavailable; detached/end-of-file anchor.

**Source line:** 521 (current end-of-file anchor only; original line unavailable).

**Highlighted text:**

Unavailable: no surviving highlighted source text is displayed.

**Comment — ajh26 — 11 September, 10:41 am:**

```text
Rewrite this section to have more of a tone that you developed a tool or a pipline for use. Get rid of the word preliminary.
```

#### C38 — mcneal — 15 September, 4:52 pm

**Section:** Original section unavailable; detached/end-of-file anchor.

**Source line:** 521 (current end-of-file anchor only; original line unavailable).

**Highlighted text:**

Unavailable: no surviving highlighted source text is displayed.

**Comment — mcneal — 15 September, 4:52 pm:**

```text
Here and in this paragraph. Vector is more than direction.

The explanation seems to move away from Fletcher's, which causes a couple of small inaccuracies. 
```

#### C39 — mcneal — 15 September, 4:52 pm

**Section:** Original section unavailable; detached/end-of-file anchor.

**Source line:** 521 (current end-of-file anchor only; original line unavailable).

**Highlighted text:**

Unavailable: no surviving highlighted source text is displayed.

**Comment — mcneal — 15 September, 4:52 pm:**

```text
position vector
```

#### C40 — mcneal — 15 September, 4:56 pm

**Section:** Original section unavailable; detached/end-of-file anchor.

**Source line:** 521 (current end-of-file anchor only; original line unavailable).

**Highlighted text:**

Unavailable: no surviving highlighted source text is displayed.

**Comment — mcneal — 15 September, 4:56 pm:**

```text
Scalars cannot cross multiply, and the scalar value is probably best represented with a subscript to clarify it's axis. Alternatively, you can define it in the text so that the relationship is clear.
```

#### C41 — mcneal — 15 September, 4:59 pm

**Section:** Original section unavailable; detached/end-of-file anchor.

**Source line:** 521 (current end-of-file anchor only; original line unavailable).

**Highlighted text:**

Unavailable: no surviving highlighted source text is displayed.

**Comment — mcneal — 15 September, 4:59 pm:**

```text
Contradicts the previous text where the relationship was a simplification. Unclear how d and r are related.
```

#### C42 — mcneal — 15 September, 5:51 pm

**Section:** Original section unavailable; detached/end-of-file anchor.

**Source line:** 521 (current end-of-file anchor only; original line unavailable).

**Highlighted text:**

Unavailable: no surviving highlighted source text is displayed.

**Comment — mcneal — 15 September, 5:51 pm:**

```text
This relies on several conditions that should probably be stated.  Off the top of my head, I can think of a few ways this would be an unsafe generalization:

- The bracket and tendon do not carry the same force.
- Their displacements aren't kinematically additive.
- The bracket deflection changes \(\hat{\mathbf u}\), so the line of action rotates.
- Rotational compliance or force/moment coupling matters.
- Either element is nonlinear.
- Tendon preload makes its tangent stiffness state-dependent.
- The bracket displacement can increase rather than decrease actuator length depending on geometry/sign convention. 

Those are the conditions I would consider explicitly stating to be safe. 
```

#### C43 — mcneal — 15 September, 5:55 pm

**Section:** Original section unavailable; detached/end-of-file anchor.

**Source line:** 521 (current end-of-file anchor only; original line unavailable).

**Highlighted text:**

Unavailable: no surviving highlighted source text is displayed.

**Comment — mcneal — 15 September, 5:55 pm:**

```text
Unclear whether this was calculated without the tendon or if its the sum of both in series.
```

### 30-results.tex

#### C44 — ajh26 — 14 October, 3:56 pm

**Section:** Results → Balance Platform and Inverted-Pendulum Demonstrations

**Source line:** 378.

**Highlighted text:**

```text
\section{Balance Platform and Inverted-Pendulum Demonstrations}\label{sec:results_balance}
```

**Comment — ajh26 — 14 October, 3:56 pm:**

```text
If these are the results, definitely remove the balance platform/inverted pendulum things.
```

#### C45 — bbolen83 — 18 September, 12:28 am

**Section:** Original section unavailable; detached/end-of-file anchor.

**Source line:** 599 (current end-of-file anchor only; original line unavailable).

**Highlighted text:**

Unavailable: no surviving highlighted source text is displayed.

**Comment — bbolen83 — 18 September, 12:28 am:**

```text
Put it on a page with table 4.1 and rotate them
```

### 94-AppendixC.tex

#### C46 — bbolen83 — 18 September, 12:26 am

**Section:** Original section unavailable; detached/end-of-file anchor.

**Source line:** 249 (current end-of-file anchor only; original line unavailable).

**Highlighted text:**

Unavailable: no surviving highlighted source text is displayed.

**Comment — bbolen83 — 18 September, 12:26 am:**

```text
No no, this has one particular knee angle but this should, in fact, be a more general T_k = [cos theta_k -sin theta_k.  etc
```

## Pending tracked changes

These entries remain pending; no change was accepted or rejected.

### T01 — 10-introduction.tex — Deleted

**Author:** bbolen83 — 27 September, 5:40 am

**Section:** Original section(s) unavailable; grouped at the current end of file.

**Source line:** 50 (current end-of-file anchor only; original fragment locations unavailable).

**Deleted text (full review entry):**

```text


\section{Motivation}

Biomimetic robots provide a platform for investigating how mechanical structure and control interact to generate movement, offering complementary insights to simulation-based studies \citep{andrews_general_1985, millard_flexing_2013, arnold_model_2010}. Wheeled robots remain the platform of choice in structured environments, where their speed, efficiency, and payload capacity are difficult to match; in the unstructured terrain of the natural and built world, legged machines hold the advantage, and biomimetic legs that reproduce human force capabilities have the potential to excel where wheels cannot \citep{frey_legged_2026, ha_learning_2025}. Musculoskeletal humanoids and limb systems driven by tendon-like actuators have been shown to support biologically meaningful experiments on motion generation, internal force transmission, and neural control hypotheses \citep{asano_design_2017, asano_musculoskeletal_2019, hitzmann_anthropomorphic_2018, liu_robotic_2018, chen_realizing_2019}. The dual perspective of using robots to learn about biology, and using biology to design more adaptive robots, has been emphasized in both biomimetic robotics and biomechanics \citep{ijspeert_amphibious_2020, asano_design_2017}. With such robotic systems, experiments can be performed that would not be practical or ethical with human test subjects, such as removing a specific muscle from the system or blocking a sensory feedback pathway. Experiments of this kind allow direct testing of how the components of the neuromechanical system interact, and they offer new opportunities for testing neural control theories that are faster and cheaper than human observation studies \citep{shin_constructive_2018, goldsmith_investigating_2021, schilling_adaptive_2022}.

Artificial muscles are a key enabling technology for this class of robot. Among the available actuation technologies, McKibben-style braided pneumatic actuators (BPAs) combine low mass, high force-to-weight ratio, and intrinsic compliance with a force--length curve that qualitatively resembles that of skeletal muscle \citep{liang_comparative_2020, hunt_modeling_2017}. BPAs have been used in quadrupedal robots \citep{aschenbeck_design_2006, hunt_modeling_2017, scharzenberger_design_2019}, musculoskeletal humanoids \citep{asano_design_2017, hitzmann_anthropomorphic_2018}, and dexterous manipulation systems \citep{chen_realizing_2019}. However, BPAs also differ from biological muscle in important ways: their maximum tensile force occurs at resting length rather than at an optimal fiber length, their maximum contraction ratio is substantially lower than that of human muscle, and their force output is a highly nonlinear function of internal pressure, contraction, and loading state \citep{hunt_modeling_2017, sarosi_comparative_2017, liang_comparative_2020}. Designing a biomimetic robot around BPAs therefore requires quantitative, predictive models of the actuator and of the joint torques it can produce.

\section{The Biomimetic Humanoid Robot Project}\label{sec:project}

This dissertation is part of a long-running effort at Portland State University to build a biomimetic humanoid robot from the trunk down, actuated by pneumatic artificial muscles and controlled by a synthetic neural network (SNN). The project's guiding hypothesis is that a robot with human-like biomechanics --- human-like muscle numbers, architecture, and routing, human-like joint ranges of motion, and human-like isometric torque profiles --- controlled by a biologically grounded neural architecture, will be capable of biomimetic movement and will serve as a testbed for neuromechanical hypotheses, including loss-of-function experiments that are impossible in living animals. A by-product of this work is expected to be improved prosthetic and orthotic device design.

The present work builds directly on a series of studies by the author and collaborators. In the author's master's project, the kinematics and dynamics of a biomimetic robot leg were matched against the Gait2392 OpenSim model, artificial muscle attachment and wrapping locations were determined for the robot legs, and MATLAB tooling was developed for evaluating muscle lengths, moment arms, and isometric joint torques over the full range of motion \citep{bolen_bipedal_2019, bolen_determination_2019}. That work showed that simply transplanting human muscle attachment locations onto a BPA-actuated robot does not reproduce human torque profiles. Hand placement of the muscle origin and insertion locations for 27 muscle-actuators per side was required. In parallel, a biomimetic four-bar knee linkage with a migrating instantaneous center of rotation was designed and validated by Steele \citep{steele_development_2017, steele_biomimetic_2018, steele_biomimetic_2018-1}, and a canine-inspired quadruped platform (Muscle Mutt, formerly  called DoggyDeux) with antagonistic BPA pairs and synthetic neural network control was developed by \citet{scharzenberger_design_2019}, extending an earlier demonstration by the Agile and Adaptive Robotics Laboratory (AARL) of central pattern generator (CPG) control of a pneumatically actuated robot \citep{hunt_development_2017}.


\section{Problem Statement}\label{sec:problem}

Within this context, three gaps stand between the existing robot design and a robot that can walk under neural control.

First, the actuator force models available at the start of this work were inadequate for design purposes. The manufacturer's simulation tool over-predicts the force of short BPAs by large margins, and the published empirical models are fitted to a single resting length each, so they do not accurately predict how force changes when a designer chooses a different actuator length \citep{sarosi_comparative_2017, martens_modeling_2017, hunt_modeling_2017}. Because muscle resting length is one of the primary degrees of freedom available to the robot designer, this knowledge-gap reduces confidence that a given design will perform as expected once it is built. 

Second, actuator force alone does not determine joint torque. The robot knee applies BPA force through brackets, artificial tendons, and wrapping paths. It was unknown how much these effects reduce usable torque relative to the rigid-body prediction, or how to best correct for them. Addtionally, the biomimetic knee's optimized muscle routing \citep{morrow_optimization_2020} had never been validated experimentally or compared against human isometric torque data.

Third, no simulation pipeline existed for closing the loop between a musculoskeletal robot model of this kind and a biologically realistic neural controller. The Neuromechanical Modeling literature provides the ingredients --- two-level spinal CPG architectures \citep{rybak_modelling_2006, mccrea_organization_2008}, neuromechanical simulation frameworks \citep{markin_afferent_2010}, and synthetic nervous system tooling \citep{nourse_sns_2023} --- but they had not been assembled for a BPA-driven bipedal robot with human-like routing.

\section{Research Objectives}\label{sec:objectives}

This dissertation addresses the first two gaps experimentally and lays the groundwork for the third. The objectives are:

\begin{enumerate}
    \item \textbf{Characterize the isometric force of Festo BPAs} with uninflated diameters of \qtylist{10;20}{\mm} over a wide range of resting lengths, contractions, and pressures; derive improved maximum-force and normalized force--strain--pressure models; and compare them against existing models (Chapters~\ref{ch:methods} and~\ref{ch:results}, Section~\ref{sec:results_force}).
    
    \item \textbf{Validate and improve the prediction of isometric knee torque} produced by BPAs on a pinned joint and on a biomimetic four-bar knee; identify and model the losses due to constant length offsets, series compliance, and actuator wrapping; and compare the resulting torque envelopes with human capabilities (Methods Sections~\ref{sec}--\ref{sec}, Results Section~\ref{sec}, and Discussion Section~\ref{sec})

    \item \textbf{Develop tools and a framework for neuromechanical control} of the humanoid robot by converting its OpenSim musculoskeletal model to MuJoCo, implementing an SNS-Toolbox two-level spinal CPG with proprioceptive and heel/toe contact feedback, and establishing a staged tuning and verification workflow (Methods Section~\ref{sec}, Results Section~\ref{sec}, and Discussion Section~\ref{sec}).

\end{enumerate}

\section{Dissertation Organization}\label{sec:organization}

Chapter~\ref{ch:background} reviews the relevant background: biomimetic and musculoskeletal robots, artificial muscle technologies and BPA modeling, human lower-limb biomechanics and the OpenSim benchmark models, the neural control of locomotion with an emphasis on two-level CPG architectures and sensory afferent feedback, synthetic nervous systems and the Sensory Afferent Database, and the simulation tooling that connects OpenSim models to MuJoCo. 

Chapter~\ref{ch:methods} presents the materials and methods: the actuator force model and its characterization experiment, the pinned-knee and biomimetic knee test stands, the hybrid torque calculation method that isolates the loss mechanisms, the correction terms fitted to them, and the optimization algorithms used to redesign actuator placements. 

Chapter~\ref{ch:results} first summarizes the completed foundation studies (the placement optimization, the Sensory Afferent Database, the BPA balance platform, and the BPA pressure-dynamics model), then reports the actuator characterization, the identification of the torque correction terms on the pinned and biomimetic knees, the redesigned actuator routes, and the experimental comparison of the biomimetic knee with human isometric torque. 

Chapter~\ref{ch:discussion} discusses the findings, their limitations, and their implications for real-time control. 

Chapter~\ref{ch:futurework} assembles the validated model and the neural controller into a working neuromechanical simulation pipeline, an OpenSim-to-MuJoCo conversion coupled to a two-layer CPG with full sensory feedback; with these tools in place, researchers can run lesion, stimulation, and gait-modulation experiments on the simulated humanoid leg today, and the same tools define the path to a treadmill robot with embedded real-time control. 

Chapter~\ref{ch:conclusion} concludes the report.

Appendix~\ref{app:muscle_calcs} details the specifics of calculating calculating homogenous transformation matrices, BPA strain, relative strain, force vectors, displacement vectors, moment arm calculations, and torque about a joint.

Appendix~\ref{app:neuron_equations} details neuron model equations used by the tools.

Appendix~\ref{app:model_geometry} gives geometry of the models used and paths of BPAs and muscle-actuator
```

### T02 — 15-background.tex — Changed

**Author:** ajh26 — 14 October, 1:31 pm

**Section:** Background → Braided Pneumatic Actuators

**Source line:** 38.

**Original text:**

```text
requirement
```

**Replacement text:**

```text
improvement
```

### T03 — 15-background.tex — Changed

**Author:** ajh26 — 14 October, 2:05 pm

**Section:** Background → Simulation Tooling: From OpenSim to MuJoCo

**Source line:** 158.

**Original text:**

```text
these tools were 
```

**Replacement text:**

```text
I 
```

### T04 — 15-background.tex — Added

**Author:** ajh26 — 14 October, 2:05 pm

**Section:** Background → Simulation Tooling: From OpenSim to MuJoCo

**Source line:** 158.

**Added text (full review entry):**

```text
these tools 
```

### T05 — 15-background.tex — Added

**Author:** ajh26 — 14 October, 2:48 pm

**Section:** Background → Testing Neuromechanical Hypotheses on Robots: The Ethical Case

**Source line:** 170.

**Added text (full review entry):**

```text
also 
```

### T06 — 20-methods.tex — Deleted

**Author:** ajh26 — 14 October, 3:13 pm

**Section:** Materials and Methods → Balance Platform and Inverted-Pendulum Testbeds

**Source line:** 389.

**Deleted text (full review entry):**

```text
 walking
```

### T07 — 20-methods.tex — Changed

**Author:** ajh26 — 14 October, 4:02 pm

**Section:** Materials and Methods → Prior Neuromechanical Simulation: The AnimatLab Bipedal Walker

**Source line:** 410.

**Original text:**

```text
increasing foot compliance and conserving energy
```

**Replacement text:**

```text
simplifying the system
```

### T08 — 20-methods.tex — Deleted

**Author:** bbolen83 — 27 September, 3:09 pm

**Section:** Original section(s) unavailable; grouped at the current end of file.

**Source line:** 521 (current end-of-file anchor only; original fragment locations unavailable).

**Deleted text (full review entry):**

```text
01_testJigs1.pdf\includegraphics[width=0.75\textwidth]{testJigs.pdf}KneeICR.epsThe bracket stiffness matrix is modeled as a diagonal matrix to reflect independent stiffness components along its principal axes:_{br}'^^\begin{figure}[htbp]
\begin{center}
\includegraphics[width=1\textwidth]{bktFrame.pdf}
\end{center}
\caption{Bracket frame orientations for: \textbf{(A)} Flexor configuration on the pinned knee, \textbf{(B)} extensor configuration on the pinned knee, \textbf{(C)} flexor and extensor configurations on the biomimetic knee with $\phi$\qty{10}{\mm} BPAs, and \textbf{(D)} flexor configuration on the biomimetic knee with the $\phi$\qty{20}{\mm} BPA.}
\label{fig:bktFrame}
\end{figure} only$K_$K_ongoing}.

The goodness-of-fit metrics used to judge the optimization results are the root mean square error (RMSE), the fraction of variance unexplained (FVU), and the maximum absolute residual, each computed between the measured and predicted isometric torque over the tested range of motion.count.

\section{Prior Neuromechanical Simulation: The AnimatLab Bipedal Walker}\label{sec:prior_walker}The course report describing that model \citep{bolen_design_2020} is summarized here because its biomechanics and network architecture are the direct ancestors of both the current AnimatLab walker model, whose neuron and synapse parameters are documented in Appendix~\ref{app:neuron_equations}, and the two-layer CPG controller planned in Chapter~\ref{ch:futurework}.plant; its contribution to push-off would make an interesting study in its own right once the model walks on the ground, so lateral stabilization falls to the vestibular layer and hardware design of Chapter~\ref{ch:futurework}plantplantplant
```

### T09 — 30-results.tex — Deleted

**Author:** bbolen83 — 27 September, 5:35 am

**Section:** Original section(s) unavailable; grouped at the current end of file.

**Source line:** 599 (current end-of-file anchor only; original fragment locations unavailable).

**Deleted text (full review entry):**

```text
\chapter{Results}\label{ch:results}\par}

\qty{620}{\kPa} .$\chi_{0}$, $\chi_{1}$, and $\chi_{1}$\textbf{B)}acheivable\section{Preliminary Neuromechanical Simulation Recordings}\label{sec:preliminary_sim_results}

The AnimatLab phase-1 run described in Section~\ref{sec:preliminary_sim_methods} produced synchronized joint, neural, contact-sensory, and body-motion recordings over \qty{10}{\second}. Basic checks confirmed finite values and a uniform \qty{0.2}{\milli\second} sampling interval in both exported charts. Each chart contained 50,001 samples from 0 through \qty{10}{\second}; nine trailing zero-filled rows beyond the configured chart end were excluded. Figure~\ref{fig:animatlab_phase1_preliminary} shows bilateral hip, knee, and ankle angles. Neural and contact-sensory recordings remain archived with the run; the contact-sensory channels record membrane voltage rather than contact force. These recordings demonstrate signal acquisition and alternating activity in this run, rather than validated locomotor performance. The model included an externally applied horizontal force during the first second, and neither stable walking nor effective ground-contact gating was established.

\begin{figure}[htbp]
\centering
\includegraphics[width=\textwidth,height=0.60\textheight,keepaspectratio]{figs/Preliminary/animatlab_phase1_joint_angles.pdf}
\caption[Preliminary AnimatLab hip, knee, and ankle recordings]{Preliminary AnimatLab phase-1 joint motion: \textbf{(A)} hip, \textbf{(B)} knee, and \textbf{(C)} ankle angles in the model's recorded joint conventions. Blue solid and orange dashed traces denote left and right sides. All 50,001 samples within the configured 0--\qty{10}{\second} chart window are retained without filtering. The model includes an external horizontal force during the first second. These traces document joint motion in an instrumented run; they do not establish stable walking or effective ground-contact gating.}
\label{fig:animatlab_phase1_preliminary}
\end{figure}

\section{Follow-Up Identification and Route Redesign}\label{sec:ongoing}

The identification work outlined at the close of the experimental chapters has since been completed, and its outcomes are summarized here; the coefficient values shown with the figures above (Table~\ref{tab:OptCoeff}) belong to the initial single-bracket identification of those experiments.

The correction terms were re-identified over the full archived test sets with the two-stage procedure of Section~\ref{sec:improved_model}: a two-bracket flexor fit cross-validated over the five pinned-flexor tests (the encoder-corrected \qty{47}{\cm} test held out), followed by an extensor refit with the stiffness pair locked to the flexor front. The adopted values are $\chi_{0} = \qty{+8.9}{\mm}$, $\chi_{1} = \qty{5.62e4}{\N\per\m}$, and $\chi_{2} = \qty{1.85e4}{\N\per\m}$ for the flexor, and $\chi_{0} = \qty{-10.1}{\mm}$ with $\chi_{3} = 0.621$ for the extensor on the locked pair ($4.354\times10^{4}$, $1.701\times10^{4}$~N/m); full provenance is given in Appendix~\ref{app:model_geometry}, Table~\ref{tab:app_xi_values}. At these values, the corrected flexor model predicts the five pinned-flexor tests with per-test RMSE of 1.94/1.51/2.29/1.63/1.21~N{\,}m (48, 46, 47, 40-tendon, 41\,cm; FVU at most 0.17), the \qty{47}{\cm} test being pure held-out validation. On the six-test extensor validation pool, per-test RMSE spans 0.56--1.88~N{\,}m against 0.95--3.62~N{\,}m for the uncorrected rigid model --- an improvement on every pool test. On the biomimetic extensor test ($l_{0} = \qty{51.8}{\cm}$) the corrected prediction improves substantially on the baseline (RMSE 2.20 vs.\ 3.03~N{\,}m) but leaves FVU above unity, so the shared-pair model still under-predicts that configuration's measured torque variation; this remains the open gap in the cross-configuration validation. The solver question raised during the campaign was settled empirically: the genetic multiobjective solver was retained for the identification, where the Pareto front is required for holdout validation, and the surrogate-plus-pattern-search pipeline was reserved for the route redesign where it is already the method of record (Section~\ref{sec:optimization}).

The route redesign built on these values is likewise far advanced. The two-BPA flexor arrangement --- motivated by the single-muscle shortfall at high flexion angles documented above --- reached a feasible solution meeting the flexor torque target under the modeling assumption that both actuators share one attachment route rather than a physically mirrored pair; making the mirrored pair explicit is the next design refinement, and the flexor re-run under the newly adopted correction terms was still in progress at the time of writing. The extensor route was redesigned directly under the adopted correction terms: it meets or exceeds its human-scaled target at every angle in the range of motion, with the thinnest torque margin below $0.1\%$ at a single angle, while respecting the wrap-release schedule and contraction limits of the route-elimination rule (Appendix~\ref{app:model_geometry}). These redesign solutions are the plant definitions carried into the neuromechanical simulation work of Sections~\ref{sec:preliminary_sim_methods}--\ref{sec:preliminary_sim_results} and Chapter~\ref{ch:futurework}.


```

### T10 — 40-discussion.tex — Deleted

**Author:** bbolen83 — 27 September, 3:40 pm

**Section:** Original section(s) unavailable; grouped at the current end of file.

**Source line:** 171 (current end-of-file anchor only; original fragment locations unavailable).

**Deleted text (full review entry):**

```text
it is also with theWhere our measured data differs significantly from the Festo model, is that as\N} Therefore, the errorfactors such as:The reason for BPA models predicting higher force than what is measured experimentally is a recurrent theme across many models: they do not account for the end effects of the pressure vessel where end effects dominate.Near the boundary conditions are where axial, radial, and circumferential strains deviate from analytic solutions to a cylindrical model.play an important factorModeling work of ideal McKibben actuators describeenable for fasterThe diameter adjustment can be impliedwill benefit more from improved modeling of airflow dynamics than dynamics involving actuator hysteresis, damping, and other dynamic behavior.For example, the describes the passiveHowever, there are some distinct differences between how these artificial muscles behave when compared to biological muscles.For example, the biological muscle optimum fiber length $l_{OFL}$ is not equivalent to $l_0$ in the muscles. The maximum active force that the artificial muscles can produce is at the resting length, also the maximum length, of the artificial muscles.the maximum active force a biological muscle length can produceposition location toleranceSimplifying the model and testing it allowed us to see how we were deficient in our previous analysis.The bracket frame of reference was placed where it starts to cantilever and is rotated about $Z$.Optimization found results forvalidated with secondDuring greater magnitude knee flexion angles, shorter $l_{0}$ BPAs were noted as being stretched and compressed (see Fig.~\ref{fig:testJigs}\textbf{C}).The same force that caused the stretch and compression would also sometimes cause the BPA to slip off of the anterior bolt head of the bolt that holds the femur to the knee.This would result in a $\pm \; 20\,\text{mm}$ displacement in the $Z$ direction, which means a change the force vector and a change in torque.However, the  hybrid calculated torque matching or exceeding predicted torque, except for higher degrees of flexion in the $l_{0}=\qty{41.5}{\cm}$ BPA, indicates a loss of functional $l_m$; this is why the $\chi_{3}$ term was created in \eqref{eq:deltaL}.The predicted line more closely matches the torque values. However, the predicted torque is still higher than experimentally measured torque.In the presented nonlinear biomimetic robotic system, rigid body mechanics simplifications, often encountered in other robotic and biomechanical systems, are removed.Instead, working from the simple to the complex elucidated mechanisms of system compliance that will lead to much more accurate first iteration designs.plant model
```

### T11 — 50-futurework.tex — Deleted

**Author:** bbolen83 — 27 September, 5:30 am

**Section:** Original section(s) unavailable; grouped at the current end of file.

**Source line:** 59 (current end-of-file anchor only; original fragment locations unavailable).

**Deleted text (full review entry):**

```text
characterizes the plant:mechanical plantThe suspended walker is the CPG acting as a feedforward controller:That premise is itself well supported: biomechanics and sensory feedback have been shown to regularize the behavior of different locomotor CPG architectures \citep{deng_biomechanical_2022}, so the contribution of the present work lies in applying this principle to the suspended walker rather than in the principle.The AnimatLab walker model, the direct descendantfollowing the models of left--right coordinationThree classes of proprioceptive feedback will close the loop between the MuJoCo plant and the CPG,in which tonic stimulus is appliedoscillate at an approximately \qty{0.9}{\second} periodvisual measurements of the horizon and surrounding(e.g., walking, running, turning)on a embedded computing stackcontact sensing), and an Nvidia Jetson-class companion computere.g., run, walk, stimulate a named neural pathwayfeedforward plant model
```

### T12 — 60-conclusion.tex — Deleted

**Author:** bbolen83 — 28 September, 9:30 am

**Section:** Original section(s) unavailable; grouped at the current end of file.

**Source line:** 12 (current end-of-file anchor only; original fragment locations unavailable).

**Deleted text (full review entry):**

```text
\qtylist{10;20}{\mm} diameter BPAspredicts with the wrong trendselect actuator diameter and resting length, and to evaluatea constant length offset, series compliance represented by an equivalent stiffness, and a wrapping loss termis realizable in hardware, and identifyingMuJoCo plant convertedneural lesions, sensory deletions, pathway stimulation
```

### T13 — 92-AppendixA.tex — Deleted

**Author:** bbolen83 — 15 September, 4:22 pm

**Section:** Original section(s) unavailable; grouped at the current end of file.

**Source line:** 249 (current end-of-file anchor only; original fragment locations unavailable).

**Deleted text (full review entry):**

```text
, as they exist in the project's MATLAB codeFcompressive
```

### T14 — 94-AppendixC.tex — Deleted

**Author:** bbolen83 — 27 September, 3:31 pm

**Section:** Original section(s) unavailable; grouped at the current end of file.

**Source line:** 249 (current end-of-file anchor only; original fragment locations unavailable).

**Deleted text (full review entry):**

```text
Regenerating the results files updates this table.
```

### T15 — main.tex — Deleted

**Author:** bbolen83 — 14 September, 1:44 am

**Section:** Original section(s) unavailable; grouped at the current end of file.

**Source line:** 192 (current end-of-file anchor only; original fragment locations unavailable).

**Deleted text (full review entry):**

```text
Creatingwith
```

## Resolved comment archive

These threads were already resolved when inspected. Their status was not changed. Historical source line numbers and section anchors are not available in this archive.

### R01 — 20-methods.tex

**Section:** Historical section unavailable; see file above.

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

Unavailable: the archive displays an empty quoted-text field.

**Comment — ajh26 — 10 September, 2:24 pm:**

```text
Also used to design attachment locations for a robotic knee that meets the torque capabilities of a human knee.
```

**Resolution record — bbolen83 — 3 October, 12:53 am:**

```text
Marked as resolved
```

### R02 — 20-methods.tex

**Section:** Historical section unavailable; see file above.

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

Unavailable: the archive displays an empty quoted-text field.

**Comment — ajh26 — 11 September, 9:12 am:**

```text
With a quick read, I believe all this. Have Stu double check the equations.
```

**Resolution record — bbolen83 — 3 October, 12:53 am:**

```text
Marked as resolved
```

### R03 — 15-background.tex

**Section:** Historical section unavailable; see file above.

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

Unavailable: the archive displays an empty quoted-text field.

**Comment — ajh26 — 10 September, 1:43 pm:**

```text
it is not clear who has done this advancing
```

**Resolution record — bbolen83 — 3 October, 12:52 am:**

```text
Marked as resolved
```

### R04 — 15-background.tex

**Section:** Historical section unavailable; see file above.

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

Unavailable: the archive displays an empty quoted-text field.

**Comment — ajh26 — 10 September, 1:43 pm:**

```text
model?
```

**Resolution record — bbolen83 — 3 October, 12:52 am:**

```text
Marked as resolved
```

### R05 — 15-background.tex

**Section:** Historical section unavailable; see file above.

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

Unavailable: the archive displays an empty quoted-text field.

**Comment — ajh26 — 10 September, 1:45 pm:**

```text
Most of this information is elsewhere in your dissertation, and I don't really like this section. What isn't elsewhere can probably moved into the areas that talk about related items.
```

**Resolution record — bbolen83 — 3 October, 12:52 am:**

```text
Marked as resolved
```

### R06 — 15-background.tex

**Section:** Historical section unavailable; see file above.

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

Unavailable: the archive displays an empty quoted-text field.

**Comment — ajh26 — 10 September, 1:39 pm:**

```text
How about the rat stuff from myself, Deng, Young, and Jackson?
```

**Resolution record — bbolen83 — 3 October, 12:52 am:**

```text
Marked as resolved
```

### R07 — 15-background.tex

**Section:** Historical section unavailable; see file above.

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

Unavailable: the archive displays an empty quoted-text field.

**Comment — ajh26 — 10 September, 1:38 pm:**

```text
not the AARL
```

**Resolution record — bbolen83 — 3 October, 12:52 am:**

```text
Marked as resolved
```

### R08 — 10-introduction.tex

**Section:** Historical section unavailable; see file above.

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

Unavailable: the archive displays an empty quoted-text field.

**Comment — ajh26 — 10 September, 11:39 am:**

```text
Write a whole paragraph for each chapter. Could probably be combined with the Research Objectives section.
```

**Resolution record — bbolen83 — 2 October, 12:37 pm:**

```text
Marked as resolved
```

### R09 — 10-introduction.tex

**Section:** Historical section unavailable; see file above.

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

Unavailable: the archive displays an empty quoted-text field.

**Comment — ajh26 — 10 September, 11:37 am:**

```text
This is a whole paper, seems weird it is just a section while the previous paper is a chapter and a half.
```

**Resolution record — bbolen83 — 25 September, 2:27 am:**

```text
Marked as resolved
```

### R10 — 10-introduction.tex

**Section:** Historical section unavailable; see file above.

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

Unavailable: the archive displays an empty quoted-text field.

**Comment — ajh26 — 10 September, 11:32 am:**

```text
reword
```

**Resolution record — bbolen83 — 25 September, 2:12 am:**

```text
Marked as resolved
```

### R11 — 02-abstract.tex

**Section:** Abstract (file-level identification).

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

```text
manufacturer's tool predicts with the wrong trend
```

**Comment — bbolen83 — 14 September, 2:01 am:**

```text
Should I cut this or make it more vague?
```

**Resolution record — bbolen83 — 25 September, 2:03 am:**

```text
Marked as resolved
```

### R12 — 02-abstract.tex

**Section:** Abstract (file-level identification).

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

```text
\chapter*{\centerline{Abstract}}
\markboth{\MakeUppercase{Abstract}}{}
\iftoggle{fulltoc}{
  \addcontentsline{toc}{chapter}{Abstract}
}{}
```

**Comment — ajh26 — 10 September, 11:19 am:**

```text
Pull your abstract down to one page. 'Extra' material can be moved to the introduction and 'dissertation organization' section.
```

**Reply — bbolen83 — 14 September, 2:55 pm:**

```text
Ok. This was harder than I thought with the double-space.
```

**Resolution record — bbolen83 — 25 September, 2:03 am:**

```text
Marked as resolved
```

### R13 — 20-methods.tex

**Section:** Historical section unavailable; see file above.

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

Unavailable: the archive displays an empty quoted-text field.

**Comment — mcneal — 15 September, 5:25 pm:**

```text
Unclear if this is within the bracket frame or the bracket stiffness frame. Also, this implies that the deflection is purely translational. If that is true and there is no wrench, it might make sense to state that. 
```

**Resolution record — bbolen83 — 17 September, 11:55 pm:**

```text
Marked as resolved
```

### R14 — 20-methods.tex

**Section:** Historical section unavailable; see file above.

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

Unavailable: the archive displays an empty quoted-text field.

**Comment — mcneal — 15 September, 5:16 pm:**

```text
Feels like a rotation matrix is missing. It's unclear how the relationship shown connects the two frames. 
```

**Resolution record — bbolen83 — 17 September, 11:54 pm:**

```text
Marked as resolved
```

### R15 — 20-methods.tex

**Section:** Historical section unavailable; see file above.

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

Unavailable: the archive displays an empty quoted-text field.

**Comment — mcneal — 15 September, 5:25 pm:**

```text
Hooke's law
```

**Resolution record — bbolen83 — 17 September, 11:51 pm:**

```text
Marked as resolved
```

### R16 — 20-methods.tex

**Section:** Historical section unavailable; see file above.

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

Unavailable: the archive displays an empty quoted-text field.

**Comment — mcneal — 15 September, 5:23 pm:**

```text
F'
```

**Resolution record — bbolen83 — 17 September, 11:49 pm:**

```text
Marked as resolved
```

### R17 — 20-methods.tex

**Section:** Historical section unavailable; see file above.

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

Unavailable: the archive displays an empty quoted-text field.

**Comment — ajh26 — 11 September, 9:06 am:**

```text
I'm not sure what this means
```

**Reply — ajh26 — 11 September, 9:07 am:**

```text
OK, I think I get it now, but it should be made a little more clear where this value comes from.
```

**Resolution record — bbolen83 — 17 September, 11:47 pm:**

```text
Marked as resolved
```

### R18 — 20-methods.tex

**Section:** Historical section unavailable; see file above.

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

Unavailable: the archive displays an empty quoted-text field.

**Comment — mcneal — 15 September, 5:01 pm:**

```text
unit vector is not the direction
```

**Resolution record — bbolen83 — 17 September, 11:45 pm:**

```text
Marked as resolved
```

### R19 — 20-methods.tex

**Section:** Historical section unavailable; see file above.

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

Unavailable: the archive displays an empty quoted-text field.

**Comment — bbolen83 — 15 September, 1:20 pm:**

```text
Why not unit vector? -Stu
```

**Resolution record — bbolen83 — 17 September, 11:40 pm:**

```text
Marked as resolved
```

### R20 — 15-background.tex

**Section:** Historical section unavailable; see file above.

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

Unavailable: the archive displays an empty quoted-text field.

**Comment — ajh26 — 10 September, 1:50 pm:**

```text
I think you should pull back a little on how you refer to the 'future work'. 2 reasons. 1) YOu need to not confuse your committee on what is completed vs what is not. 2) You need to frame your dissertation as providing real benefits as is. The benefits include the cool experiments that will be performed in the future. The way this sentence is written kind of muddies the water.
```

**Resolution record — bbolen83 — 14 September, 3:50 pm:**

```text
Marked as resolved
```

### R21 — 15-background.tex

**Section:** Historical section unavailable; see file above.

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

Unavailable: the archive displays an empty quoted-text field.

**Comment — ajh26 — 10 September, 11:49 am:**

```text
Again, at the end, connect to your work.
```

**Resolution record — bbolen83 — 14 September, 3:33 pm:**

```text
Marked as resolved
```

### R22 — 15-background.tex

**Section:** Historical section unavailable; see file above.

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

Unavailable: the archive displays an empty quoted-text field.

**Comment — ajh26 — 10 September, 11:48 am:**

```text
plural?
```

**Resolution record — bbolen83 — 14 September, 3:33 pm:**

```text
Marked as resolved
```

### R23 — 15-background.tex

**Section:** Historical section unavailable; see file above.

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

Unavailable: the archive displays an empty quoted-text field.

**Comment — ajh26 — 10 September, 11:49 am:**

```text
Doesn't follow from other sentences
```

**Resolution record — bbolen83 — 14 September, 3:33 pm:**

```text
Marked as resolved
```

### R24 — 15-background.tex

**Section:** Historical section unavailable; see file above.

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

Unavailable: the archive displays an empty quoted-text field.

**Comment — ajh26 — 10 September, 11:47 am:**

```text
Doesn't follow from topic sentence/previous sentences.
```

**Resolution record — bbolen83 — 14 September, 3:33 pm:**

```text
Marked as resolved
```

### R25 — 10-introduction.tex

**Section:** Historical section unavailable; see file above.

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

Unavailable: the archive displays an empty quoted-text field.

**Comment — ajh26 — 10 September, 11:31 am:**

```text
This does not match with the previous paragraphs which is all preliminary work. This is work you did during the PhD, as such it should either be moved to a chapter (included in Dessertation Organization) or the discussions section.
```

**Reply — bbolen83 — 14 September, 3:07 pm:**

```text
Moved to chapter 30-results.
```

**Resolution record — bbolen83 — 14 September, 3:07 pm:**

```text
Marked as resolved
```

### R26 — 10-introduction.tex

**Section:** Historical section unavailable; see file above.

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

Unavailable: the archive displays an empty quoted-text field.

**Comment — ajh26 — 10 September, 11:24 am:**

```text
Put in a clear sentence here describing the gap(s) in knowledge.
```

**Resolution record — bbolen83 — 14 September, 2:56 pm:**

```text
Marked as resolved
```

### R27 — main.tex

**Section:** Historical section unavailable; see file above.

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

Unavailable: the archive displays an empty quoted-text field.

**Comment — ajh26 — 21 July, 3:53 pm:**

```text
To be removed.
```

**Resolution record — bbolen83 — 14 September, 1:44 am:**

```text
Marked as resolved
```

### R28 — main.tex

**Section:** Historical section unavailable; see file above.

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

Unavailable: the archive displays an empty quoted-text field.

**Comment — ajh26 — 21 July, 3:51 pm:**

```text
I think you should actually format this more like a big paper. Pull both your introductions into 1 chapter. Your Methods into 1 chapter. Your Results into 1 chapter. And your discussion into 1 chapter. A bit like the 'original' paper.
```

**Resolution record — bbolen83 — 14 September, 1:44 am:**

```text
Marked as resolved
```

### R29 — main.tex

**Section:** Historical section unavailable; see file above.

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

Unavailable: the archive displays an empty quoted-text field.

**Comment — ajh26 — 21 July, 3:46 pm:**

```text
You also need a Background chapter. This is where you flex your knowledge for your committee.
```

**Resolution record — bbolen83 — 14 September, 1:43 am:**

```text
Marked as resolved
```

### R30 — main.tex

**Section:** Historical section unavailable; see file above.

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

Unavailable: the archive displays an empty quoted-text field.

**Comment — ajh26 — 21 July, 3:46 pm:**

```text
It is fine to have a personal motivation section, but you need an introduction that motivates the work from a scientific standpoint. Take about 2-3 pages to build a good introduction motivation, and then create a paragraph summarizing each of your subsequent chapters.
```

**Resolution record — bbolen83 — 14 September, 1:43 am:**

```text
Marked as resolved
```

### R31 — 20-methods.tex

**Section:** Historical section unavailable; see file above.

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

Unavailable: the archive displays an empty quoted-text field.

**Comment — ajh26 — 11 September, 8:57 am:**

```text
Sometimes 2 muscles now?
```

**Resolution record — bbolen83 — 11 September, 12:35 pm:**

```text
Marked as resolved
```

### R32 — 20-methods.tex

**Section:** Historical section unavailable; see file above.

**Source line:** Unavailable in the resolved-comment archive.

**Saved highlighted / quoted text:**

Unavailable: the archive displays an empty quoted-text field.

**Comment — ajh26 — 11 September, 8:55 am:**

```text
Machined parts now too?
```

**Resolution record — bbolen83 — 11 September, 12:27 pm:**

```text
Marked as resolved
```
