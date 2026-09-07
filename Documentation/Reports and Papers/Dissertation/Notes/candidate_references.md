# Candidate references — verified 2026-09-07

All entries below were verified against publisher pages, PubMed, or DOI records on
2026-09-07 (no unverified citations). Each entry includes the exact insertion point
in the current monograph chapters. After you approve, I will: (1) add the entries to
`Bibliography/thesis.bib`, (2) edit the chapter files, (3) paste into Overleaf on your
go-ahead. You may also prefer to add them to Zotero yourself and re-export — both work.

---

## 1. Webb 2001 — robots as scientific models of biological behaviour

**Citation:** Webb, B. (2001). Can robots make good models of biological behaviour?
Behavioral and Brain Sciences, 24(6), 1033–1050. DOI: 10.1017/S0140525X01000127
(verified: PubMed PMID 12412325, Cambridge Core)

**Use:** The canonical methodological defense of biorobots as scientific models.
Directly supports the dissertation's framing claim that a robot can be a legitimate
testbed for neuromechanical hypotheses.

**Insert into:** Chapter 1, §1.1 Motivation, this sentence:
> "With such robotic systems, experiments can be performed that would not be
> practical or ethical with human test subjects, such as removing a specific muscle
> from the system or blocking a sensory feedback pathway."

becomes:
> "With such robotic systems, experiments can be performed that would not be
> practical or ethical with human test subjects, such as removing a specific muscle
> from the system or blocking a sensory feedback pathway
> \citep{webb_robots_2001, russell_principles_1959, tannenbaum_russell_2015}."

## 2. Webb 2002 — robots in invertebrate neuroscience

**Citation:** Webb, B. (2002). Robots in invertebrate neuroscience. Nature,
417(6886), 359–363. DOI: 10.1038/417359a (verified: Nature)

**Use:** Concrete precedent that neuroscientists accept robot models as carriers of
neural hypotheses (cricket phonotaxis, insect navigation). Pairs with Webb 2001.

**Insert into:** Chapter 2, §2.1 Biomimetic and Musculoskeletal Robots, end of the
first paragraph (after "...faster and cheaper than human observation studies" — or
the equivalent §2.1 sentence about neuromechanical experiments):
> "...such as removing a specific muscle from the system or blocking a sensory
> feedback pathway \citep{webb_robots_2001}. Robotic models have a track record of
> generating and testing neural hypotheses in neuroscience itself, from insect
> navigation to phonotaxis \citep{webb_robots_2002}."

## 3. Russell & Burch 1959 — the Three Rs (Replacement, Reduction, Refinement)

**Citation:** Russell, W. M. S., & Burch, R. L. (1959). The Principles of Humane
Experimental Technique. London: Universities Federation for Animal Welfare.
(reissued 1992; verified via Johns Hopkins CAAT online edition and Tannenbaum &
Bennett 2015)

**Use:** Origin of the 3Rs framework. The dissertation reframes a synthetic nervous
system robot as a *Replacement* platform: lesion, deletion, and stimulation
experiments that would require animals (or are impossible in humans) run at zero
moral cost.

**Insert into:** same sentence as Webb 2001 above, plus the ethics section
(see Notes/neuromechanics_state_of_art_and_ethics.tex).

## 4. Tannenbaum & Bennett 2015 — modern restatement of the 3Rs

**Citation:** Tannenbaum, J., & Bennett, B. T. (2015). Russell and Burch's 3Rs then
and now: The need for clarity in definition and purpose. Journal of the American
Association for Laboratory Animal Science, 54(2), 120–132.
(verified: PMC4382615)

**Use:** Authoritative modern reading of the 3Rs; strengthens the ethics argument
beyond the 1959 original.

**Insert into:** ethics passage (Ch. 1 §1.1 anchor sentence, and the full ethics
draft).

## 5. Morton & Bastian 2006 — cerebellum and locomotor adaptation

**Citation:** Morton, S. M., & Bastian, A. J. (2006). Cerebellar contributions to
locomotor adaptations during splitbelt treadmill walking. The Journal of
Neuroscience, 26(36), 9107–9116. DOI: 10.1523/JNEUROSCI.2617-06.2006
(verified: PubMed PMID 16957067)

**Use:** Experimental support for the cerebellar-layer claim in Future Work: the
cerebellum is essential for adapting locomotor patterns (split-belt walking),
which is exactly the "gateway between gaits" role proposed for the simulated
cerebellar layer.

**Insert into:** Chapter 6, §6.3, this sentence:
> "...and tune the phase relationships among them, in line with the cerebellum's
> established role in adaptation and coordination."

becomes:
> "...and tune the phase relationships among them, in line with the cerebellum's
> established role in locomotor adaptation, demonstrated in humans by impaired
> split-belt adaptation after cerebellar damage
> \citep{morton_cerebellar_2006}."

## 6. Qin, Li & Shen 2018 — visual-inertial drift correction (ocular feedback)

**Citation:** Qin, T., Li, P., & Shen, S. (2018). VINS-Mono: A robust and versatile
monocular visual-inertial state estimator. IEEE Transactions on Robotics, 34(4),
1004–1020. DOI: 10.1109/TRO.2018.2853729 (verified: IEEE, arXiv:1708.03852)

**Use:** Established robotics practice for exactly the mechanism Ben proposes:
vision correcting IMU drift via tight visual-inertial fusion. Grounds the ocular
feedback layer in prior art.

**Insert into:** Chapter 6, §6.3, this sentence:
> "...by treating visual measurements of the horizon and surrounding as a
> low-frequency correction to the vestibular estimate."

becomes:
> "...by treating visual measurements of the horizon and surroundings as a
> low-frequency correction to the vestibular estimate, mirroring the tight
> visual-inertial fusion that keeps drift bounded in state-of-the-art robot state
> estimators \citep{qin_vinsmono_2018}."

## 7. Song et al. 2021 — deep RL in neuromechanical simulation (state of the art)

**Citation:** Song, S., Kidziński, Ł., Peng, X. B., Ong, C., Hicks, J., Levine, S.,
Atkeson, C. G., & Delp, S. L. (2021). Deep reinforcement learning for modeling human
locomotion control in neuromechanical simulation. Journal of NeuroEngineering and
Rehabilitation, 18, 126. DOI: 10.1186/s12984-021-00919-y (verified: PMC8365920)

**Use:** State-of-the-art neuromechanical simulation built on OpenSim; the leading
alternative neural-control paradigm (deep RL) to our synthetic-nervous-system
approach. Citing it lets the dissertation position SNS control as biologically
interpretable where RL is a black box.

**Insert into:** Chapter 2, §2.8 Simulation Tooling, after the MyoConverter
sentence:
> "...and run the full neuromechanical loop faster than real time
> \citep{song_deep_2021, degroote_perspective_2021}."

## 8. De Groote & Falisse 2021 — perspective on predictive musculoskeletal simulation

**Citation:** De Groote, F., & Falisse, A. (2021). Perspective on musculoskeletal
modelling and predictive simulations of human locomotion. Proceedings of the Royal
Society B: Biological Sciences, 288(1946), 20202432. DOI: 10.1098/rspb.2020.2432
(verified: Royal Society listing)

**Use:** Recent state-of-the-art perspective on musculoskeletal simulation of human
locomotion — frames where the field is and what remains open. Good for the new
"State of the Art in Neuromechanical Simulation" section (drafted in
Notes/neuromechanics_state_of_art_and_ethics.tex).

**Insert into:** Chapter 2, §2.8 (with Song et al., above) and/or the new §2.9.
