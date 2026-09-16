# Overleaf review notes for zcode

This file records every unresolved review comment and tracked change authored by `ajh26` that was visible in Overleaf's **Review → Overview** panel. Dates are copied as Overleaf displayed them; specific times are omitted. Items authored by `bbolen83` are excluded.

Inventory: **66 individual `ajh26` items** — 49 comments and 17 tracked changes — across `02-abstract.tex`, `10-introduction.tex`, `15-background.tex`, `20-methods.tex`, and `main.tex`.

## 02-abstract.tex

### Comment

Date: 10 September  
Section: abstract  
Lines: 1–5  
Text highlighted:

```latex
\chapter*{\centerline{Abstract}}
\markboth{\MakeUppercase{Abstract}}{}
\iftoggle{fulltoc}{
  \addcontentsline{toc}{chapter}{Abstract}
}{}
```

ajh26 commented: "Pull your abstract down to one page. 'Extra' material can be moved to the introduction and 'dissertation organization' section."

### Tracked deletion

Date: 10 September  
Section: abstract  
Line: 13  
Text deleted:

```latex
Planned extensions add vestibular (IMU), ocular, and cerebellar layers, culminating in a treadmill-walking robot with fully embedded real-time neural control on Teensy and Jetson-class hardware.
```

## 10-introduction.tex

### Comment

Date: 10 September  
Section: 10-introduction  
Line: 7  
Text highlighted:

```latex
. 
```

This is the period and following space between the sentence ending

```latex
\citep{hunt_modeling_2017, sarosi_comparative_2017, liang_comparative_2020}.
```

and the sentence beginning

```latex
Designing a biomimetic robot around BPAs therefore requires quantitative, predictive models of the actuator and of the joint torques it can produce.
```

ajh26 commented: "Put in a clear sentence here describing the gap(s) in knowledge."

### Comment

Date: 10 September  
Section: 10-introduction  
Line: 15  
Text highlighted:

```latex
Two
```

ajh26 commented: "This does not match with the previous paragraphs which is all preliminary work. This is work you did during the PhD, as such it should either be moved to a chapter (included in Dessertation Organization) or the discussions section."

### Tracked replacement

Date: 10 September  
Section: 10-introduction  
Line: 21  
Changed:

```latex
cannot
```

to:

```latex
do not accurately
```

### Tracked deletion

Date: 10 September  
Section: 10-introduction  
Line: 21  
Text deleted:

```latex
maximum 
```

### Comment

Date: 10 September  
Section: 10-introduction  
Line: 21  
Text highlighted:

```latex
blocks rational actuator selection
```

ajh26 commented: "reword"

### Tracked replacement

Date: 10 September  
Section: 10-introduction  
Line: 23  
Changed:

```latex
s
```

to:

```latex
 force
```

This changes `BPAs` to `BPA force` in the sentence.

### Tracked deletion

Date: 10 September  
Section: 10-introduction  
Line: 23  
Text deleted:

```latex
, and the isometric test system itself is compliant
```

### Tracked addition

Date: 10 September  
Section: 10-introduction  
Line: 23  
Text added:

```latex
 best
```

### Tracked replacement

Date: 10 September  
Section: 10-introduction  
Line: 23  
Changed:

```latex
, and
```

to:

```latex
. Addtionally,
```

### Tracked addition

Date: 10 September  
Section: 10-introduction  
Line: 23  
Text added:

```latex
or compared 
```

### Comment

Date: 10 September  
Section: 10-introduction  
Line: 33  
Text highlighted:

```latex
Section
```

The highlighted word occurs in `(Section~\ref{sec:results_torque})`.

ajh26 commented: "This is a whole paper, seems weird it is just a section while the previous paper is a chapter and a half."

### Comment

Date: 10 September  
Section: 10-introduction  
Line: 39  
Text highlighted:

```latex
Chapter
```

This is the first word of the `Dissertation Organization` paragraph.

ajh26 commented: "Write a whole paragraph for each chapter. Could probably be combined with the Research Objectives section."

## 15-background.tex

### Comment

Date: 10 September  
Section: 15-background  
Line: 1  
Text highlighted:

```latex
\chapter{Background}\label{ch:background}
```

ajh26 commented: "My opinion is to feel not guilty about borrowing images from other papers for a dissertation background section (Just make sure to cite it!). That way it won't take a ton of time and is helpful in getting your readers up to speed without sending them directly to the papers."

### Comment

Date: 10 September  
Section: 15-background  
Line: 5  
Text highlighted:

```latex
of freedom and joint torques of the animal limbs they mimic.
```

ajh26 commented: "At the end of each section/paragraph, you can add a sentence directly connecting it to your PhD work. Example: The work in this dissertation will help develop biomimetic robots that more accurately captures animal joint torque profiles."

### Comment

Date: 10 September  
Section: 15-background  
Line: 9  
Text highlighted:

```latex
Compliant
```

ajh26 commented: "Doesn't follow from topic sentence/previous sentences."

### Comment

Date: 10 September  
Section: 15-background  
Line: 9  
Text highlighted:

```latex
Recent
```

ajh26 commented: "Doesn't follow from other sentences"

### Comment

Date: 10 September  
Section: 15-background  
Line: 9  
Text highlighted:

```latex
reviews
```

ajh26 commented: "plural?"

### Comment

Date: 10 September  
Section: 15-background  
Line: 9  
Text highlighted:

```latex
}.
```

This is the end of `\citep{zhagiparova_recent_2025}.`.

ajh26 commented: "Again, at the end, connect to your work."

### Comment

Date: 10 September  
Section: 15-background  
Line: 24  
Text highlighted:

```latex
\begin{multline}
```

ajh26 commented: "formatting"

### Comment

Date: 10 September  
Section: 15-background  
Line: 30  
Text highlighted:

```latex
and the relative contraction $\epsilon^*$ used in that model are defined in Section~\ref{sec:force_model}.
```

ajh26 commented: "there is no e* above."

### Comment

Date: 10 September  
Section: 15-background  
Line: 32  
Text highlighted:

```latex
Chapter~\ref{ch:results} compares these models with data collected for this dissertation and develops an improved characterization.
```

ajh26 commented: "Yes, like this!"

### Comment

Date: 10 September  
Section: 15-background  
Line: 36  
Text highlighted:

```latex
OpenSim
```

ajh26 commented: "What is opensim?"

### Comment

Date: 10 September  
Section: 15-background  
Line: 38  
Text highlighted:

```latex
Steele's biomimetic knee joint, a 1-DoF four-bar linkage whose crossed links produce a migrating instantaneous center of rotation similar to the human knee, provides the physical knee used in this dissertation \citep{steele_development_2017, steele_biomimetic_2018, steele_biomimetic_2018-1}.
```

ajh26 commented: "Turn this into a whole paragraph with figure"

### Tracked deletion

Date: 10 September  
Section: 15-background  
Line: 42  
Text deleted:

```latex
spinal 
```

### Comment

Date: 10 September  
Section: 15-background  
Line: 50  
Text highlighted:

```latex
steering
```

ajh26 commented: "wasn't this balance?"

### Comment

Date: 10 September  
Section: 15-background  
Line: 50  
Text highlighted:

```latex
inter-leg coordination analyses in stick-insect-inspired controllers \citep{nourse_analyzing_2019}, reinforcement-reflex hybrids in insects \citep{goldsmith_investigating_2021}
```

ajh26 commented: "Not AARL, but Quinn lab"

### Comment

Date: 10 September  
Section: 15-background  
Line: 50  
Text highlighted:

```latex
adaptive hindlimb walking controllers \citep{schilling_adaptive_2022}
```

ajh26 commented: "not the AARL"

### Comment

Date: 10 September  
Section: 15-background  
Line: 50  
Text highlighted:

```latex
.
```

This is the period after `\citep{scharzenberger_design_2019}`.

ajh26 commented: "How about the rat stuff from myself, Deng, Young, and Jackson?"

### Tracked replacement

Date: 10 September  
Section: 15-background  
Line: 52  
Changed:

```latex
planned
```

to:

```latex
discussed
```

### Comment

Date: 10 September  
Section: 15-background  
Line: 58  
Text highlighted:

```latex
\section{State of the Art in Neuromechanical Simulation}\label{sec:bg_sota}
```

ajh26 commented: "Most of this information is elsewhere in your dissertation, and I don't really like this section. What isn't elsewhere can probably moved into the areas that talk about related items."

### Comment

Date: 10 September  
Section: 15-background  
Line: 60  
Text highlighted:

```latex
plant
```

ajh26 commented: "model?"

### Comment

Date: 10 September  
Section: 15-background  
Line: 60  
Text highlighted:

```latex
has advanced
```

ajh26 commented: "it is not clear who has done this advancing"

### Comment

Date: 10 September  
Section: 15-background  
Line: 70  
Text highlighted:

```latex
and the neuromechanical controller can be verified against published tonic stimulation and deletion phenomena \citep{selionov_tonic_2009, rybak_modelling_2006} before any locomotion tuning is attempted.
```

ajh26 commented: "I think you should pull back a little on how you refer to the 'future work'. 2 reasons. 1) YOu need to not confuse your committee on what is completed vs what is not. 2) You need to frame your dissertation as providing real benefits as is. The benefits include the cool experiments that will be performed in the future. The way this sentence is written kind of muddies the water."

### Tracked replacement

Date: 10 September  
Section: 15-background  
Line: 70  
Changed:

```latex
will
```

to:

```latex
can
```

### Tracked deletion

Date: 10 September  
Section: 15-background  
Line: 72  
Text deleted:

```latex
 at all
```

## 20-methods.tex

### Comment

Date: 10 September  
Section: 20-methods  
Line: 10  
Text highlighted:

```latex
wrapping
```

ajh26 commented: "Also used to design attachment locations for a robotic knee that meets the torque capabilities of a human knee."

### Comment

Date: 11 September  
Section: 20-methods  
Line: 173  
Text highlighted:

```latex
independent
```

ajh26 commented: "x and z are not independent. Make sure you discuss why. Perhaps set them up as independent, then discuss why they equal each other and make the new equation eliminated x3"

### Comment

Date: 11 September  
Section: 20-methods  
Line: 192  
Text highlighted:

```latex
Force vector in a body's frame can be written as $\mathbf{F}$.
```

ajh26 commented: "I'm not sure what this means"

### Comment

Date: 11 September  
Section: 20-methods  
Line: 192  
Text highlighted:

```latex
Force vector in a body's frame can be written as $\mathbf{F}$.
```

ajh26 commented: "OK, I think I get it now, but it should be made a little more clear where this value comes from."

### Comment

Date: 11 September  
Section: 20-methods  
Line: 254  
Text highlighted:

```latex
This
```

ajh26 commented: "With a quick read, I believe all this. Have Stu double check the equations."

### Comment

Date: 11 September  
Section: 20-methods  
Line: 274  
Text highlighted:

```latex
Newton-Raphson method to solve the nonlinear force-balance equation
```

ajh26 commented: "Stop conditions?"

### Comment

Date: 11 September  
Section: 20-methods  
Line: 280  
Text highlighted:

```latex
to update the BPA attachment location
```

ajh26 commented: "I'm not sure I understand the point of this. Is this for the optimization later? Or are you modeling the compliance change in the attachment? You then next talk about optimized stiffness terms, but it isn't clear yet that those values were optimized."

### Comment

Date: 11 September  
Section: 20-methods  
Line: 280  
Text highlighted:

```latex
In
```

ajh26 commented: "New section start here to talk about otpimization of parameters."

### Comment

Date: 11 September  
Section: 20-methods  
Line: 280  
Text highlighted:

```latex
In
```

ajh26 commented: "Make sure you do a good job justifying the use of optimization instead of analytical or theoretical methods."

### Comment

Date: 11 September  
Section: 20-methods  
Line: 292  
Text highlighted:

```latex
stiffness
```

ajh26 commented: "Some random thoughts I just had. Since muscles have similar force-length curves, do we not see a lot of large squishy creatures because of the loss in force? Also, not a lot of 'protruding' attachment points? Might be worth pointing out in the discussion after more rumination."

### Comment

Date: 11 September  
Section: 20-methods  
Line: 292  
Text highlighted:

```latex
Because the bracket forms an arch shape and extends for \qty{75}{\mm} in the $Z$ direction,
```

ajh26 commented: "This is confusing. Do you just mean because it wraps around and makes contact with the knee?"

### Comment

Date: 11 September  
Section: 20-methods  
Line: 292  
Text highlighted:

```latex
Because the bracket forms an arch shape and extends for \qty{75}{\mm} in the $Z$ direction,
```

ajh26 commented: "Also, are you starting a new subject here...new paragraph?"

### Tracked replacement

Date: 11 September  
Section: 20-methods  
Line: 313  
Changed:

```latex
 only
```

to:

```latex
,
```

### Comment

Date: 11 September  
Section: 20-methods  
Line: 317  
Text highlighted:

```latex
motion
```

ajh26 commented: "This whole section about optimization is quite confusing and without much context."

### Comment

Date: 11 September  
Section: 20-methods  
Line: 321  
Text highlighted:

```latex
MATLAB
```

ajh26 commented: "This is a much better way to introduce the optimization."

### Comment

Date: 11 September  
Section: 20-methods  
Line: 323  
Text highlighted:

```latex
described above
```

ajh26 commented: "?"

### Comment

Date: 11 September  
Section: 20-methods  
Line: 330  
Text highlighted:

```latex
prediction fails
```

ajh26 commented: "?"

### Comment

Date: 11 September  
Section: 20-methods  
Line: 334  
Text highlighted:

```latex
count
```

ajh26 commented: "Things in this section need to be explained more clearly and clearly justified."

### Comment

Date: 11 September  
Section: 20-methods  
Line: 336  
Text highlighted:

```latex
\section{Prior Neuromechanical Simulation: The AnimatLab Bipedal Walker}\label{sec:prior_walker}
```

ajh26 commented: "Should probably be in the background."

### Tracked deletion

Date: 11 September  
Section: 20-methods  
Line: 338  
Text deleted:

```latex
The course report describing that model \citep{bolen_design_2020} is summarized here because its biomechanics and network architecture are the direct ancestors of both the current AnimatLab walker model, whose neuron and synapse parameters are documented in Appendix~\ref{app:neuron_equations}, and the two-layer CPG controller planned in Chapter~\ref{ch:futurework}.
```

### Tracked deletion

Date: 11 September  
Section: 20-methods  
Line: 340  
Text deleted:

```latex
; its contribution to push-off would make an interesting study in its own right once the model walks on the ground
```

### Tracked addition

Date: 11 September  
Section: 20-methods  
Line: 340  
Text added:

```latex
actively
```

### Tracked deletion

Date: 11 September  
Section: 20-methods  
Line: 340  
Text deleted:

```latex
, so lateral stabilization falls to the vestibular layer and hardware design of Chapter~\ref{ch:futurework}
```

### Comment

Date: 11 September  
Section: 20-methods  
Line: 342  
Text highlighted:

```latex
limitations that directly motivated the systematic, verification-first tuning approach of Chapter~\ref{ch:futurework}
```

ajh26 commented: "These types of sentences can maybe be rephrased to stand alone better and be more like the structure of a background."

### Comment

Date: 11 September  
Section: 20-methods  
Line: 346  
Text highlighted:

```latex
Preliminary
```

ajh26 commented: "Rewrite this section to have more of a tone that you developed a tool or a pipline for use. Get rid of the word preliminary."

## main.tex

### Comment

Date: 21 July  
Section: main  
Lines: 145–146  
Text highlighted:

```latex
% Motivation (no chapters/08-motivation.tex yet --- uncomment when written)
% \include{chapters/08-motivation}
```

ajh26 commented: "It is fine to have a personal motivation section, but you need an introduction that motivates the work from a scientific standpoint. Take about 2-3 pages to build a good introduction motivation, and then create a paragraph summarizing each of your subsequent chapters."

### Comment

Date: 21 July  
Section: main  
Line: 148  
Text highlighted:

```latex
%-------------------------
```

This is the beginning of `%-------------------------Listings-----------------------------`.

ajh26 commented: "You also need a Background chapter. This is where you flex your knowledge for your committee."

### Comment

Date: 21 July  
Section: main  
Line: 168  
Text highlighted:

```latex
\include{chapters/20-methods}
```

ajh26 commented: "I think you should actually format this more like a big paper. Pull both your introductions into 1 chapter. Your Methods into 1 chapter. Your Results into 1 chapter. And your discussion into 1 chapter. A bit like the 'original' paper."

### Comment

Date: 21 July  
Section: main  
Line: 170  
Text highlighted:

```latex
\include{chapters/40-discussion}
```

ajh26 commented: "To be removed."

