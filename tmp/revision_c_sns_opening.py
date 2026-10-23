from pathlib import Path
R=Path(r'D:\Github\Bipedal_Robot');C=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C'
p=C/'chapters/60-sns.tex';s=p.read_text(encoding='utf-8-sig');end=s.index(r'\section{A Neuromechanical Simulation Toolchain}')
opening=r'''\chapter{Synthetic Nervous System Simulation}\label{ch:sns}
\graphicspath{{figs/60-sns/}{figs/60-sns/AnimatLab_reference/}}

The actuator and torque models of Chapters~\ref{ch:muscles} and~\ref{ch:torque} describe the robot's mechanics. Studying their interaction with neural control also requires a simulation that connects the mechanical model to a synthetic nervous system (SNS). I first attempted this integration in AnimatLab by combining a published walker with a spinal controller. The model produced stepping-like motion in air, but hand tuning did not produce walking. This chapter explains how that attempt motivated the development of software tools for assembling, testing, and comparing neuromechanical models.

\section{The AnimatLab Bipedal Walker}\label{sec:prior_walker}

Morrow and I combined the bipedal mechanics of \citet{li_bipedal_2017} with the two-level central pattern generator (CPG) of \citet{deng_neuromechanical_2019} in AnimatLab \citep{cofer_animatlab_2010}. The mechanical model had a pelvis and bilateral thigh, shank, foot, and toe segments, with sagittal-plane joints and antagonist Hill-type muscles at the hip, knee, and ankle. I added a passive toe joint. The controller separated rhythm generation from pattern formation and included muscle-length, velocity, and force feedback. Figure~\ref{fig:animatlab_model_architecture} shows the model and controller organization.

I adjusted applied currents, synaptic weights, and muscle parameters by hand. The suspended model developed alternating muscle activity and stepping-like motion, presented at the 2021 NIH BRAIN Initiative Investigators Meeting \citep{bolenTorqueValuesArtificialKnee2021June15-17}. It did not achieve weight-bearing walking. The muscle-feedback gains were not systematically tuned, and the combined model did not include the Li walker's ground-contact feedback.

The limitation was therefore more specific than an inability to generate a rhythm. The model generated alternating activity, but I had no systematic way to determine which parameter changes would preserve that activity as feedback and mechanical loading changed. A useful development environment needed to support repeatable tests of the controller, the actuator interface, and their coupling to the mechanical model. Further model schematics and recordings are provided in Appendix~\ref{app:simulation_verification}.

\section{Requirements for the New Tools}\label{sec:sns_gaps}

The existing workflow lacked four connections needed for the BPA-actuated robot: an SNS coupled to the human-scale and BPA-routed musculoskeletal models; transferable implementations of the archived walker controllers; graphical circuit construction linked to a maintained mechanical simulator; and interfaces joining CAD-derived bodies, BPA forces, and neural commands. These were integration requirements of this project, not capabilities absent from computational neuroscience as a whole.

The software reviewed in Chapter~\ref{ch:background} supplied the individual components. OpenSim provided musculoskeletal geometry, MyoConverter translated that geometry to MuJoCo, and SNS-Toolbox provided neuron and synapse models \citep{delp_opensim_2007, caggiano_myosuite_2022, todorov_mujoco_2012, nourse_sns_2023}. Simulink provided graphical circuit construction, and Simscape Multibody provided a route for importing mechanical geometry. My contribution was to connect these components for the models and actuators used here and to establish tests of the resulting implementations. Different environments served different tasks: scripted simulation and parameter search in MuJoCo, graphical model assembly in Simulink, and reference-model evaluation in AnimatLab.

'''
p.write_text(opening+s[end:],encoding='utf-8')
(C/'review/01_SNS_opening.md').write_text('Replaced the opening, prior-walker account, and gap framing. Removed universal novelty claims, parameter-tuning/cell-identity causation, repeated roadmaps, and undocumented claims about the extent of hand-tuning effort. Retained the Deng–Li combination, toe addition, air stepping, hand tuning, and failure to achieve weight-bearing walking. Supporting schematics/recordings are assigned to the supplementary appendix in the next section pass. No recorded test values were changed.\n')
print('SNS opening and requirements saved')
