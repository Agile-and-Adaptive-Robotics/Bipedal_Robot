import io

F = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\60-sns.tex"
t = io.open(F, encoding="utf-8").read()

start = t.index("The preceding chapters validate the robot's mechanics")
end = t.index("\\section{The AnimatLab Bipedal Walker}")
old = t[start:end]

new = r"""The preceding chapters establish how braided pneumatic actuator force and knee torque can be predicted. To study how neural control interacts with these mechanics, I also needed a simulation environment that connects anatomically organized spinal circuits to the robot's muscle routing and BPA force model. The existing tools reviewed in Chapter~\ref{ch:background} (Section~\ref{sec:bg_neuromech}) provide the neural, biomechanical, and physics components separately; using them for this robot required implemented interfaces and tests of the behavior preserved through conversion. The gap this chapter fills is that integration: no existing environment coupled an anatomically organized two-level controller, with Ia, Ib, group II, and contact feedback, to a BPA-actuated robot model with validated actuator forces.

In this chapter I develop and evaluate those interfaces in MuJoCo with SNS-Toolbox and in a graphical Simulink workflow, and I verify them rather than assume them: a controller ported between environments cannot be presumed equivalent, so each port is checked against its source implementation and, where the lineage's archived models allow, against the published behavior of the models it descends from. The complementary roles of the environments reflect the complementary needs of the program, fast scripted closed-loop simulation for tuning and evaluation, graphical circuit construction for assembly and inspection, and CAD-derived mechanics for the physical assemblies; Section~\ref{sec:sns_new_tools} states which capability each environment was missing before this work. I also examine whether selected behaviors persist when graded neural components are replaced with spiking components. The contribution is the integration and its verification, with the distinction from existing Simulink libraries established in Section~\ref{sec:sns_new_tools} by direct comparison. The neural vocabulary, the two-level central pattern generator with rhythm-generator and pattern-formation layers, the graded neuron and synapse models, and the afferent channels, is defined in Sections 2.5--2.8 and used here without redefinition.

The demonstrations establish specific component and closed-loop capabilities, and they identify failures of behavior transfer and limits of the assembled models with equal weight. Supported stepping and progression are distinguished from unsupported walking, which has not been demonstrated. Section~\ref{sec:prior_walker} records the AnimatLab walker whose hand-tuning motivated this verification-first procedure: manual tuning left the circuits' parameter sensitivity unmeasured, so the present procedure records performance under parameter and input perturbations instead. Section~\ref{sec:preliminary_sim_methods} is the methods section, Section~\ref{sec:preliminary_sim_results} the results, Section~\ref{sec:spiking_results} the spiking mirrors, and Section~\ref{sec:sns_discussion} the discussion.

"""

t = t[:start] + new + t[end:]
io.open(F, "w", encoding="utf-8").write(t)
print("SNS opening replaced:", len(old), "chars ->", len(new), "chars")
