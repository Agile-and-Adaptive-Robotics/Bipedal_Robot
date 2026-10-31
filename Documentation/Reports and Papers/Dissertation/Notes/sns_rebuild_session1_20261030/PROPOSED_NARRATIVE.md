# Proposed SNS narrative passages — Session 1

Review text only. Cross-references, citation keys, and the exact Hunt slide wording must be reconciled before insertion. This proposal has not changed the active manuscript.

## Replacement for 6.2: AnimatLab Results and Lessons

The suspended AnimatLab model provided an initial test of the connection between neural activity and joint motion. In a saved graded-neuron trial, alternating extensor and flexor activity propagated from the rhythm generator through the hip pattern-formation and motor-neuron subsystems. Hip, knee, and ankle angles varied periodically over the same interval. The oscillation frequency was approximately 2.2 Hz. A short recording of these signals shows the relationship between the neural layers and the resulting movement more clearly than a long record of rhythm-generator activity alone.

Suspended motion did not establish that the controller could support the model on the ground. In a separate unsupported ground test, the pelvis fell while the neural rhythm continued. The two tests therefore exposed different requirements: producing coordinated joint motion and regulating muscle activity during loading. The ground test also showed why an oscillating neural signal could not serve as the sole criterion for successful locomotor control.

These experiments motivated a workflow that could preserve a controller's structure while testing it against different mechanical models and sensory conditions. Graphical editing remained useful for inspecting connections and applying stimuli, but repeated controller fitting and mechanical tests required an efficient programmatic interface. The following sections describe the tools developed to combine these capabilities and the demonstrations used to check them.

*Evidence note, not manuscript prose:* The suspended and ground examples are separate archived configurations and trials. Do not write “the model was then dropped” or imply one continuous matched experiment. The proposed figure uses the original AnimatLab angle conventions; these differ from the converted Gait2392 knee convention. Add the precise source/configuration information to the appendix rather than the main narrative.

## Replacement opening for 6.3: Knowledge Gaps and New Tools

The next problem was transferring a neural controller between mechanical models. A controller that coordinates a simplified walker must be adapted when muscle geometry, joint motion, inertia, or sensory signals change. This issue becomes more consequential when moving to a detailed human musculoskeletal model such as Gait2392 and then to a robot actuated by braided pneumatic actuators. The controller must be evaluated with the mechanics and actuator model that it will control.

I therefore needed a workflow that supported both graphical circuit development and repeated numerical experiments. Graphical blocks make it possible to inspect a neural subsystem, change its connections, and connect it to a mechanical model. Headless execution runs an experiment without opening its graphical interface, allowing many tests to be automated. An application programming interface (API) exposes the commands needed to construct, control, and record those experiments; a software development kit (SDK) provides supporting libraries and tools for using such interfaces.

The tools in Table 6.1 address these requirements in stages. The Li–Deng and W2L models provide initial controller tests, and the Gait2392-based models add detailed musculoskeletal mechanics. The SNS–Simulink library supports graphical neural development, while Simscape connects those controllers to mechanical components and CAD-derived assemblies. The MuJoCo–Simulink interface allows tested control and sensor signals to pass between the simulation environments. Experiment drivers then support repeated training, perturbation, and curriculum tests, with literature-derived circuits and spiking variants providing additional controller options.

*Scope note, not manuscript prose:* Avoid claiming that these tools automatically translate every network or mechanical model between programs. Explain the demonstrated interfaces and their limitations where their results are presented. Avoid the blanket assertion that Python has no graphical neural-development tools.

## Proposed biarticular transition

Biarticular muscles introduce an additional question when organizing the pattern-formation layer. Their activation affects more than one joint, so grouping them solely by one anatomical joint can conceal their mechanical role. Deng et al. (2023) examined the contribution of biarticular muscles to the stability of a rat-hindlimb neuromechanical model. This work motivates evaluating biarticular coordination explicitly when adapting the controller to another musculoskeletal model. The present implementations and their connection choices must still be evaluated on their own mechanical models.

*Scope note:* Describe the implemented four-cell RF/hamstring subsystem as an experimental architecture when that is the configuration shown. Keep sartorius and gracilis separate, as Ben directed. Do not imply that the four-cell design or its particular connection strengths were established by the cited study without a direct circuit comparison.

## Proposed tools-to-curriculum transition

The available tools allow individual components to be checked before the full model is asked to step. Hunt's training method provides the organizing framework for fitting neural and mechanical responses, while the five-task curriculum defines the physical conditions under which the controller is evaluated. These serve different purposes: training determines controller parameters, and the curriculum tests the resulting configuration under progressively more demanding conditions. Harness-supported ground stepping requires both support of the stance leg and measurable unloading and clearance of the swing leg. An oscillating suspended model does not meet those criteria.

*Citation note:* Cite Hunt's February 15, 2019 INRC2 presentation, *Real-time control of dynamic walking in a dog-like robot*. Final wording of the four operations requires checking that presentation. Do not retain the 2017 citation as the source of this lecture method.

## Proposed closing direction

End with a specific CAD-derived robot demonstration: identify the neural controller, BPA force model, joint/contact conditions, and external support, then state the observed behavior. Connect the demonstration to the next chapter's robot-design requirements. Keep this short and positive about the implemented coupling, while distinguishing supported standing from locomotion.

Do not finalize that paragraph around an unselected recording. The existing MuJoCo hopping GIFs use ideal reference tensions and an engineering allocator; the author-confirmed Simscape/BPA coupling is a separate contribution. Select and identify the matching standing demonstration before adding a numerical result or caption.
