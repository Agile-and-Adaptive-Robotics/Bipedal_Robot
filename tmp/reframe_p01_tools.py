from pathlib import Path
import json,shutil
R=Path(r'D:\Github\Bipedal_Robot');E=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_E';V=E/'review/targeted_reconciliation';W=V/'tools_inventory';W.mkdir(exist_ok=False)
old=V/'P01_replacement.tex';shutil.copy2(old,W/'P01_previous_REJECTED.tex')
text=r'''The missing capability in this project was a reusable way to assemble neural controllers with the human musculoskeletal reference and artificial-muscle models, then evaluate their behavior under specified support, stimulation, and perturbation conditions. I addressed it by developing the following models and tools from the existing simulation packages:
\begin{enumerate}
    \item A custom Li--Deng walker and controller, developed in AnimatLab and implemented in MuJoCo with SNS-Toolbox, with matched implementations of its contact-driven neural core in NumPy and Simulink.
    \item The SNS-Simscape block library and GUI for block programming neural pathways, including sensory, neuron, synapse, activation, and custom BPA-force blocks, together with coupled reflex and CPG demonstrations.
    \item Gait2392-based neuromechanical models with alternative neural controllers, and a robot counterpart with BPA routing and the actuator-force model from Chapter~\ref{ch:muscles}.
    \item Integration of SNS-Simscape controllers with MuJoCo through the existing MathWorks bidirectional bridge, and model-import tools for bringing musculoskeletal and CAD-derived mechanisms into the MATLAB environment.
    \item An automated curriculum and experiment workflow that combines staged parameter searches with adjustable support, neural stimuli, pathway interventions, mechanical perturbations, and recorded neural and mechanical outputs.
    \item Implementations of literature-derived circuits and selected spiking-neuron counterparts for comparing neural mechanisms and coding assumptions.
\end{enumerate}
OpenSim, MuJoCo, SNS-Toolbox, Optuna, Simulink, and the MathWorks bridge provide the underlying software; the contributions here are the model implementations, custom blocks, integrations, and evaluation procedures. Table~\ref{tab:sns_capability} summarizes representative functional-verification tests and proof-of-concept demonstrations of these tools.
'''
old.write_text(text,encoding='utf-8')
table=r'''\begin{table}[htbp]
\caption[Tools developed and representative demonstrations]{Tools developed in this project and representative functional-verification tests and proof-of-concept demonstrations. Each result refers to the named circuit or mechanical model.}
\label{tab:sns_capability}
\small
\begin{tabularx}{\textwidth}{@{}>{\raggedright\arraybackslash}p{0.23\textwidth}*{2}{>{\raggedright\arraybackslash}X}@{}}
\toprule
\textbf{Tool developed} & \textbf{Function} & \textbf{Representative demonstration} \\
\midrule
Li--Deng/W2L model implementations & Combine the walker mechanics and neural architecture; transfer selected controller components between environments & Matched contact-driven rhythm cores reproduce \qty{1.000}{\second} periods with trace correlation at least 0.99996 \\
SNS-Simscape library and GUI & Construct neural pathways graphically and connect sensory, neural, and actuator blocks & The loaded-arm reflex demo gives \ang{2.00} sag with feedback and \ang{7.85} without it; actuator blocks agree with reference force calculations \\
Gait2392/SNS model implementations & Evaluate alternative controllers on the human reference and custom actuator/routing models & Individual-actuator tests check torque directions; selected human-model controllers produce bilateral stepping-like activity under support \\
W2L BPA-force integration & Apply the Festo force model within the coupled simulation & The 12-actuator BPA model demonstrates air stepping and retained-support standing, with peak pressure \qty{372}{\kPa} \\
SNS-Simscape/MuJoCo integration & Exchange activation commands and mechanical measurements through the MathWorks bridge & Single-muscle and control-rate tests reproduce the Python-reference sensor readback over the scored windows \\
Curriculum and experiment drivers & Automate staged tuning, support changes, neural stimuli, and mechanical perturbations & Saved stage-specific trials, candidate replays, and fixed-controller pathway and body-change comparisons demonstrate the workflow \\
Literature-circuit and spiking implementations & Test selected published neural behaviors and compare graded and spiking formulations & The reduced Rybak circuit passes 13 drive/deletion checks; selected spiking counterparts have paired implementation tests \\
\bottomrule
\end{tabularx}
\end{table}
'''
(V/'P01_companion_table.tex').write_text(table,encoding='utf-8')
inventory='''# Tools created and representative demonstrations — P01

## Editorial decision

Use **functional verification and proof-of-concept demonstrations**. “Functional verification” means checking the implemented function; “proof-of-concept demonstration” means using it in a representative model or experiment. Avoid “fundamental proof,” which implies a stronger mathematical or scientific result. P01 is a contribution inventory, not a list of unfinished work. Ben rejected the former incomplete-integration framing; its text is archived here, not recommended.

The paragraph and capability table are a paired targeted change. Replace the live paragraph beginning “I developed several software tools” through its capability-table sentence with `P01_replacement.tex`. Replace only the table environment bearing `tab:sns_capability` with `P01_companion_table.tex`. Preserve that label and existing review threads. If Ben’s current paragraph differs, reconcile its wording before applying. Do not upload an entire chapter. No Overleaf edits were made.

## Evidence-supported inventory

| Contribution | What can be claimed | Small representative demonstration and record |
|---|---|---|
| Custom Li–Deng/W2L model and controller implementations | Connor–Ben’s original shared-RG AnimatLab combination, Ben’s later per-side RG/feedback extension, and MuJoCo/SNS-Toolbox implementation; matched neural components in Simulink | Contact-driven core: 1.000 s periods, r≥0.99996, matching contact-removal latch. `spinal/campaigns/20260930/w2l_trace_compare.json`, `sns_simscape.md`. The Simulink result names the neural core, not a complete mechanical twin. |
| SNS-Simscape library and custom GUI | Editable sensory, neural, synaptic, activation, and biological/BPA force blocks; palette, parameter editing, and network import/export | `gui/SNSGui.m`, README GUI smoke record; actuator-reference checks; beer reflex ON/OFF 2.00°/7.85°. These are direct tool-function demonstrations. |
| Gait2392/SNS controller/model implementations | Repaired OpenSim-to-MuJoCo human reference, alternative s3k, w2lvar, and syn6 neural organizations, and custom routing/actuation models | Single-actuator sign/path tests; selected supported human stepping and recorded reference comparisons. Build census and saved parameters distinguish each neural variant. Do not imply identical outcomes for all controllers. |
| BPA-force integration | Custom Festo force computation coupled to mechanical routes; this is not functionality created by MuJoCo or SNS-Toolbox | `spinal/campaigns/20260930/bpa_run/bpa_gate.log`: all 12 W2L muscles use BPA force, bilateral air excursions and supported standing, peak 372 kPa within 620 kPa ceiling. This is distinct from the 92-entry Gait2392 human model and robotbody route. |
| Bidirectional SNS-Simscape/MuJoCo integration | Pairing the new neural library with the existing MathWorks bridge | `mujoco_bridge/BRIDGE_REPORT.md`, tests a/b/c and sensor readback; matched integration data. Credit the bridge authors. Bridge coupling does not convert a Simscape mechanical model into MuJoCo. |
| CAD/musculoskeletal import helpers | Imported mechanisms and model-building/import workflow in MATLAB | Saved SolidWorks XML/URDF imports, inventory logs, and OpenSim-derived Simscape model. Claim a complete imported robot with controller and BPA force only when its specific runnable saved model/result is identified. This is an evidence lookup, not a demand for new experiments. |
| Automated curriculum, optimization, and experimental controls | Five-stage variant curriculum, bounded Optuna searches, support scheduling, neural stimuli, pathway interventions, timed pushes, body-parameter changes, and output recording | `_curriculum_w2lvar.py`, `_curriculum_syn6.py`, saved stage trials and fresh replays; corrected W2L fixed-controller interventions. Optuna is reused; the staged objectives, drivers, coupling, and protocols are the project’s work. |
| Literature-derived circuit implementations | Reduced Rybak and Shevtsova neural models, Deng-derived circuits, Li-derived adaptation, and Di Russo-derived controller coupled to Gait2392 | Rybak `verify_rybak2006.log`: 13 checks including resetting/non-resetting deletions. Shevtsova has a reduced SNS implementation and recorded neural checks. Di Russo has stage-1 fit and completed stage-2 controller evaluations. List the actual demonstrated function, not successful reproduction of every paper’s locomotion. |
| Selected spiking neural counterparts | Spiking library blocks, selected demos, production-network mirror, and archived AnimatLab conversions | `sns_build_spiking.m`, `sns_units_test_spiking.m`, `build_network_spiking.py`, Appendix B’s implementation/calibration and paired comparisons. “Spiking mirrors of everything above” is broader than the identified records. |

## The five-stage sequence

The variant drivers explicitly implement (1) deafferented air stepping, (2) afferented air stepping, (3) standing balance, (4) walking-pattern tuning without ground contact, and (5) walking with ground contact. The stock `_curriculum.py` has a different historical fifth stage (PF-layer contact variant). Name the driver/sequence used; do not merge their stage labels. “Implemented curriculum” is the tool contribution and does not require claiming every candidate completed every stage.

## Wording of the gap

The project lacked a reusable way to assemble alternative neural controllers with a many-actuator human reference and BPA-actuated models, specify support and sensory conditions, and automate controlled evaluations. The contribution is the implemented models, block-programming environment, custom force integration, staged tuning, and experimental controls. Do not claim Python is inherently difficult or no comparable platform exists. The demonstrated operations explain the gap and its resolution without a subjective language comparison.

Use “many-actuator musculoskeletal model” or “muscle-redundant model.” The latter denotes more muscle actuators than joint degrees of freedom; “over-actuated” may suggest independently controllable actuator directions and needs a more specific control definition.

## Claim attribution and precise placements

- Rybak’s verified implementation located here is Python, not a demonstrated Simulink Rybak model. Put it under the literature-model toolset until the claimed Simulink filename is supplied.
- Shevtsova’s standalone/SNS lumbar circuit is not evidence of a completed Shevtsova–Gait2392 coupling. Shinohara rules/specification are not a runnable Shinohara–Gait2392 model. Di Russo–Gait2392 coupling exists and can be listed as an implemented model.
- A Simulink controller communicating with a MuJoCo model through a bridge does not establish a CAD→Simscape→MuJoCo mechanical export route. List CAD import helpers and neural/physics integration separately unless the exporter/result is identified.
- The new ablation report is a contact+force ablation with Ia/II retained, not full deafferentation. It also retains stale zero-sensitivity interpretations and calls a W2L model BPA without an actuator manifest; do not import those claims into P01 without matching the actual model. Its saved metrics can be evaluated separately.

These are claim-boundary and attribution notes for the editors. The proposed P01 prose and table lead with delivered tools and their demonstrations, not outstanding tasks. No new computation was launched.
'''
(W/'TOOLS_AND_DEMONSTRATIONS.md').write_text(inventory,encoding='utf-8')
patches=json.loads((V/'patches.json').read_text(encoding='utf-8'))
for p in patches:
 if p['id']=='P01':p.update(action='Paired paragraph-and-table replacement: list tools created and relevant demonstrations',old='Whole current paragraph beginning I developed several software tools through the capability-table sentence, plus the table environment labeled tab:sns_capability.',new=text,status='Supersedes the rejected incomplete-integration wording. Apply the companion table with this paragraph; exact current author wording must be checked before application.')
(V/'patches.json').write_text(json.dumps(patches,indent=2),encoding='utf-8')
md=V/'OVERLEAF_TARGETED_PATCHES.md';prior=md.read_text(encoding='utf-8');shutil.copy2(md,W/'PATCH_SHEET_before_P01_update.md')
header='''# P01 revised by Ben — tools and demonstrations

The former P01 incomplete-integration wording is **rejected and superseded**. Use the updated `P01_replacement.tex` and paired `P01_companion_table.tex`. The paragraph lists created models/tools and the table lists representative functional-verification tests and proof-of-concept demonstrations. See `tools_inventory/TOOLS_AND_DEMONSTRATIONS.md` for evidence, attribution, and precise scope. No new work list is inserted into P01. The historical patch sheet below remains for the other patches; its old P01 quotation is not to be applied.

---

'''
md.write_text(header+prior,encoding='utf-8')
p=R/'CHATGPT_REPORT.md';prior=p.read_text(encoding='utf-8');p.write_text('# P01 reframed: tools created and representative demonstrations\n\nBen rejected the incomplete-work framing. The current P01 fragments list delivered model/tool contributions and functional-verification/proof-of-concept demonstrations, with credit to the reused packages. Evidence and naming boundaries are in Revision_E/review/targeted_reconciliation/tools_inventory/TOOLS_AND_DEMONSTRATIONS.md. The paragraph and capability table are paired targeted replacements, not chapter uploads. No Overleaf edits or new model runs were made.\n\n---\n\n'+prior,encoding='utf-8')
print('Revised P01 paragraph + companion table; evidence inventory saved. No manuscript or Overleaf source changed.')
