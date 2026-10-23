from pathlib import Path
import re,json
C=Path(r'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022\Revision_C')
archive=C/'review/clarification_before';archive.mkdir(exist_ok=True)
for name in ['60-sns.tex','30-torque.tex','94-AppendixC.tex']:
    p=C/'chapters'/name;(archive/name).write_bytes(p.read_bytes())
p=C/'chapters/60-sns.tex';s=p.read_text(encoding='utf-8')
old='The limitation was therefore more specific than an inability to generate a rhythm. The model generated alternating activity, but I had no systematic way to determine which parameter changes would preserve that activity as feedback and mechanical loading changed. A useful development environment needed to support repeatable tests of the controller, the actuator interface, and their coupling to the mechanical model. Further model schematics and recordings are provided in Appendix~\\ref{app:simulation_verification}.'
new=r'''I subsequently expanded the neural model into left- and right-side circuits with interconnected rhythm generators, extending the controller to coordinate both legs. The right-side subsystems were constructed from the left-side organization, with commissural pathways linking the two rhythm generators. This bilateral construction was an extension of the combined model, rather than simply a reproduction of the published components. Further circuit schematics and implementation records are provided in Appendix~\ref{app:simulation_verification}.

The earlier air-stepping result did not establish whether the bilateral coupling and sensory pathways would coordinate the legs under load. When a configuration failed, I needed to distinguish a limitation of the controller design from an implementation error or a problem in its coupling to the mechanical model. Manual tuning alone did not resolve those questions. Addressing them required repeatable tests of the neural circuits, actuator interface, and closed-loop model, motivating the tools described below.'''
assert old in s;s=s.replace(old,new);p.write_text(s,encoding='utf-8')
p=C/'chapters/30-torque.tex';s=p.read_text()
start=s.index('Nine isometric holds have been recorded');end=s.find('\n\n',start)
if end<0:end=len(s)
s=s[:start]+r'''The extensor and flexor isometric measurements were reduced using the test-specific analysis scripts, which transform the measured load-cell wrench to the moving knee center. The reduced torques are compared with the corrected model and the human reference envelopes. The redesigned configurations meet or exceed the reference torques over most of the tested range of motion.'''+s[end:]
p.write_text(s,encoding='utf-8')
p=C/'chapters/94-AppendixC.tex';s=p.read_text();start=s.index(r'\section{Initial Cam-Route Isometric Holds}');removed=s[start:]
assert r'\label{tab:cam_holds}' in removed
(archive/'removed_excel_check_table.tex').write_text(removed,encoding='utf-8')
p.write_text(s[:start].rstrip()+'\n',encoding='utf-8')
p=C/'review/04_Torque.md';s=p.read_text();s=s.replace('New cam table transcribes only workbook values from ExtTest20mm_2; cached row 16 torque was not substituted for row 17 actual torque.','The spreadsheet-derived cam table was subsequently removed: Ben clarified that both Excel torque rows are sanity checks, and the test-specific reduction scripts supply the final torque data.');p.write_text(s,encoding='utf-8')
p=C/'review/06_Completed_study_wording.md';s=p.read_text();s=s.replace('The measurement table still contains only available workbook values.','The spreadsheet-derived measurement table has been removed; Excel rows 16 and 17 are sanity checks, not presentation torque data.').replace('and update the table/plots with the final data.','and use the final analysis-script outputs for reported values and plots.');p.write_text(s,encoding='utf-8')
(C/'review/07_Bilateral_and_torque_source.md').write_text('''# Clarifications from Ben

The main SNS opening now explicitly credits Ben's expansion of the combined walker into left- and right-side neural circuits with interconnected RGs. This follows the earlier Li/Deng construction and 2021 air-stepping account, and leads into the unresolved question of bilateral coordination under load and the need to separate circuit, implementation, and mechanical-interface effects. The later bilateral network is not presented as the source of the earlier air-stepping result. Appendix D retains its detailed assembly/test status.

Excel means the measurement spreadsheet. Both row 16 and row 17 are sanity checks, not the final presentation torque data. Removed the added spreadsheet-derived cam table, its six-torque range, and the main-text reference. Kept a copy of the removed material in clarification_before solely to document the correction. The existing torque_sheet_evidence.json is a superseded spreadsheet inspection, not a results source. Final reported torques and plots must come from Ben's test-specific reduction scripts and their outputs. No script was run and no replacement torque values were calculated in this correction.

The completed-study wording remains Ben's drafting instruction while final flexor measurements are collected tonight. New bipulley optimization on easteregg2 is in progress; the current dissertation figures describe the earlier recorded runs, not that active campaign. No new optimization outcomes are asserted.
''',encoding='utf-8')
p=C/'review/torque_sheet_evidence.json';j=json.loads(p.read_text());j['status']='SUPERSEDED: Excel torque rows 16/17 are sanity checks, not presentation data. Use test-specific analysis-script outputs.';p.write_text(json.dumps(j,indent=2),encoding='utf-8')
print('Clarified bilateral contribution and removed spreadsheet-derived result table')
