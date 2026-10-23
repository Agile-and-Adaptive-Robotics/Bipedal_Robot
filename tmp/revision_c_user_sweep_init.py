from pathlib import Path
import shutil,json,hashlib
R=Path(r'D:\Github\Bipedal_Robot');C=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C';S=C/'review/user_sweep_20261023';S.mkdir(exist_ok=True)
assert not (S/'before').exists(),'Already initialized'
shutil.copytree(C/'chapters',S/'before/chapters');shutil.copytree(C/'figs',S/'before/figs')
items=[
('01','Remove the quoted stiffness-pair nonuniqueness/external-overshoot-selection paragraph; do not substitute a defensive Pareto justification.'),
('02','Replace Fig. 2.3 CAD BPA with supplied inflated/deflated photo; preserve originals, give both camera photos descriptive filenames.'),
('03','Fig. 2.4: uppercase panel letters.'),('04','Fig. 2.5: panel letters at upper left.'),('05','Reorganize Figs. 2.5–2.6: skeleton, robot PAM routes, human Gait2392 paths together; remove repeated models and awkward letters.'),
('06','Resize fig:Fmax_fit; fit Figs. 3.3 and 3.4 within page text area without altering MATLAB plot data.'),('07','Center Table 3.3 and increase text size; preserve every value.'),
('08','Add diagram connecting Eq. 4.2 dot-product moment-arm construction to cross-product torque method; reuse appropriate Fig. 2.6 geometry where useful.'),('09','Fig. 4.4 panel letters bold.'),('10','Locate and show bending schematic with Eq. 4.24.'),
('11','Shorten captions of Tables 4.1–4.3; move essential scoring/fold context to prose, preserve values.'),('12','Table 4.4 chi0 in meters, with correct unit conversion.'),('13','Explain anatomical origin/insertion search bounds, biomimetic routing/RoM tradeoff, and subsequent force/part-failure tradeoff.'),
('14','Fig. 4.11 panel C showing torque magnitude across RoM: ZCode-owned; do not duplicate concurrent plot work.'),('15','Add supplied cam-knee hardware photograph to Fig. 4.12.'),
('16','SADb gap: structured database rather than structured account.'),('17','SADb source libraries: Zotero libraries of AARL and Ben, plus papers already in Airtable.'),('18','Read SfN24_Bolen.pdf; replace vague AI claims with evidence-based account of limited help and categorization difficulties.'),
('19','SNS 6.1: initial Connor/Ben controller had ONE RG controlling both sides; air stepping at NIH2021; Ben added afferent and heel/toe feedback; distinguish later bilateral RG expansion.'),('20','State AnimatLab unsupported, based on old physics engine, and slow in this workflow.'),
('21','ZCode verification: Table 6.1 split verification/limits; identify actual authored models versus existing MuJoCo/SNS tools, external MuJoCo–Simulink bridge attribution and direction, controller ports/demos/curriculum scope, resolved welded joint, actual muscle integration status. Do not guess answers.'),
('22','Remove unsupported Hunt/Young neuromechanical Simulink-library comparison; flag prior advisor reply for correction.'),('23','Define W2L at first use, grounded in actual model naming.'),('24','Gap includes tested neuromechanical models with Gait2392-scale actuator counts and degrees of freedom; scope novelty honestly.'),('25','Fresh full build, render changed pages, conservation/line-evidence audit, clean upload bundle, sync authorized finished changes to Overleaf. No report uploads.')]
data={k:{'status':'ZCODE' if k in ['14','21'] else 'TODO','request':v,'evidence':''} for k,v in items}
(S/'checklist.json').write_text(json.dumps(data,indent=2),encoding='utf-8')
block='\n\n## Ben\'s figure/prose correction sweep — 2026-10-23\n\nWorking copy: `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C`. Originals in ProofFinal remain unchanged. Ben explicitly authorizes Chapter 3 layout changes in this pass. ZCode-owned items are requests to verify, not facts to invent. Progress is updated in this checklist and `review/user_sweep_20261023/checklist.json`.\n\n'+'\n'.join(f'- [ ] U{k} [{data[k]["status"]}] {v}' for k,v in items)+'\n'
with (R/'CHATGPT_REPORT.md').open('a',encoding='utf-8') as f:f.write(block)
print('Saved 25-item handoff checklist and complete before snapshot')
