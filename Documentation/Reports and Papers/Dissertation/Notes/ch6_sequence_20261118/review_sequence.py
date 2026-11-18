"""Apply the historical reviewer to the revised text and extend the single dated review table."""
import argparse
import html
import json
from pathlib import Path
import re
import sys
from zipfile import ZipFile

HERE = Path(__file__).resolve().parent
NOTES = HERE.parent
DIS = NOTES.parent
PREVIOUS = NOTES/'chapters_6_7_20261118'
sys.path.insert(0,str(NOTES/'advisor_feedback_history'))
import review


# Current line, prior lines, earlier reader question, source of this finding,
# observed problem, and the independent readback or verification.
CASES = [
    ('Opening',9,[7,11],['M-022','M-129'],'Ben','The old opening defined neuromechanics before readers knew why the AnimatLab workflow had to change.',
     'The first paragraph names the slow trial and the changes made by hand. Section 6.5 answers it with bounded searches and replay. No paired AnimatLab benchmark was recorded.'),
    ('Model sequence',44,[52,58,59,60,61,62],['M-022','M-065'],'Ben','Table 6.1 started with the reduced knee and arm before the walkers that motivated the work.',
     'Seven rows follow Li–Deng, bilateral W2L, the human reference, its rerouted robot counterpart, the CAD robot, the knee, and the arm.'),
    ('Source workflows',38,[42,46],['M-097','M-127'],'Ben','The model-preparation paragraph and old workflow did not show the two routes as the author presents them in the defense.',
     'Original slides 64 and 67 identify the supplied OpenSim and SOLIDWORKS models, Ben’s additions, and the two MuJoCo interfaces.'),
    ('Contact signal',31,[11,59],['M-023'],'Independent readers','An intermediate draft described foot signals from a treadmill experiment, which could be read as imported measurements.',
     'The bilateral circuit receives heel and toe feedback from the simulated feet during the supported trial.'),
    ('Prescribed activation',52,[60],['M-023','M-126'],'Independent dissertation reader','The model table could imply that a neural controller supplied the force-validation activations.',
     'The table and subsection distinguish prescribed OpenSim static-optimization activations from neural drive in movement trials.'),
    ('Human model',63,[88,90],['M-126','INTRO-T06'],'Ben and independent dissertation reader','The former generic description omitted the OpenSim source file, actuator repairs, and directly adjacent converter credit.',
     'The text identifies gait2392_simbody.osim, restored force scales, the 92 entries and six disabled rotators, and the MyoConverter citation.'),
    ('Force agreement',65,[154],['M-023'],'Independent dissertation reader','The reported lower end of the eight-muscle R² range was 0.87.',
     'The saved force comparison gives 0.8836081344–0.9769681745. Chapter 6 and Appendix D report 0.88–0.98; Appendix D’s per-muscle table was already correct.'),
    ('Robot routes',69,[88,94],['M-127','INTRO-T04'],'Ben and independent readers','The transition from a human muscle path to a robot route did not name the source, generated file, or the law used in the suspended result.',
     'The text follows gait2392_robotbody.osim through bilateral route preparation to gait2392_robot.osim. It identifies the Hill-type suspended example and the separate pressure-dependent BPA interface.'),
    ('CAD loading',79,[62,223],['M-023','M-105'],'Independent technical reader','The intermediate draft implied a causal test of whether loading itself cleared the swing foot.',
     'The supported CAD experiment checks stance loading and contralateral unloading and clearance together. Its 54 ideal tensions are distinguished from a pressure-derived force.'),
    ('Controller scale',83,[96,98],['M-097','M-129'],'Ben and independent readers','The GUI was described before the scale of the engineering problem and which circuit posed it.',
     'The saved full human W2Lvar topology has 888 nodes and 7,868 nonzero directed pairs. The text explains how to select the inputs to a knee motor pool and inspect their synaptic parameters.'),
    ('Stimulus and response',125,[169],['M-023','M-105'],'Ben','The earlier figure magnified only the stimulus. Its horizontal position made the stimulus look later than the response.',
     'Five source signals were checked in both views. Within each column every trace has the same zero and time scale. The 30 ms pulse precedes the 111 ms knee peak.'),
    ('Retest after CAD',141,[46,92],['M-023','M-127'],'Ben and both independent readers','“Repeat the same controller input” could imply clamped sensory feedback or a completed matched CAD experiment.',
     'The proposed test holds circuit, external drive, and load fixed, recomputes feedback from each body, and compares force and angle. Stretch applies in models with length feedback.'),
    ('Retest and optimization diagrams',154,[117,142],['M-097','M-105'],'Ben','The graphical controller and parameter search needed the workflows used in the defense.',
     'Slide 70 shows the saved W2L connections and CAD retest. Slide 71 supplies the three optimization operations. Both citation contexts point to Figure 6.6.'),
    ('Search sequence',162,[117,142],['M-129','M-105'],'Ben','The earlier headings separated the design problem from the procedure that answers it.',
     'The reader can select a function, identify a target, bound a search, inspect the measurements, and replay a fitted candidate. The chapter preserves the unsuccessful calibration cases.'),
    ('Conclusion',298,[282],['INTRO-T33','M-105'],'Ben','The earlier assistant draft called the reader “A colleague.” Ben had already written “A future researcher” in the supplied export.',
     'The conclusion uses Ben’s term, identifies repaired human muscle scales and the BPA force calculation, scopes the wall-time pair, and explains the test a future researcher can repeat.'),
]
AUTHOR_MARKS={'Opening':-2,'Model sequence':-2,'Stimulus and response':-3,'Conclusion':-1}


def write_json(name,value):
    (HERE/name).write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--export',type=Path,default=DIS/'Bolen_Dissertation.zip')
    args=parser.parse_args()
    with ZipFile(args.export) as z:
        before=z.read('chapters/60-sns.tex').decode().splitlines()
        before_d=z.read('chapters/95-AppendixD.tex').decode().splitlines()
    old=review.read_json(PREVIOUS/'PASSAGE_REVIEW.json')
    paths=review.DEFAULT_CHAPTERS+[review.CHAPTERS/'95-AppendixD.tex']
    current=list(review.passages(paths[0]))
    current_d=list(review.passages(paths[-1]))
    active={p.name:{review.digest(s) for _,_,_,s in review.passages(p)} for p in paths}
    kept=[p for p in old['passages'] if p['after_sha256'] in active[p['file']]]
    superseded=[p['id'] for p in old['passages'] if p not in kept]
    rows=[]
    for title,line,old_lines,examples,author,problem,answer in CASES:
        source=next(s for start,_,_,s in current if start==line)
        rows.append({'id':'SEQ-'+str(len(rows)+1).zfill(2),'title':title,'file':'60-sns.tex','line':line,
                     'source_of_finding':author,'original':'\n'.join(before[n-1] for n in old_lines),
                     'baseline_export_lines':old_lines,'finding':problem,'proposed_and_applied':source,
                     'after_sha256':review.digest(source),'historical_examples':examples,'reviewer_reason':answer,
                     'local_readiness':3,'advisor_score':None,'author_prompted_editor_deduction':AUTHOR_MARKS.get(title)})
    updated_d=next(s for n,_,_,s in current_d if n==74)
    rows.append({'id':'SEQ-D1','title':'Appendix force agreement','file':'95-AppendixD.tex','line':74,
                 'source_of_finding':'Independent dissertation reader','original':before_d[73],
                 'finding':'The same rounded lower bound appeared in Appendix D.','proposed_and_applied':updated_d,
                 'after_sha256':review.digest(updated_d),'historical_examples':['M-023'],
                 'reviewer_reason':'Matched to the muscle-specific values in Table D2 and the saved comparison CSV.',
                 'local_readiness':3,'advisor_score':None,'author_prompted_editor_deduction':None})
    decisions={'model':old['model'],'prompt':old['prompt']+' For each figure and table, compare its order and signal timing with the narrative, native data, and printed page. Test a proposed CAD change independently of a completed experiment.',
               'advisor_event_ids_at_review':old['advisor_event_ids_at_review'],'original_export':review.read_json(HERE/'TRANSFER_VERIFICATION.json')['baseline_sha256'],
               'passages':kept+rows,'superseded_prior_case_ids':superseded,
               'editorial_deductions_against_prior_assistant_draft':AUTHOR_MARKS,
               'interpretation':'The four negative marks apply Ben’s severity rubric to his direct feedback on the prior assistant draft. His own export had already corrected the Conclusion wording. The verified slide preference is a presentation decision. No new advisor event was observed.'}
    write_json('PASSAGE_REVIEW.json',decisions)
    screen=review.make_review(paths,review.signals(),decisions)
    assert not screen['saved_judgments_stale_or_absent'] and not screen['reopened_judgments']
    write_json('CHAPTER_SCREEN.json',screen)
    table=(PREVIOUS/'REVIEW_TABLE.html').read_text(encoding='utf-8')
    rendered=[]
    for row in rows:
        mark=row['author_prompted_editor_deduction']
        fields=(f"{row['source_of_finding']} · {row['title']}",f"{row['file']}, line {row['line']}",row['original'],
                row['proposed_and_applied'],row['finding'],row['reviewer_reason']+' Locally verified, 3/4.'+
                (f' Prior assistant draft, {mark} from Ben’s feedback.' if mark else ''))
        rendered.append('<tr data-group="Chapter 6 sequence">'+''.join('<td>'+html.escape(x)+'</td>' for x in fields)+'</tr>')
    headline=('<h1>Dissertation feedback and revisions</h1><p>Chapter 6 now follows the route from hand tuning to repeatable parameter trials. The first rows contain the author’s new directions and the two independent readers’ findings. The four new deductions total −8 against the previous assistant draft. The 429 recovered advisor signals retain their −551 total, including −42 from November.</p>'
              '<p><a href="Chapter_6_reading_pages.pdf">Bookmarked reading pages</a> · <a href="REVISION_REVIEW.md">Changes and evidence</a> · <a href="PASSAGE_REVIEW.json">Passage scores</a> · <a href="../chapters_6_7_20261118/SCORE_SUMMARY.json">Historical scores</a> · <a href="../advisor_feedback_history/README.md">Scoring rules</a></p>')
    table=re.sub(r'<h1>.*?</p>\s*<p>.*?</p>(?=\s*<nav>)',lambda _:headline,table,count=1,flags=re.S)
    table=table.replace('<option>Chapter reread','<option>Chapter 6 sequence</option><option>Chapter reread',1)
    count=table.count('<tr data-group=')
    table=table.replace('<tbody>','<tbody>\n'+'\n'.join(rendered),1)
    assert table.count('<tr data-group=')==count+len(rendered) and headline in table
    (HERE/'REVIEW_TABLE.html').write_text(table,encoding='utf-8')
    print('Table rows',count+len(rendered),'locally verified cases',len(decisions['passages']),
          'superseded',superseded,'screened',len(screen['passages']),'observed advisor score',screen['observed_advisor_score'])


if __name__=='__main__':
    main()
