"""Apply the historical reader to the current text and extend the cumulative review table."""
import html
import json
from pathlib import Path
import re
import sys
from zipfile import ZipFile


HERE = Path(__file__).resolve().parent
NOTES = HERE.parent
PRIOR = NOTES/'ch6_sequence_20261118'
sys.path.insert(0,str(NOTES/'advisor_feedback_history'))
import review


# Each row anchors its original passage to the preceding dated transfer and its
# answer to the new manuscript. A single author observation can affect many rows.
# Current line, earlier line, title, reader, missing context, verified response
CASES = [
    (9,9,'AnimatLab','Ben','The opening used AnimatLab as a familiar program name.',
     'The opening now identifies graphical software coupling a neural circuit to a moving body and describes the variables Morrow and Bolen adjusted.'),
    (11,11,'Li–Deng and CPG','Independent dissertation reader','The study names, RG, and PF identified a source without explaining which parts produced the rhythm or timed the muscles.',
     'Li supplied the sagittal-plane biped. Deng supplied the two-level circuit. RG sets its period and PF times motor recruitment. The chapter introduces both before Figure 6.1.'),
    (31,31,'MuJoCo and SNS-Toolbox','Ben','The reader had to infer whether MuJoCo and SNS-Toolbox were competing simulators or separate parts of the test.',
     'The physics engine calculates motion and contact. The Python library simulates the neurons and synapses coupled to it.'),
    (33,31,'Walker Two-Layer','Ben, earlier directive U23','The earlier chapter wrote “In a later walker model, W2L” and had lost an already documented author correction.',
     'The first named appearance spells out Walker Two-Layer. The shared rhythm generator becomes one per leg. The simulated heels and toes feed contact information back to the circuit.'),
    (35,33,'OpenSim, Gait2392, and BPAs','Independent dissertation reader','The software, the human model, and the robot muscle routes arrived in the same sentence.',
     'OpenSim is software. Gait2392 is its 23-degree-of-freedom human model with 92 musculotendon actuators. The robot counterpart has relocated routes for muscles that shorten under pressure.'),
    (37,33,'CAD workflow and tools','Independent technical reader','The CAD program, neural block library, simulator, and bridge were named around a figure without distinct roles.',
     'The paragraph above Figure 6.2 gives SOLIDWORKS, SNS-Simscape, Simulink, and the MathWorks bridge separate jobs. The bridge sends commands and returns motion and forces.'),
    (55,51,'Fixed ground and treadmill','Independent technical reader','Table 6.1 placed the earlier W2L pathway removal on the treadmill.',
     'The stationary plane and zero ground displacement in the three saved conditions establish the first setting. The W2L108 evaluator imposes belt motion in the later test. Table 6.1 gives both in order.'),
    (60,56,'BPA20','Independent dissertation reader','The arm row used a block name without identifying its size.',
     'The first BPA20 cell gives the nominal 20 mm diameter. The preceding paragraph has already defined how a BPA acts.'),
    (67,63,'MyoConverter and neural signals','Independent technical reader','Converter credit and the graded and spiking circuits needed an operating description.',
     'The converter translates geometry and paths. Graded neurons carry a continuously varying signal. Spiking neurons emit events at threshold.'),
    (69,65,'Static optimization','Independent technical reader','The source of the reference activations was named before the reader knew what it estimated.',
     'The force comparison now identifies how static optimization estimates muscle load sharing at individual poses. Its citation is shared with Appendix D.'),
    (79,75,'Simscape Multibody','Independent dissertation reader','The second mechanical environment needed a description.',
     'The CAD paragraph identifies the rigid-body and joint library and its relationship to Simulink.'),
    (87,83,'W2Lvar and motor pools','Ben and independent readers','A code-like name and topology counts did not tell the reader which body or muscle organization had changed.',
     'Bolen adapted W2L to Gait2392. Each leg has one pair of hip pattern channels and another for knee and ankle. The same paragraph defines the motor pool the researcher inspects.'),
    (98,94,'Graphical editor','Independent dissertation reader','The command appeared as a tool label in the table.',
     'The first table entry calls SNSGui a circuit editor. The executable example and entry point follow.'),
    (109,105,'LIF, Ia, and Ib','Independent technical reader','A diagram required readers to recognize labels for a neuron and two sensory pathways.',
     'Leaky integrate-and-fire neurons generate a spike at threshold. Ia spindle fibers signal stretch, and Ib tendon-organ fibers signal tension. The definitions precede Figure 6.3.'),
    (118,114,'SADb','Independent dissertation reader','The abbreviation alone could leave a reader outside the database chapter guessing where the experimental evidence was stored.',
     'Sensory Afferent Database appears in full in the sentence that describes the researcher tracing a stimulated fiber and its motor target.'),
    (151,147,'NumPy','Independent technical reader','The name of a numerical library could be read as an additional neural model.',
     'The paragraph identifies NumPy as the Python numerical array library and describes which time histories were compared.'),
    (153,149,'Cross-program W2L test','Independent technical reader','The two W2L implementations and the source of the input needed to remain distinct from the intact moving walker.',
     'The 26-neuron Simulink circuit and 95-neuron reference share imposed heel signals and the same drive in this test. All four neural comparison signals and the earlier reported metrics are preserved.'),
    (174,170,'Sensory return in suspension','Independent dissertation reader','A stage label could be mistaken for a change in mechanical support.',
     'The stage description defines deafferented and afferented by whether the controller receives sensory return, before the table of tests.'),
    (203,199,'BiW2Lvar+ and hunt','Independent technical reader','The controller label did not explain the hip–knee extension of W2Lvar or the relationship to Hunt.',
     'Bolen added pathways for muscles crossing both joints. Hunt supplied the calibration method. The lowercase selector identifies the saved configuration.'),
    (222,218,'s3k and syn6','Ben and independent readers','The first occurrence reported a failure before explaining the controllers or their intended outputs.',
     'The comparison identifies each as a Gait2392 controller. Four pattern cells per leg stagger flexor and extensor signals in s3k. Six channels per leg recruit muscle groups in syn6.'),
    (224,218,'Original timing failure','Independent technical reader','The original grouping could imply s3k had been evaluated against the same acceptance windows as syn6.',
     'The 0/12 timing result belongs to syn6. Its components, target, and stuck channel are identified. The s3k trace had two components per side and no predefined channel windows.'),
    (232,226,'W2L pathway removal','Independent technical reader','The original results paragraph also called the supported surface a treadmill.',
     'The supported, fixed-ground description now matches the plane, evaluator, and three saved trial inputs. The 40 N perturbation, phase, distance, and scoring window are preserved.'),
    (241,235,'W2L108','Independent technical reader','The wording assigned the 108 neurons and 216 connections to an unspecified robot model.',
     'The circuit has 108 neurons and 216 connections. It controlled the twelve-muscle walker on a moving treadmill. The right and left knee values remain attached to that test.'),
    (258,252,'Timing reference','Independent technical reader','“Native SNS” left the timed reference implementation unidentified.',
     'The 1,025-cell Gait2392 comparison used the controller’s NumPy reference implementation in SNS-Toolbox and its sparse counterpart.'),
    (260,256,'Hunt-PT+RC','Independent technical reader','A compact-backend label carried assumptions about both the balance pathway and the recurrent motor pathway.',
     'The center-of-mass balance and Renshaw recurrent-inhibition pathways are described before Table 6.4. Its 1,033-cell comparison uses a separate pair of implementations.'),
    (264,256,'Timing caption','Independent technical reader','The meaning of warm timings was implicit in the earlier caption.',
     'The third row times three initialized blocks and gives the median. The fourth row times a coupled evaluator. The two timing scopes retain their original values.'),
    (298,290,'Transfer experiment','Dissertation editor','The longer discussion pushed the second half of the conclusion onto an otherwise empty page in an intermediate proof.',
     'The two proposed experiments keep the prescribed and closed-loop cases, their observations, and the role of pressure and tendon compliance. The full conclusion now fits on printed page 162.'),
]


def write(name,data):
    (HERE/name).write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')


def main():
    previous = review.read_json(PRIOR/'PASSAGE_REVIEW.json')
    chapters = review.DEFAULT_CHAPTERS+[review.CHAPTERS/'95-AppendixD.tex']
    active = {p.name:{review.digest(s) for _,_,_,s in review.passages(p)} for p in chapters}
    current = {n:s for n,_,_,s in review.passages(chapters[0])}
    retained = [p for p in previous['passages'] if p['after_sha256'] in active[p['file']]]
    superseded = [p['id'] for p in previous['passages'] if p not in retained]
    with ZipFile(PRIOR/'Chapter_6_and_Appendix_D_Overleaf.zip') as z:
        old = z.read('chapters/60-sns.tex').decode().splitlines()
    rows = []
    for i,(line,baseline,title,source,problem,answer) in enumerate(CASES,1):
        rows.append({'id':f'FIRST-{i:02}', 'title':title, 'file':'60-sns.tex','line':line,
                     'source_of_finding':source,'original':old[baseline-1],'prior_delivery_line':baseline,
                     'finding':problem,'proposed_and_applied':current[line],'after_sha256':review.digest(current[line]),
                     'historical_examples':['BEN-NOV18-FIRST-USE','M-022','M-023'],
                     'reviewer_reason':answer,'local_readiness':3,'advisor_score':None,
                     'author_prompted_editor_deduction':-2 if title=='Walker Two-Layer' else None})
    decisions = {'model':previous['model'], 'prompt':previous['prompt']+
                 ' Read a printed chapter as a first-time reader. For each program or model, identify it and its role from the first encounter. For controller variants, identify body, circuit, and experimental condition before interpreting a result. Check captions, diagrams, and table order; verify the code name and its meaning against original model files.',
                 'advisor_event_ids_at_review':previous['advisor_event_ids_at_review'],
                 'baseline':'../ch6_sequence_20261118/Chapter_6_and_Appendix_D_Overleaf.zip',
                 'prior_assistant_chapter_sha256':'523614733eb2041f09bd8e8139cee0fddc605891ca6008a2f9d75c2dc8e453ab',
                 'passages':retained+rows,'superseded_prior_case_ids':superseded,
                 'author_deductions_this_pass':{'BEN-NOV18-FIRST-USE':-2},
                 'previous_author_deductions':previous['editorial_deductions_against_prior_assistant_draft'],
                 'interpretation':'One recurring substantive author observation covered several definitions. Each affected name is reviewed, and the one observation gets one deduction. Readiness 3/4 records source and context verification by the local reviewers.'}
    write('PASSAGE_REVIEW.json',decisions)
    screen = review.make_review(chapters,review.signals(),decisions)
    assert not screen['saved_judgments_stale_or_absent'] and not screen['reopened_judgments']
    assert screen['observed_advisor_signals']==429 and screen['observed_advisor_score']==-551
    write('CHAPTER_SCREEN.json',screen)
    table = (PRIOR/'REVIEW_TABLE.html').read_text(encoding='utf-8')
    old_count = table.count('<tr data-group=')
    headline = ('<h1>Dissertation feedback and revisions</h1><p>The first rows show where Chapter 6 now explains each model and program, how the fixed-ground and moving-belt tests differ, and what each controller timing test measured. Earlier rows retain their dated wording and dispositions. Ben’s first-use observation adds one −2 against the previous assistant draft. The five author marks now total −10. The 429 attributed advisor signals total −551, including the November subtotal of −42.</p>'
                '<p><a href="Chapter_6_terms_reading_pages.pdf">Seven reading pages</a> · <a href="FIRST_USE_REVIEW.md">Changes and evidence</a> · <a href="PASSAGE_REVIEW.json">Passage decisions</a> · <a href="../advisor_feedback_history/README.md">Historical scoring rules</a></p>')
    table = re.sub(r'<h1>.*?</p>\s*<p>.*?</p>(?=\s*<nav>)',lambda _:headline,table,count=1,flags=re.S)
    table = table.replace('<option>Chapter 6 sequence','<option>First-use review</option><option>Chapter 6 sequence',1)
    rendered = []
    for row in rows:
        status = row['reviewer_reason']+' Applied locally, 3/4.'
        if row['author_prompted_editor_deduction']:
            status += ' One author-prompted deduction of −2 for the recurring issue.'
        cols = (row['source_of_finding']+' · '+row['id'], f"60-sns.tex, line {row['line']} · {row['title']}",
                row['original'],row['proposed_and_applied'],row['finding'],status)
        rendered.append('<tr data-group="First-use review">'+''.join('<td>'+html.escape(c)+'</td>' for c in cols)+'</tr>')
    table = table.replace('<tbody>','<tbody>\n'+'\n'.join(rendered),1)
    assert table.count('<tr data-group=') == old_count+len(rows)
    (HERE/'REVIEW_TABLE.html').write_text(table,encoding='utf-8')
    print(len(retained),'earlier decisions preserved,',len(superseded),'superseded,',len(rows),'new decisions,',
          len(screen['passages']),'screened units,',old_count+len(rows),'table rows')


if __name__=='__main__':
    main()
