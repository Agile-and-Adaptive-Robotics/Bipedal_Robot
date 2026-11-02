from pathlib import Path
import json,hashlib

W=Path(__file__).resolve().parent
def read(name):return json.loads((W/name).read_text(encoding='utf-8-sig'))
pass2={x['slide']:x['sha256'].lower() for x in read('build/pass2_all_hashes.json')}
pass3={x['slide']:x['sha256'].lower() for x in read('build/pass3_all_hashes.json')}
changed=[3,9,23,38,39,42,44,46,47,50,57,74,82,105,127]
records=[]
for n in range(1,128):
    path=W/f'build/render/slide_{n:03d}.png'
    sha=hashlib.sha256(path.read_bytes()).hexdigest()
    assert sha==pass3[n]
    reviewer='supervisor' if n<=12 else 'Turcic-perspective reviewer' if n<=53 else 'Reddy-perspective reviewer' if n<=105 else 'root'
    checked={'slide':n,'sha256':sha,'image':str(path),'pass2_individual_reviewer':reviewer,'changed_from_pass2':sha!=pass2[n]}
    assert checked['changed_from_pass2']==(n in changed)
    if n in changed:checked['pass3_individually_rechecked_by']=['root']
    if n in [38,39,42,44,46,74,82]:checked['pass3_individually_rechecked_by'].append('supervisor')
    if n in [57,74,82,105]:checked['pass3_individually_rechecked_by'].append('Reddy-perspective reviewer')
    records.append(checked)
report={'final_slide_count':127,'individual_inspection_coverage':127,'main':74,'hidden_questions':53,'final_export_images_match_reviewed_pass3':True,'changed_images_reinspected':changed,'method':'All second-pass images individually inspected; every changed third-pass image individually reinspected; all final images verified byte-identical to the inspected third pass. No contact-sheet-only review.','native_geometry_warnings_disposition':{'36':'Leader terminates at the apparatus bracket; no text obscured in native render.','41':'All table rows, units, and source fit above the footer in native render.','60':'Text boxes intersect geometrically but rendered title, model qualifier, and phase labels do not obscure each other.'},'slides':records}
(W/'build/visual_review_coverage.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print('Recorded individual visual-review coverage for all127 final slides.')
