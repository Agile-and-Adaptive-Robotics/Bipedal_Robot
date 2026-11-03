from pathlib import Path
import hashlib,json,shutil,zipfile
W=Path(__file__).resolve().parent;B=W/'build/plot_rotation_180'
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def save(p,x):Path(p).write_text(json.dumps(x,indent=2,ensure_ascii=False),encoding='utf-8')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
d=read(W/'build/delivery_checks.json');p=read(B/'package_check.json');v=read(B/'validation.json')
assert d['passed'] and d['review_pdf_pages']==145
assert d['sha256']==p['candidate_sha256']==v['finalSha256']==sha(d['deck'])
with zipfile.ZipFile(d['deck']) as z:
    changed=[n for n in z.namelist() if hashlib.sha256(z.read(n)).hexdigest()!=p['original_part_hashes'][n]]
assert changed==p['changed_parts'] and len(changed)==1
coverage_path=W/'build/visual_review_coverage.json'
shutil.copyfile(coverage_path,B/'prior_visual_review_coverage.json')
c=read(coverage_path);changed_renders=[]
for row in c['slides']:
    n=row['slide'];h=sha(W/f'build/render/slide_{n:03d}.png')
    if row['sha256']!=h:changed_renders.append(n)
    if n in [8,33]:
        row.update(sha256=h,individual_review='Parent editor: original-view comparison and native PowerPoint re-render',final_changed_slide_recheck=True)
assert changed_renders==[8,33]
c.update(deck_sha256=d['sha256'],final_changed_slides=[8,33],unchanged_image_bytes_from_prior_individually_reviewed_final=143)
c['latest_edit_scope']='Only the shared surface image changed. Slides 8 and 33 individually inspected after native rendering; the other 143 slide package contents and existing renders are unchanged.'
save(coverage_path,c)
receipt_path=W/'build/final_delivery_receipt.json'
shutil.copyfile(receipt_path,B/'prior_delivery_receipt.json')
r=read(receipt_path)
edit={'date':'2026-11-02','description':'Rotated both comparison surface views by 180 degrees about their z axes, restoring the original axis orientation. No data, scale, color, slide text, or layout changes.','slides':[8,33],'deck_sha256':d['sha256'],'validation':str(B/'validation.json'),'all_other_package_parts_identical':True,'rendered_and_individually_inspected':[8,33]}
r.update(sha256=d['sha256'],file_bytes=d['bytes'],validation=str(B/'validation.json'),status='v21 updated in place: comparison plots rotated and inspected on slides 8 and 33; prior complete-deck review remains applicable to unchanged slides.')
r.setdefault('in_place_edits',[]).append(edit)
save(receipt_path,r);save(W.parent/'current_deck.json',r)
save(B/'completed.json',edit)
removed=[]
for name in ['artifact_images.pptx','candidate.pptx','prior.pptx']:
    f=(B/name).resolve();assert f.is_relative_to(B.resolve())
    if f.exists():
        removed.append({'path':str(f),'bytes':f.stat().st_size});f.unlink()
save(B/'cleanup.json',removed)
print(json.dumps({'revision':21,'changed_slides':[8,33],'changed_package_parts':changed,'unchanged_slides':143,'sha256':r['sha256'],'temporary_presentations_removed':len(removed)}))
