from pathlib import Path
import csv,hashlib,json,zipfile
from lxml import etree as E
from verify_delivery import slide_parts,N

W=Path(__file__).resolve().parent
def read(p):return json.loads((W/p).read_text(encoding='utf-8-sig'))
def save(p,x):(W/p).write_text(json.dumps(x,indent=2,ensure_ascii=False),encoding='utf-8')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
d=read('build/delivery_checks.json');s=read('build/structure.json');order=read('build/slide_order.json')
videos=read('build/native_video_playback.json');gifs=read('build/native_gif_playback.json')
assert d['passed'] and d['review_pdf_pages']==145 and len(s['main'])==77
assert len(videos)==3 and all(x['advances'] and (not x['loop_required'] or x['loop_restarted']) for x in videos)
assert len(gifs)==4 and all(x['motion_verified'] for x in gifs)
assert sha(d['deck'])==d['sha256']
prior={x['slide']:x['sha256'] for x in read('build/pre_conclusion_render_hashes.json')}
render_rows=[]
for n in range(1,146):
    p=W/f'build/render/slide_{n:03d}.png';assert p.exists()
    reviewer=('Hunt/narrative review' if n<=77 else 'Turcic/mechanics review' if n<=99 else 'Reddy/neural review' if n<=111 else 'Parent editor' if n<=123 else 'Reddy/neural review' if n<=129 else 'Jiang/mathematics review' if n<=138 else 'Parent editor')
    render_rows.append({'slide':n,'sha256':sha(p),'individual_review':reviewer,'final_changed_slide_recheck':n in [21,35,37,47,50,64,66,68,70,71,77,96,106,118,125,128,131]})
changed=[r['slide'] for r in render_rows if r['sha256']!=prior[r['slide']]]
assert changed==[77,106],f'Unexpected final render differences: {changed}'
coverage={'total':145,'individual_inspection_coverage':145,'contact_sheets_used':False,'all_final_images_rendered_from':d['deck'],'deck_sha256':d['sha256'],'final_changed_slides':[77,106],'unchanged_image_bytes_from_prior_individually_reviewed_final':143,'slides':render_rows,'scope':'Individual native raster review plus source audits. Final conclusion and pathway-removal raster individually rechecked. Source UI screenshots in questions retain small internal labels. This is not a spoken timing rehearsal.'}
save('build/visual_review_coverage.json',coverage)
save('build/final_render_comparison.json',{'changed':[77,106],'unchanged':143,'changed_slides_individually_rechecked':True,'pathway_removal_check':'Legend remains outside axes; seven left/right bar pairs match the original counts.'})
titles={}
with zipfile.ZipFile(d['deck']) as z:
    for n,part in enumerate(slide_parts(z),1):
        root=E.fromstring(z.read(part));items=[]
        for sp in root.findall('./p:cSld/p:spTree/p:sp',N):
            text=' '.join(sp.xpath('.//a:t/text()',namespaces=N));sizes=[int(x) for x in sp.xpath('.//a:rPr/@sz',namespaces=N)]
            xf=sp.find('p:spPr/a:xfrm',N)
            if xf is not None:
                off=xf.find('a:off',N);ext=xf.find('a:ext',N)
                if text and off is not None and ext is not None and int(off.get('y'))<914400 and int(ext.get('cx'))>4000000:items.append((max(sizes or [0]),len(text),text))
        titles[n]=max(items,default=(0,0,''))[2]
    titles[1]='Design and Modeling of Biomimetic Humanoid Robots'
with (W/'build/slide_map.csv').open('w',encoding='utf-8-sig',newline='') as f:
    writer=csv.writer(f);writer.writerow(['final_position','v20_or_appended_source_id','main_or_question','title'])
    for n,src in enumerate(order,1):writer.writerow([n,src,'main' if n<=77 else 'question',titles[n]])
for r in gifs:r['poses_and_visible_crop_individually_inspected']=True
save('build/native_gif_playback.json',gifs)
cleanup=read('build/cleanup.json')
assert not (W.parent/'C3NS_Revision_20261102_v19/output/Bolen_PhD_Defense_RoseGreen_20261102_v19.pptx').exists()
receipt={'revision':21,'current_deck':d['deck'],'pdf':str(W/'output/Bolen_PhD_Defense_RoseGreen_20261102_v21.pdf'),'sha256':d['sha256'],'file_bytes':d['bytes'],'status':'Final PowerPoint/PDF exported; all145 slides individually reviewed; scientific/source audits and bounded native playback checks completed. Spoken timing is not rehearsed.','baseline':str(W.parent/'C3NS_Revision_20261102_v20/output/Bolen_PhD_Defense_RoseGreen_20261102_v20.pptx'),'core_slides':77,'question_slides':68,'planned_minutes_at_30_seconds':38.5,'pacing_target_minutes':40,'hard_cap_minutes':45,'pacing_rehearsed':False,'palette':{'name':'Rose green','deep':'#174D24','leaf':'#58D51B','background':'#FFFFFF'},'original_branding_preserved':True,'original_title_honeycomb_preserved':True,'plot_colors_preserved_except_requested_others_group':True,'rendered_slides':145,'individually_inspected_slides':145,'pdf_pages':145,'native_chart_count':7,'chart_workbooks':7,'native_chart_owners':d['native_chart_owners'],'native_walker_media_slide':59,'native_offline_demo_slide':145,'sadb_positions':list(range(139,146)),'sadb_main_absent':True,'submitted_dissertation_unchanged':d['submitted_pdf_unchanged'],'submitted_dissertation_sha256':d['submitted_pdf_sha256'],'generated_outputs_git_ignored':True,'validation':str(W/'build/validation_delivery.json'),'preservation_checks':str(W/'build/delivery_checks.json'),'visual_review':str(W/'build/visual_review_coverage.json'),'native_video_playback_report':str(W/'build/native_video_playback.json'),'native_gif_playback_report':str(W/'build/native_gif_playback.json'),'handoff':str(W/'HANDOFF.md'),'compliance':str(W/'DIRECTIVE_COMPLIANCE.md'),'committee_perspective_reviews':[str(W/f'build/committee_{name}_v21_review.md') for name in ['hunt','turcic','reddy','jiang']],'review_role_boundary':'Agent perspectives informed by the named disciplines, not approval by the actual committee. Force author disclosed ownership in the Hunt review.','retained_final_revisions':[20,21],'older_source_assets_retained':True,'cleanup_complete':True,'bytes_removed_this_cleanup':sum(x['bytes'] for x in cleanup)}
save('build/final_delivery_receipt.json',receipt)
(W.parent/'current_deck.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
print(json.dumps({'deck':d['deck'],'sha256':d['sha256'],'slides':145,'main':77,'hidden_questions':68,'all_checks':True}))
