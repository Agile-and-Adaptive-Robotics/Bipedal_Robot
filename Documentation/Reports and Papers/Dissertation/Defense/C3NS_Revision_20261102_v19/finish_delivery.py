from pathlib import Path
import json, hashlib
W=Path(__file__).resolve().parent
pointer=W.parent/'current_deck.json'
state=json.loads(pointer.read_text(encoding='utf-8'))
checks=json.loads((W/'build/delivery_checks.json').read_text())
state.update({
 'revision':19,'current_deck':checks['deck'],'pdf':str(W/'output/Bolen_PhD_Defense_RoseGreen_20261102_v19.pdf'),
 'sha256':checks['sha256'],'file_bytes':checks['bytes'],
 'baseline':str(W.parent/'C3NS_Revision_20261031_v18/output/Bolen_PhD_Defense_RoseGreen_20261031_v18.pptx'),
 'baseline_sha256':'d7e0b9521e94ca391ff8470fb565ff9823acbcf2d0487ff3a5c122546fde738c',
 'status':'v19 delivered; advisor narrative and independent scientific/visual review completed',
 'submitted_dissertation':r'C:\Users\Ben Bolen\Documents\ME603 Dissertation\Defense\Bolen_Dissertation.pdf',
 'planned_minutes':35,'pacing_rehearsed':False,'core_slides':70,'question_slides':53,
 'retained_final_revisions':[18,19],'rendered_slides':123,'individually_inspected_slides':123,
 'contact_sheets_inspected':0,'pdf_pages':123,
 'validation':'C3NS_Revision_20261102_v19/build/validation_delivery.json',
 'native_playback_verification':'C3NS_Revision_20261102_v19/build/native_video_playback.json',
 'visual_review':['C3NS_Revision_20261102_v19/build/supervisor_review.md','C3NS_Revision_20261102_v19/build/mechanics_review.md','C3NS_Revision_20261102_v19/build/neural_review.md'],
 'handoff':'C3NS_Revision_20261102_v19/HANDOFF.md',
 'github_100mb_note':'Final deck is 94,941,742 bytes and remains Git-ignored; no review ZIP created.',
 'latest_reference_commit':'37b558f1','latest_reference_includes_uncommitted_nov2_handoff_update':True,
 'cleanup_complete':True,
})
pointer.write_text(json.dumps(state,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
refs=[W.parent/'SESSION1_NEURAL_CONTROLLER_HANDOFF.md',Path(r'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\GAIT_DISSERTATION_RECIPE_REVIEW.md'),Path(r'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\DEFENSE_DESIGN_DECISIONS.md')]
(W/'build/reference_checkpoint.json').write_text(json.dumps([{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size} for p in refs],indent=2))
print('Current deck points to revision19. Reference hashes recorded.')
