from pathlib import Path
W=Path(__file__).resolve().parent
P=W.parent/'C3NS_Revision_20261102_v19'
for name in ['render_review.ps1','verify_delivery.py','check_video_playback.ps1']:
 text=(P/name).read_text(encoding='utf-8').replace('20261102_v19','20261102_v20')
 if name=='verify_delivery.py':text=text.replace('20261031_v18','20261102_v19')
 (W/name).write_text(text,encoding='utf-8')
text=(P/'finalize.mjs').read_text(encoding='utf-8').replace('20261102_v19','20261102_v20').replace('20261031_v18','20261102_v19')
(W/'finalize.mjs').write_text(text,encoding='utf-8')
