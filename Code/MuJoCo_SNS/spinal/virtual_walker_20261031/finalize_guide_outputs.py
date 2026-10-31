import json
from import_virtual_walker import HERE
from assess_guides import analyze
from render_trial import render
from draw_trial import draw
for filename in ('GUIDE_STUDY_V1.json','GUIDE_EXTREMES_V1.json'):
 path=HERE/'runs'/filename;rows=json.loads(path.read_text())
 for row in rows:row['assessment']=analyze(row['assessment']['source'])
 path.write_text(json.dumps(rows,indent=2))
folder=HERE/'runs/W2L_VirtualWalker/guides_v1_t2_c0.01_h0_p5';render(folder);draw(folder)
