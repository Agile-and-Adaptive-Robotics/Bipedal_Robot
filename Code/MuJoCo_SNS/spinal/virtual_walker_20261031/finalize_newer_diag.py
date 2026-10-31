import json
from import_virtual_walker import HERE
from assess_guides import analyze
from render_trial import render
from draw_trial import draw
folder=HERE/'runs/W2L_VirtualWalker_Newer_GUI1nA/newer_v3_diag_v1_ground_guides';out=analyze(folder);out['variant']='Diagnostic tied0/0 tonic, V3 output0, c1 retained, GUI1nA source; not original source configuration.';(folder/'GUIDE_ASSESSMENT.json').write_text(json.dumps(out,indent=2));render(folder,5.,8.2);draw(folder,5.,8.2)
