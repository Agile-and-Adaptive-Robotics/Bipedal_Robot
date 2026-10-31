"""Compose existing native-fall and supported-W2L captures; no simulation."""
from pathlib import Path
import json,hashlib
from PIL import Image
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(__file__).parent
a=Path(r'D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\C3NS_Revision_20261031_v17\assets\free_drop_poster.png')
b=Path(r'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer\canonical108_ground_hip0p6_flex0p5_fresh_v1\SUPPORTED_STEPPING_STILL.png')
crop_a=(45,65,400,310);crop_b=(120,115,790,625)
images=[Image.open(a).convert('RGB').crop(crop_a),Image.open(b).convert('RGB').crop(crop_b)]
plt.rcParams.update({'font.family':'Arial','font.size':10,'pdf.fonttype':42,'ps.fonttype':42})
fig,axes=plt.subplots(1,2,figsize=(7.2,3.7));fig.patch.set_facecolor('white')
for ax,image,title in zip(axes,images,[r'$\mathbf{(A)}$ Initial AnimatLab model',r'$\mathbf{(B)}$ Later MuJoCo W2L']):
 ax.imshow(image);ax.axis('off');ax.set_anchor('N');ax.set_title(title,fontsize=12,loc='left',pad=9)
fig.text(.255,.14,'Root released\nNo vertical support',ha='center',va='top',fontsize=11.6,linespacing=1.3)
fig.text(.745,.14,'Passive rails and 30 mm guides\n108 neurons; 40.5% BW net assistance',ha='center',va='top',fontsize=11.6,linespacing=1.3)
fig.subplots_adjust(left=.025,right=.985,bottom=.22,top=.90,wspace=.10)
fig.savefig(out/'FALL_TO_SUPPORTED.pdf',facecolor='white');fig.savefig(out/'FALL_TO_SUPPORTED.png',dpi=350,facecolor='white');plt.close(fig)
caption='(A) The initial Li–Deng/W2L adaptation after support release in AnimatLab. (B) A later 108-neuron W2L model during alternating ground stepping in MuJoCo, with passive rails and 30 mm guides. The selected run produced nine left and eight right qualified unloaded-clear swings in 10 s, with 40.5% mean net vertical assistance. These are different development and support conditions, rather than a matched comparison or independent-walking demonstration.'
(out/'FALL_TO_SUPPORTED_caption.txt').write_text(caption,encoding='utf-8')
record=dict(machine='easteregg2',panelA=dict(source=str(a),sha256=hashlib.sha256(a.read_bytes()).hexdigest(),crop_pixels=crop_a,condition='Existing poster from Ben\'s free_drop.mp4; source-time not inferred from a different trial.'),panelB=dict(source=str(b),sha256=hashlib.sha256(b.read_bytes()).hexdigest(),crop_pixels=crop_b,time_s=4.84,run='canonical108_ground_hip0p6_flex0p5_fresh_v1',condition='108-neuron graded SNS/Hill; near-flat IK starting pose; hipE60%/hipF30%; original30mm guides/rails−20mm; 40.50%BW whole-trial net assistance.'),composition='Existing pixels cropped to remove original header/empty background; no model image generated or dynamics rerun.',font_floor_pt=10,figure_inches=[7.2,3.7],caption=caption)
(out/'FALL_TO_SUPPORTED_PROVENANCE.json').write_text(json.dumps(record,indent=2),encoding='utf-8');print(out/'FALL_TO_SUPPORTED.pdf')
