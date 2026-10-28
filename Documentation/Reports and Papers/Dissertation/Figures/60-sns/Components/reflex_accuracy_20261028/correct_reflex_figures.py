from pathlib import Path
import fitz,shutil,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,FancyArrowPatch
r=Path('D:/Github/Bipedal_Robot');d=r/'Documentation/Reports and Papers/Dissertation';f=d/'Figures/60-sns';a=d/'Archive/reflex_accuracy_20261028';a.mkdir(exist_ok=True);comp=f/'Components/reflex_accuracy_20261028';comp.mkdir(exist_ok=True)
for n in ['knee_reflex_current_layout','beverage_current_layout']:
 if not (a/(n+'.pdf')).exists():shutil.copy2(f/(n+'.pdf'),a/(n+'.pdf'))

if not (a/'60-sns.tex').exists():shutil.copy2(d/'ProofFinal/chapters/60-sns.tex',a/'60-sns.tex')
def circuit(name,labels,title):
 fig,ax=plt.subplots(figsize=(6.53,3));fig.subplots_adjust(0,0,1,1);ax.set(xlim=(0,2),ylim=(0,1));ax.axis('off')
 ax.text(.02,.96,title,fontsize=12,weight='bold',va='top')
 for j,(own,ant,mn) in enumerate(labels):
  xs=[j+.12,j+.48,j+.84];ys=[.65]*3;target=(j+.48,.24)
  for x,label,exc in zip(xs,['Ia '+own,'Ib '+own,'Ia '+ant],[True,False,False]):
   ax.text(x,.79,label,ha='center',fontsize=11)
   ax.add_patch(Circle((x,.65),.035,fill=False,lw=1.1))
   start=(x,.615);end=(target[0]+(x-target[0])*.13,.29)
   ax.add_patch(FancyArrowPatch(start,end,arrowstyle='-|>' if exc else '-',mutation_scale=12,color='#0072B2',linestyle='--',lw=1.3))
   if not exc:ax.plot(*end,'o',color='#0072B2',ms=5)
  ax.text(*target,mn,ha='center',va='center',fontsize=11,bbox=dict(boxstyle='round,pad=.3',fc='white',ec='black'))
 ax.text(1,.04,'Arrowhead: excitation     Filled circle: inhibition',ha='center',fontsize=11)
 fig.savefig(comp/(name+'.pdf'));plt.close(fig)
circuit('knee_circuit',[('ext','flex','MN ext'),('flex','ext','MN flex')],'C  Implemented knee feedback connections')
circuit('arm_circuit',[('biceps','triceps','MN biceps'),('triceps','biceps','MN triceps')],'C  Implemented elbow feedback connections')
def make(n,kind):
 src=fitz.open(a/(n+'.pdf'));
 if kind=='knee':
  for block in src[0].get_text('dict')['blocks']:
   for line in block.get('lines',[]):
    for span in line['spans']:
     if 'Matched-implementation response' in span['text']:src[0].add_redact_annot(fitz.Rect(span['bbox']),fill=(1,1,1))
  src[0].apply_redactions(images=0,graphics=0)
 out=fitz.open();p=out.new_page(width=470,height=650 if kind=='arm' else 600)
 if kind=='arm':
  # Retain mechanics and kinematic data; label in reading order.
  p.show_pdf_page(fitz.Rect(0,16,470,190),src,0,clip=fitz.Rect(0,13,494,173))
  p.insert_text((8,12),'A1  Empty cup (start)',fontsize=12);p.insert_text((252,12),'A2  Loaded cup (end)',fontsize=12)
  p.show_pdf_page(fitz.Rect(10,210,460,423),src,0,clip=fitz.Rect(0,298,260,557))
  p.insert_text((8,205),'B  Feedback on/off: forearm elevation',fontsize=12)
  c=fitz.open(comp/'arm_circuit.pdf');p.show_pdf_page(fitz.Rect(0,427,470,643),c,0)
 else:
  p.show_pdf_page(fitz.Rect(5,15,465,300),src,0,clip=fitz.Rect(0,0,507,317))
  p.insert_text((8,12),'A  Startup under constant load',fontsize=12)
  p.show_pdf_page(fitz.Rect(20,320,190,592),src,0,clip=fitz.Rect(65,337,240,595))
  p.insert_text((8,317),'B  Antagonist-BPA knee',fontsize=12)
  c=fitz.open(comp/'knee_circuit.pdf');p.show_pdf_page(fitz.Rect(194,340,470,575),c,0)
 # Keep circuit labels legible: knee circuit gets its own full-width lower strip instead.
 if kind=='knee':
  out.close();out=fitz.open();p=out.new_page(width=470,height=740)
  p.show_pdf_page(fitz.Rect(5,15,465,300),src,0,clip=fitz.Rect(0,0,507,317));p.insert_text((8,12),'A  Startup under constant load',fontsize=12)
  p.show_pdf_page(fitz.Rect(85,325,385,510),src,0,clip=fitz.Rect(65,337,240,595));p.insert_text((8,317),'B  Antagonist-BPA knee',fontsize=12)
  c=fitz.open(comp/'knee_circuit.pdf');p.show_pdf_page(fitz.Rect(0,518,470,734),c,0)
 temp=comp/(n+'.pdf');out.save(temp);out[0].get_pixmap(matrix=fitz.Matrix(1.2,1.2)).save(a/(n+'_corrected.png'));
 try:shutil.copy2(temp,f/(n+'.pdf'))
 except OSError:
  old=f/(n+'.pdf');old.rename(a/(n+'_locked.pdf'));shutil.copy2(temp,old)
make('beverage_current_layout','arm');make('knee_reflex_current_layout','knee')
p=d/'ProofFinal/chapters/60-sns.tex';s=p.read_text(encoding='utf-8').replace('feedback-on/off kinematics (A), empty and loaded cup positions (B1, B2), and reflex pathways (C)','empty and loaded cup positions (A1, A2), feedback-on/off kinematics (B), and all six implemented feedback connections (C)')
s=s.replace('and sensory pathways (C).','and all six implemented feedback connections (C).')
s=s.replace('The reduced sag demonstrates', 'As the forearm sags, biceps stretch feedback excites the biceps motoneuron and reciprocally inhibits the triceps motoneuron. The separate Ib pathways provide autogenic inhibition of each muscle; they are not an external-load command. The reduced sag demonstrates')
s=s.replace('As the forearm sags, biceps stretch feedback excites the biceps motoneuron and reciprocally inhibits the triceps motoneuron. The separate Ib pathways provide autogenic inhibition of each muscle; they are not an external-load command. '*2, 'As the forearm sags, biceps stretch feedback excites the biceps motoneuron and reciprocally inhibits the triceps motoneuron. The separate Ib pathways provide autogenic inhibition of each muscle; they are not an external-load command. ');p.write_text(s,encoding='utf-8');print('Two corrected vector figures and captions saved.')
