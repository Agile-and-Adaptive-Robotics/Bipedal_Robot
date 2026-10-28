from pathlib import Path
import fitz,shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,FancyArrowPatch
D=Path('Documentation/Reports and Papers/Dissertation');A=D/'Archive/knee_and_robot_update_20261028';comp=D/'Figures/60-sns/Components/knee_transient_reflex_20261028'
s=Path('tmp/correct_reflex_figures.py').read_text();fn=s[s.index('def circuit('):s.index("circuit('knee_circuit'")].replace('fontsize=11','fontsize=14').replace('fontsize=12','fontsize=14');exec(fn)
circuit('knee_circuit',[('ext','flex','MN ext'),('flex','ext','MN flex')],'C  Implemented knee feedback connections')
old=fitz.open(A/'previous_knee_layout.pdf');panel=fitz.open(comp/'transient_panel.pdf');circ=fitz.open(comp/'knee_circuit.pdf');out=fitz.open();p=out.new_page(width=470,height=600)
p.insert_text((8,14),'A  Transient sensory-current response',fontsize=13);p.show_pdf_page(fitz.Rect(0,20,470,270),panel,0)
p.insert_text((8,280),'B  Antagonist-BPA knee (illustrative pose)',fontsize=13);p.insert_image(fitz.Rect(70,285,400,405),stream=old[0].get_pixmap(matrix=fitz.Matrix(4,4),clip=fitz.Rect(85,325,385,510),alpha=False).tobytes('png'));p.show_pdf_page(fitz.Rect(0,410,470,600),circ,0)
p.get_pixmap(matrix=fitz.Matrix(1.3,1.3)).save(A/'knee_transient_current_layout_final.png');out.save(comp/'knee_transient_current_layout_final.pdf');shutil.copy2(comp/'knee_transient_current_layout_final.pdf',D/'Figures/60-sns/knee_transient_current_layout.pdf')
shutil.copy2('tmp/plot_knee_transient_supplement.py',D/'Figures/95-AppendixD/Components/knee_transient_reflex_20261028/plot_knee_transient_supplement.py')
p=D/'ProofFinal/chapters/60-sns.tex';s=p.read_text(encoding='utf8').replace('antagonist-BPA leg model (B), and the six','antagonist-BPA leg schematic in an illustrative pose (B), and the six');p.write_text(s,encoding='utf8')
