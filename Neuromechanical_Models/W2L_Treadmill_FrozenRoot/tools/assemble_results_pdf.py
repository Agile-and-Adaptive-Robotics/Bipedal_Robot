"""Collect the native screenshots and vector result figures in one review PDF."""
from pathlib import Path
from io import BytesIO
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader
from pypdf import PdfReader,PdfWriter
R=Path.cwd();p=R/'Neuromechanical_Models/W2L_Comparison_20261030';out=p/'output';out.mkdir(exist_ok=True)
b=BytesIO();c=canvas.Canvas(b,pagesize=(792,612));c.setTitle('W2L native AnimatLab comparison')
c.setFont('Helvetica-Bold',20);c.drawString(30,578,'W2L: native AnimatLab comparison')
c.setFont('Helvetica',11);c.drawString(30,556,'Original and saved bilateral-RG controllers. Matched mechanics. No controller retuning.')
c.drawString(30,539,'104 successful ten-second trials. Air contact counts verified zero. Both free models fall.')
for x,name,title in [(30,'legacy_drop_clean.png','Free Root, t = 2.141 s'),(411,'legacy_partial50_after_release.png','50% upward force, t = 2.226 s')]:
 c.setFont('Helvetica-Bold',12);c.drawString(x,512,title);c.drawImage(ImageReader(str(p/'screenshots'/name)),x,308,width=351,height=190,preserveAspectRatio=True,anchor='c')
c.setFont('Helvetica-Bold',12);c.drawString(30,279,'First time pelvis height falls below 0.5 m')
c.setFont('Helvetica',11)
for y,t in [(260,'Original controller: 1.103 s free; 1.790 s with 50% upward force.'),(243,'Newer controller: 0.999 s free; 2.077 s with 50% upward force.'),(219,'Partial support = 205.4031 N positive-Y force at Root, 50% of declared bodyweight.'),(202,'There is no balance constraint or rotational stabilization. No unsupported walking is demonstrated.'),(174,'Pages 2-17: motor-neuron voltages and joint angles for both controllers in all four conditions.'),(157,'Pages 18-25: Deng 2019 PF/RG pulse comparisons. Page 26: Deng 2022 phase-response transfer.'),(140,'Page 27: all-condition pelvis-height comparison. Joint angles use the native AnimatLab convention.'),(108,'The literature pulse protocols were transferred to W2L. These are not full rat-model replications.')]:c.drawString(30,y,t)
c.setFont('Helvetica',10);c.drawString(30,84,'Deng et al. (2019), Biomimetics 4(1):21. https://www.mdpi.com/2313-7673/4/1/21');c.linkURL('https://www.mdpi.com/2313-7673/4/1/21',(30,80,590,94),relative=0)
c.drawString(30,67,'Deng et al. (2022), Biomimetics 7(4):226. https://www.mdpi.com/2313-7673/7/4/226');c.linkURL('https://www.mdpi.com/2313-7673/7/4/226',(30,63,610,77),relative=0)
c.drawString(30,41,'Full protocols, source snapshots, signals, and recording projects: W2L_Comparison_20261030/README.txt')
c.showPage();c.save();b.seek(0)
w=PdfWriter();w.append(PdfReader(b));paths=[]
for model in ['legacy','bilateral']:
 for cond in ['air','treadmill','drop','partial50']:
  for kind in ['MN','joints']:paths.append(p/'figures'/f'{model}_{cond}_{kind}.pdf')
for model in ['legacy','bilateral']:
 for cond in ['air','treadmill','drop','partial50']:paths.append(p/'figures'/f'{model}_{cond}_deng2019.pdf')
paths.extend([p/'figures/deng2022_phase_response.pdf',p/'figures/support_comparison.pdf'])
for f in paths:
 page=len(w.pages);w.append(str(f),import_outline=False);w.add_outline_item(f.stem.replace('_',' '),page)
dest=out/'W2L_Native_AnimatLab_Results.pdf';w.write(dest);assert len(PdfReader(dest).pages)==27
print(str(dest.resolve()))
