from pathlib import Path
import re,json
from PIL import Image
C=Path(r'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022\Revision_C');S=C/'review/user_sweep_20261023'
for rel in ['15-background/festo_bpas_inflated_deflated.jpg','30-torque/cam_extensor_hardware.jpg']:
 im=Image.open(C/'figs'/rel);print(rel,im.size,im.getexif().get(274))
 p=C/'chapters'/('15-background.tex' if rel.startswith('15') else '30-torque.tex');s=p.read_text(encoding='utf-8')
 # PDFTeX does not rotate JPEG pixels from the camera's EXIF orientation.
 if im.getexif().get(274)==6:
  s=s.replace('keepaspectratio]{'+Path(rel).name,'angle=-90,origin=c,keepaspectratio]{'+Path(rel).name)
 elif im.getexif().get(274)==8:s=s.replace('keepaspectratio]{'+Path(rel).name,'angle=90,origin=c,keepaspectratio]{'+Path(rel).name)
 p.write_text(s,encoding='utf-8')
# Split existing combined check/limit cells without asserting answers to Ben's ZCode questions.
p=C/'chapters/60-sns.tex';s=p.read_text(encoding='utf-8')
s=s.replace(r'\begin{tabularx}{\textwidth}{@{}>{\raggedright\arraybackslash}p{0.23\textwidth}>{\raggedright\arraybackslash}p{0.34\textwidth}>{\raggedright\arraybackslash}X@{}}',r'\begin{tabularx}{\textwidth}{@{}>{\raggedright\arraybackslash}p{0.19\textwidth}*{3}{>{\raggedright\arraybackslash}X}@{}}')
s=s.replace(r'\textbf{Tool} & \textbf{Implemented capability} & \textbf{Verification and limit}',r'\textbf{Model/interface} & \textbf{Capability} & \textbf{Verification} & \textbf{Limit}')
s=s.replace('Actuator geometry and dynamic signs checked; supported closed-loop evaluations, not unsupported walking','Actuator geometry and dynamic signs checked; supported closed-loop evaluations & Unsupported walking not demonstrated')
s=s.replace('Reflex demonstrations and controller comparisons; generic BPA demonstration parameters remain provisional','Reflex demonstrations and controller comparisons & Generic BPA demonstration parameters remain provisional')
s=s.replace('Rhythm-core agreement across NumPy and Simulink; complete Li--Deng contact-feedback walking remains unverified','Rhythm-core agreement across NumPy and Simulink & Complete Li--Deng contact-feedback walking remains unverified')
s=s.replace('Mechanism import achieved; standalone reflex demonstrations do not establish control of the imported models','Mechanism import achieved & Standalone reflex demonstrations do not establish control of imported models')
p.write_text(s,encoding='utf-8')
print('Separated verification/limits; implementation questions remain assigned to ZCode')
