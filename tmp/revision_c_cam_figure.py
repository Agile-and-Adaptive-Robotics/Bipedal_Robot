from pathlib import Path
import shutil,hashlib,json,re,fitz
R=Path(r'D:\Github\Bipedal_Robot');C=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C'
p=C/'chapters/30-torque.tex';s=p.read_text();s=s.replace(r'\includegraphics[width=0.98\textwidth]{ExtensorRoute_redesigned.pdf}',r'\includegraphics[width=0.98\textwidth]{extensor_cam_route.pdf}')
for m in re.finditer(r'\\begin\{figure\}.*?\\end\{figure\}',s,re.S):
    if r'\label{fig:ExtensorRoute_redesigned}' not in m.group():continue
    b=m.group();st=b.index(r'\caption');en=b.index(r'\label',st)
    b=b[:st]+r'''\caption[Cam-extensor route geometry]{Cam-extensor route in the femur frame, from \ang{-120} flexion to \ang{10} extension. Pose panels show the optimized path, cam contacts, and clearance envelopes. Squares mark femoral origins, diamonds tibial insertions, and dashed curves bend-radius checks. The legend records the contact-release sequence.}
'''+b[en:];s=s[:m.start()]+b+s[m.end():];break
s=s.replace('The extensor route has since been redesigned once more by allowing the tibial-side route points to trace a cam profile, and this cam redesign is the current extensor design.','Allowing the tibial-side route points to trace a cam profile produced the current extensor design (Figure~\\ref{fig:ExtensorRoute_redesigned}).')
p.write_text(s,encoding='utf-8')
source=R/'Documentation/Reports and Papers/Knee_Torque_Test/Figures/ExtensorRoute_redesigned.eps';shutil.copy2(source,C/'figs/30-torque/extensor_cam_route.eps')
p=C/'review/torque_figure_sources.json';j=json.loads(p.read_text());j['extensor_cam_route.eps']={'source':str(source),'sha256':hashlib.sha256(source.read_bytes()).hexdigest()};p.write_text(json.dumps(j,indent=2))
p=C/'review/04_Torque.md';s=p.read_text();s+='The main route-geometry figure now uses the companion paper\'s current cam-route EPS (explicit cam-drum radius and contacts), rather than the older direct-route PDF. Its caption was shortened and the cam paragraph now points to it. The original direct-route figure remains in the supporting appendix.\n';p.write_text(s,encoding='utf-8')
d=fitz.open(C/'figs/30-torque/extensor_cam_route.pdf');d[0].get_pixmap(matrix=fitz.Matrix(1,1)).save(C/'review/extensor_cam_route.png')
print('Updated main mechanical route figure to the existing cam export')
