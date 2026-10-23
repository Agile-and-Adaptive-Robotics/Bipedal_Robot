from pathlib import Path
import subprocess
C=Path(r'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022\Revision_C')
for name in ['pulley_torque_resolved','pulley_margin_resolved']:
    p=C/'figs/30-torque'/name
    subprocess.run([r'D:\MiKTeX\miktex\bin\x64\mgs.exe','-q','-dNOPAUSE','-dBATCH','-dEPSCrop','-dCompatibilityLevel=1.5','-sDEVICE=pdfwrite','-sOutputFile='+str(p.with_suffix('.pdf')),str(p.with_suffix('.eps'))],check=True)
p=C/'review/04_Torque.md';lines=p.read_text().splitlines();p.write_text('\n'.join(dict.fromkeys(lines))+'\n',encoding='utf-8')
print('Normalized the two new surface PDFs to PDF 1.5')
