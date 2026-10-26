from pathlib import Path
import subprocess,os,time,json,re,shutil
E=Path(r'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261023\Revision_E');V=E/'review/targeted_reconciliation'
B=E/('build_targeted_'+time.strftime('%Y%m%d_%H%M%S'));(B/'chapters').mkdir(parents=True)
exe=Path(r'D:\MiKTeX\miktex\bin\x64');env=os.environ.copy();env['BIBINPUTS']=str(E)+';';env['BSTINPUTS']=str(E)+';'
base=[str(exe/'pdflatex.exe'),'-interaction=nonstopmode','-halt-on-error',f'-output-directory={B}']
steps=[(base+['-draftmode','main.tex'],E),([str(exe/'bibtex.exe'),'main'],B),(base+['-draftmode','main.tex'],E),(base+['main.tex'],E)]
for i,(cmd,cwd) in enumerate(steps,1):
 with (B/f'pass_{i}.log').open('w',encoding='utf-8') as out:r=subprocess.run(cmd,cwd=cwd,env=env,stdout=out,stderr=subprocess.STDOUT)
 print('pass',i,r.returncode,flush=True)
 if r.returncode:raise RuntimeError(str(B))
log=(B/'main.log').read_text(errors='replace')
j=dict(build=str(B),errors=len(re.findall(r'^!',log,re.M)),undefined_references=len(re.findall(r'undefined references|Reference .* undefined|Citation .* undefined',log)),pages=re.findall(r'Output written on .*?\((\d+) pages',log,re.S),overfull_boxes=len(re.findall(r'Overfull \\hbox',log)))
(V/'build_result.json').write_text(json.dumps(j,indent=2),encoding='utf-8')
shutil.copy2(B/'main.pdf',E/'Revision_E_reference_reconciled.pdf')
print(json.dumps(j),flush=True)
