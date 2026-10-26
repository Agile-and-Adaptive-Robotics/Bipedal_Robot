from pathlib import Path
import subprocess,os,json,re
root=Path(r'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261023\LIVE_ASTRA_20261026')
project=root/'revised_project'; build=project/'build_astra_05'
for source in (root/'after').glob('*_LIVE.tex'):
    target=project/'chapters'/source.name.replace('_LIVE.tex','.tex')
    target.write_bytes(source.read_bytes())
(build/'chapters').mkdir(parents=True,exist_ok=True)
texbin=Path(r'D:\MiKTeX\miktex\bin\x64')
env=os.environ.copy(); env['BIBINPUTS']=str(project)+';'; env['BSTINPUTS']=str(project)+';'
commands=[([str(texbin/'pdflatex.exe'),'-interaction=nonstopmode','-halt-on-error','-output-directory='+str(build),'main.tex'],project),([str(texbin/'bibtex.exe'),'main'],build)]
commands += [commands[0],commands[0]]
results=[]
for i,(cmd,cwd) in enumerate(commands):
    p=subprocess.run(cmd,cwd=cwd,env=env,capture_output=True,timeout=180)
    (build/f'pass_{i}.txt').write_bytes(p.stdout+p.stderr)
    results.append(p.returncode)
    if p.returncode: break
log=(build/'main.log').read_text(errors='replace')
res={'returncodes':results,'errors':re.findall(r'^!.*',log,re.M),'undefined':re.findall(r'[^\n]*undefined[^\n]*',log),'missing':re.findall(r'[^\n]*not found[^\n]*',log),'overfull':re.findall(r'Overfull[^\n]*',log),'pages':re.findall(r'Output written[^\n]*',log)}
(root/'build_result.json').write_text(json.dumps(res,indent=2),encoding='utf-8'); print(json.dumps(res,indent=2))
