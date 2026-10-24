from pathlib import Path
import re,json,difflib,subprocess,os,time
R=Path(r'D:\Github\Bipedal_Robot'); D=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D'
conservation={}; diff=''
for file in ['10-introduction.tex','15-background.tex','50-database.tex','60-sns.tex']:
 b=(D/'review/before'/file).read_text(encoding='utf-8'); a=(D/'chapters'/file).read_text(encoding='utf-8')
 diff+=''.join(difflib.unified_diff(b.splitlines(True),a.splitlines(True),fromfile='before/'+file,tofile='after/'+file))
 patterns={'labels':r'\\label\{[^}]*\}','citations':r'\\cite\w*(?:\[[^\]]*\])*\{[^}]*\}','figure_assets':r'\\includegraphics(?:\[[^\]]*\])?\{[^}]*\}','figure_count':r'\\begin\{figure\}'}
 conservation[file]={k:re.findall(p,b,re.S)==re.findall(p,a,re.S) for k,p in patterns.items()}
 # Recorded numerical values remain: position of the 92-actuator explanation changes; list ordinals become enumerate.
 nums=lambda s: sorted(re.findall(r'(?<![A-Za-z])\d+(?:\.\d+)?',re.sub(r'\\ref\{[^}]*\}|\\cite\w*(?:\[[^\]]*\])*\{[^}]*\}','',s)))
 from collections import Counter
 conservation[file]['removed_numeric_tokens']=list((Counter(nums(b))-Counter(nums(a))).elements())
 conservation[file]['added_numeric_tokens']=list((Counter(nums(a))-Counter(nums(b))).elements())
 assert all(conservation[file][k] for k in patterns),conservation[file]
(D/'review/conservation.json').write_text(json.dumps(conservation,indent=2),encoding='utf-8'); (D/'review/changes.diff').write_text(diff,encoding='utf-8')
B=D/('build_'+time.strftime('%Y%m%d_%H%M%S')); (B/'chapters').mkdir(parents=True)
exe=Path(r'D:\MiKTeX\miktex\bin\x64');env=os.environ.copy();env['BIBINPUTS']=str(D)+';';env['BSTINPUTS']=str(D)+';'
cmd=[str(exe/'pdflatex.exe'),'-interaction=nonstopmode','-halt-on-error',f'-output-directory={B}','main.tex']
(D/'review/latest_build.txt').write_text(str(B),encoding='utf-8')
for i,(command,cwd) in enumerate([(cmd,D),([str(exe/'bibtex.exe'),'main'],B),(cmd,D),(cmd,D)]):
 with (B/f'pass_{i+1}.log').open('w',encoding='utf-8') as log:r=subprocess.run(command,cwd=cwd,env=env,stdout=log,stderr=subprocess.STDOUT)
 print('pass',i+1,r.returncode,flush=True)
 if r.returncode: raise RuntimeError('Build failed: '+str(B))
log=(B/'main.log').read_text(errors='replace')
result={'build':str(B),'errors':len(re.findall(r'^!',log,re.M)),'undefined_references':len(re.findall(r'undefined references|Reference .* undefined|Citation .* undefined',log)),'missing_figures':len(re.findall(r'File .* not found',log)),'pages':re.findall(r'Output written on .*?\((\d+) pages',log,re.S),'overfull_boxes':len(re.findall(r'Overfull \\hbox',log))}
(D/'review/build_result.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result),flush=True)
