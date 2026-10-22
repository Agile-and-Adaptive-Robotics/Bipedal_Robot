from pathlib import Path
import subprocess,os,time,json
H=Path(__file__).resolve().parent;P=H.parent.parent/'Documentation/Reports and Papers/Dissertation/ProofFinal'
B=Path(r'D:\temp')/('diss_advisor_figures_'+time.strftime('%Y%m%d_%H%M%S'));(B/'chapters').mkdir(parents=True)
exe=Path(r'D:\MiKTeX\miktex\bin\x64');env=os.environ.copy();env['BIBINPUTS']=str(P)+';';env['BSTINPUTS']=str(P)+';'
latex=[str(exe/'pdflatex.exe'),'-interaction=nonstopmode','-halt-on-error',f'-output-directory={B}','main.tex']
for i,(cmd,cwd) in enumerate([(latex,P),([str(exe/'bibtex.exe'),'main'],B),(latex,P),(latex,P)]):
    with (B/f'pass_{i+1}.txt').open('w',encoding='utf-8') as f:r=subprocess.run(cmd,cwd=cwd,env=env,stdout=f,stderr=subprocess.STDOUT)
    print('pass',i+1,r.returncode,flush=True)
    if r.returncode:raise RuntimeError(B)
(H/'build_path.txt').write_text(str(B));print('BUILD',B)
