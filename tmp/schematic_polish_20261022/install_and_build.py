from pathlib import Path
import hashlib,json,shutil,subprocess,os,re,time

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent.parent
D=ROOT/'Documentation/Reports and Papers/Dissertation';F=D/'Figures';P=D/'ProofFinal'
base=json.loads((HERE/'baseline_hashes.json').read_text())
jobs=json.loads((HERE/'staged/jobs.json').read_text())
for ch,n in jobs:
    for p in [F/ch/(n+'.pdf'),F/ch/'components'/(n+'.pptx'),F/ch/'components'/(n+'.eps')]:
        if str(p) in base and hashlib.sha256(p.read_bytes()).hexdigest()!=base[str(p)]:
            raise RuntimeError('Changed concurrently: '+str(p))
for name in ['10-introduction.tex','60-sns.tex','70-overall.tex']:
    path=P/'chapters'/name
    if hashlib.sha256(path.read_bytes()).hexdigest()!=base[str(path)]:raise RuntimeError('Chapter changed concurrently: '+name)
for ch,n in jobs:
    for ext in ['.pdf']:
        src=HERE/'staged'/ch/(n+ext);dst=F/ch/src.name
        if str(dst) not in base and dst.exists():
            backup=HERE/'before'/dst.relative_to(ROOT);backup.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(dst,backup)
        shutil.copy2(src,dst)
    for ext in ['.pptx','.eps','.alt.txt']:
        src=HERE/'staged'/ch/'components'/(n+ext);dst=F/ch/'components'/src.name
        dst.parent.mkdir(parents=True,exist_ok=True)
        if str(dst) not in base and dst.exists():
            backup=HERE/'before'/dst.relative_to(ROOT);backup.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(dst,backup)
        shutil.copy2(src,dst)

def replace(ch,old,new):
    p=P/'chapters'/ch
    data=p.read_bytes();a=old.encode();b=new.encode()
    if a not in data:raise RuntimeError('Missing exact caption text: '+ch+' '+old[:70])
    p.write_bytes(data.replace(a,b,1))

replace('10-introduction.tex',
    'Solid blocks are built and validated in this dissertation, and the dashed block is the physical robot the tools define the path toward, which feeds biological insight back to the benchmark.',
    'Solid blocks identify the models, evidence, and software developed or used in this dissertation. Dashed connections and the dashed robot block identify future integration and its intended return of biological insight to the benchmark.')
replace('60-sns.tex',
    'the chip fill tracks the preparation\'s support level from suspended air stepping through stance to ground.',
    'the numbered rows distinguish the preparations and their objectives. The air-stepping stages are suspended, and the reported ground evaluations retain the support harness.')
replace('60-sns.tex',
    'The full library adds the synapse-sum junction, the diameter-specific BPA force blocks, and the biological muscle block. Below, the internal implementation of the non-spiking neuron block as a native Simulink export. The full library also carries the synapse-sum junction, the diameter-specific BPA force elements, and the biological muscle element described in the text.',
    'The full library adds the synapse-sum junction, the diameter-specific BPA force blocks, and the biological muscle block. Below is the internal implementation of the non-spiking neuron block as a native Simulink export.')
replace('70-overall.tex',
    '(A) Longitudinal section with loads and constraints: internal pressure acts inside-out on the elastomer bladder; the braided aramid sleeve, woven in two handed weaves, carries the tension and converts radial expansion to axial contraction; layer contact pairs act between bladder, sleeve, and caps; with the ends clamped in the Onyx or aluminum caps the pressurized actuator develops an axial force reacted at the fixed end; and the dash-dot centerline is the symmetry axis of the half-model.',
    '(A) Idealized FEA wall layers, from outside inward: fabric covering, outer rubber, an aramid helix, middle rubber, an opposite-handed inner helix, and inner rubber. Contact interfaces couple the layers, and internal pressure acts outward on the inner wall. The cap end is fixed. Longitudinal symmetry permits expansion within the symmetry plane; at the transverse half-length boundary, axial displacement can be fixed, or a movable boundary with an optional external load can replace transverse symmetry.')

# Keep editable sources and the export entry point reproducible.
factory=(HERE/'build_figures.py').read_text(encoding='utf-8')
factory=factory.replace("OUT=ROOT/'tmp/schematic_polish_20261022/staged'","OUT=FIG")
factory=factory.replace("(OUT/'jobs.json').write_text(json.dumps(JOBS))","(OUT/'components').mkdir(exist_ok=True)\n    (OUT/'components/schematic_jobs.json').write_text(json.dumps(JOBS))")
(F/'make_pptx_figures.py').write_text(factory,encoding='utf-8')

BUILD=Path(r'D:\temp')/('diss_schematic_'+time.strftime('%Y%m%d_%H%M%S'))
(BUILD/'chapters').mkdir(parents=True)
(HERE/'build_path.txt').write_text(str(BUILD))
exe=Path(r'D:\MiKTeX\miktex\bin\x64')
env=os.environ.copy();env['BIBINPUTS']=str(P)+';';env['BSTINPUTS']=str(P)+';'
cmds=[([str(exe/'pdflatex.exe'),'-interaction=nonstopmode','-halt-on-error',f'-output-directory={BUILD}','main.tex'],P),([str(exe/'bibtex.exe'),'main'],BUILD)]
cmds+= [(cmds[0][0],P),(cmds[0][0],P)]
for i,(cmd,cwd) in enumerate(cmds):
    with (BUILD/f'pass_{i+1}.txt').open('w',encoding='utf-8') as log:
        result=subprocess.run(cmd,cwd=cwd,env=env,stdout=log,stderr=subprocess.STDOUT)
    print('pass',i+1,'exit',result.returncode,flush=True)
    if result.returncode:raise RuntimeError('Build failed; see '+str(BUILD/f'pass_{i+1}.txt'))
print('BUILD',BUILD,flush=True)
