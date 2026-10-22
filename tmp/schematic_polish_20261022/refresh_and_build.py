from pathlib import Path
import json,shutil,subprocess,os,time

H=Path(__file__).resolve().parent;R=H.parent.parent
D=R/'Documentation/Reports and Papers/Dissertation';F=D/'Figures';P=D/'ProofFinal'
jobs=json.loads((H/'staged/jobs.json').read_text())
def copy_safe(src,dst):
    temp=dst.with_name(dst.name+'.replacement')
    for attempt in range(10):
        try:
            shutil.copy2(src,temp);os.replace(temp,dst);return
        except OSError:
            if attempt==9:raise
            time.sleep(.4)
native=F/'60-sns/components/simulink_neuron_detail_native_export.pdf'
if not native.exists():shutil.copy2(F/'60-sns/simulink_neuron_detail.pdf',native)
for ch,n in jobs:
    copy_safe(H/'staged'/ch/(n+'.pdf'),F/ch/(n+'.pdf'))
    for e in ['.pptx','.eps','.alt.txt']:copy_safe(H/'staged'/ch/'components'/(n+e),F/ch/'components'/(n+e))
factory=(H/'build_figures.py').read_text(encoding='utf-8').replace("OUT=ROOT/'tmp/schematic_polish_20261022/staged'","OUT=FIG")
factory=factory.replace("(OUT/'jobs.json').write_text(json.dumps(JOBS))","(OUT/'components').mkdir(exist_ok=True)\n    (OUT/'components/schematic_jobs.json').write_text(json.dumps(JOBS))")
(F/'make_pptx_figures.py').write_text(factory,encoding='utf-8')
(F/'components').mkdir(exist_ok=True)
shutil.copy2(H/'staged/jobs.json',F/'components/schematic_jobs.json')
p=P/'chapters/60-sns.tex';data=p.read_bytes()
data=data.replace(b'width=\\textwidth,height=0.21\\textheight,keepaspectratio]{simulink_neuron_detail.pdf}',b'width=\\textwidth]{simulink_neuron_detail.pdf}',1)
data=data.replace(b'which also shows the internal implementation of the non-spiking neuron block as a native Simulink export (bottom)',b'which also shows a functional schematic of the non-spiking neuron implementation (bottom)',1)
data=data.replace(b'Below is the internal implementation of the non-spiking neuron block as a native Simulink export.',b'Below is a functional schematic of the non-spiking neuron implementation: synaptic, leak, and applied currents feed membrane integration, and a threshold, inverse-slope gain, and saturation produce the bounded drive output. The voltage is fed back to both current calculations.',1)
p.write_bytes(data)
p=P/'chapters/70-overall.tex';data=p.read_bytes();data=data.replace(b'(B) Exploded half-cutaway of the layered construction,',b'(B) Illustrative exploded half-cutaway of the layered construction,',1);p.write_bytes(data)
data=data.replace(b'{bpa_fea_exploded_half_cutaway_v2.png}',b'{bpa_exploded_cutaway.pdf}',1)
start=data.index(b'\\caption[Construction of the braided pneumatic actuator]')
end=data.index(b'\\label{fig:bpa_construction}',start)
caption=r'''\caption[Construction of the braided pneumatic actuator]{Construction represented by the proposed finite-element model. (A) Six idealized wall layers with contact interfaces, outward internal pressure, and a fixed cap end. Longitudinal symmetry permits in-plane expansion; at the transverse half-length boundary, axial displacement is fixed, or a movable boundary with an optional external load replaces transverse symmetry. (B) Illustrative exploded cutaway showing the fabric covering, rubber, and opposite-handed aramid reinforcement anchored at the cap. The plane callout identifies longitudinal symmetry; the rubber callout identifies retained material above it. The middle rubber layer is explicit in (A) but not separated in (B). Layer spacing and cap details are schematic. Not to scale.}
'''
data=data[:start]+caption.encode()+data[end:];p.write_bytes(data)
B=Path(r'D:\temp')/('diss_schematic_'+time.strftime('%Y%m%d_%H%M%S'));(B/'chapters').mkdir(parents=True)
(H/'build_path.txt').write_text(str(B))
exe=Path(r'D:\MiKTeX\miktex\bin\x64');env=os.environ.copy();env['BIBINPUTS']=str(P)+';';env['BSTINPUTS']=str(P)+';'
latex=[str(exe/'pdflatex.exe'),'-interaction=nonstopmode','-halt-on-error',f'-output-directory={B}','main.tex']
for i,(cmd,cwd) in enumerate([(latex,P),([str(exe/'bibtex.exe'),'main'],B),(latex,P),(latex,P)]):
    with (B/f'pass_{i+1}.txt').open('w',encoding='utf-8') as log:result=subprocess.run(cmd,cwd=cwd,env=env,stdout=log,stderr=subprocess.STDOUT)
    print('pass',i+1,result.returncode,flush=True)
    if result.returncode:raise RuntimeError('Build failed; see '+str(B))
print('BUILD',B,flush=True)
