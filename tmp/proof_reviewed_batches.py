from pathlib import Path
import subprocess,re,json,shutil
root=Path(r'D:\Github\Bipedal_Robot')
base=root/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023'
out=base/'REVIEW_INBOX/REVIEWED_20261026'
build=out/'proof_build_01'; build.mkdir(exist_ok=True)
names=['P01_replacement.tex','P01_companion_table.tex','P07_conversion_summary.tex','P06_appendix_conversion_section.tex','EG_controlled_tests.tex','EB_feedback_ablation.tex','DENG_original_engine.tex','EG_pulley_discussion.tex','EG_extensor_merged_figure.tex','EG_physical_figures_captionless_OPTION.tex']
text='\n\\clearpage\n'.join((out/'patches'/n).read_text(encoding='utf-8') for n in names)
labels=set(re.findall(r'\\label\{([^}]+)\}',text))
external=set(re.findall(r'\\ref\{([^}]+)\}',text))-labels
preamble=r'''\documentclass[11pt]{report}
\usepackage[margin=1in]{geometry}
\usepackage{graphicx,siunitx,booktabs,tabularx,array,amsmath,mathabx,natbib,url}
\graphicspath{{'''+(out/'figures').as_posix()+r'''/}}
\begin{document}
\chapter*{Targeted revision fragments: review proof}
This proof checks proposed fragments and figure placement. It is not the assembled dissertation. Cross-references to material outside the packet are represented by the placeholder X.
'''
for key in sorted(external): preamble+='\\makeatletter\\newlabel{'+key+'}{{X}{X}}\\makeatother\n'
mock=preamble[preamble.index('\\makeatletter'):]
preamble=preamble[:preamble.index('\\makeatletter')]
preamble=preamble.replace('\\begin{document}',mock+'\\begin{document}')
tex=preamble+text+'\n\\clearpage\n\\bibliographystyle{plainnat}\n\\bibliography{proof_refs}\n\\end{document}\n'
(build/'proof.tex').write_text(tex,encoding='utf-8')
bib=(base/'Revision_E/thesis.bib').read_bytes()
if b'{deng_biarticular_2023,' not in bib:
    bib+=b'\n'+(out/'deng_biarticular_2023.bib').read_bytes()
(build/'proof_refs.bib').write_bytes(bib)
bin=Path(r'D:\MiKTeX\miktex\bin\x64')
commands=[[str(bin/'pdflatex.exe'),'-interaction=nonstopmode','-halt-on-error','proof.tex'],[str(bin/'bibtex.exe'),'proof'],[str(bin/'pdflatex.exe'),'-interaction=nonstopmode','-halt-on-error','proof.tex'],[str(bin/'pdflatex.exe'),'-interaction=nonstopmode','-halt-on-error','proof.tex']]
codes=[]
for i,cmd in enumerate(commands):
    p=subprocess.run(cmd,cwd=build,capture_output=True,timeout=100)
    (build/f'command_{i}.log').write_bytes(p.stdout+p.stderr)
    codes.append(p.returncode)
    if p.returncode: break
log=(build/'proof.log').read_text(errors='replace') if (build/'proof.log').exists() else ''
result={'returncodes':codes,'errors':re.findall(r'^!.*',log,re.M),'overfull':re.findall(r'Overfull[^\n]*',log),'undefined':re.findall(r'[^\n]*undefined[^\n]*',log),'external_reference_placeholders':sorted(external),'scope':'Standalone fragment proof; not full dissertation build; existing refs outside packet mocked X.'}
(out/'proof_result.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
if len(codes)==4 and not any(codes): shutil.copy2(build/'proof.pdf',out/'REVIEWED_FRAGMENTS.pdf')
print(json.dumps(result,indent=2))
