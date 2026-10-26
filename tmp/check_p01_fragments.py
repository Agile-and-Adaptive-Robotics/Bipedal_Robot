from pathlib import Path
import subprocess,json,re
E=Path(r'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261023\Revision_E');V=E/'review/targeted_reconciliation';W=V/'tools_inventory/tex_check';W.mkdir(exist_ok=True)
src=r'''\documentclass[12pt]{article}
\usepackage[margin=1in]{geometry}
\usepackage{tabularx,booktabs,array,siunitx}
\renewcommand{\thetable}{6.\arabic{table}}
\makeatletter\newlabel{ch:muscles}{{3}{1}}\makeatother
\begin{document}
\section*{Proposed P01: tools and demonstrations}
'''+(V/'P01_replacement.tex').read_text(encoding='utf-8')+'\n'+(V/'P01_companion_table.tex').read_text(encoding='utf-8')+'\n\\end{document}\n'
(W/'P01_check.tex').write_text(src,encoding='utf-8')
cmd=[r'D:\MiKTeX\miktex\bin\x64\pdflatex.exe','-interaction=nonstopmode','-halt-on-error','P01_check.tex']
for n in (1,2):
 with (W/f'pass{n}.log').open('w',encoding='utf-8') as f:r=subprocess.run(cmd,cwd=W,stdout=f,stderr=subprocess.STDOUT)
 if r.returncode:raise RuntimeError(str(W))
log=(W/'P01_check.log').read_text(errors='replace');j={'scope':'Standalone fragment compile, not a replacement dissertation build','errors':len(re.findall(r'^!',log,re.M)),'undefined_references':len(re.findall(r'undefined references|Reference .* undefined',log)),'float_too_large':bool(re.search('Float too large',log)),'overfull_boxes':len(re.findall(r'Overfull \\hbox',log))}
(W/'result.json').write_text(json.dumps(j,indent=2),encoding='utf-8');print(json.dumps(j))
