from pathlib import Path
import re,json,hashlib,shutil,difflib,fitz
from PIL import Image,ImageOps,ImageDraw
C=Path(r'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022\Revision_C');S=C/'review/user_sweep_20261023';B=Path((C/'review/latest_build.txt').read_text(encoding='utf-8'));before=S/'before'
log=(B/'main.log').read_text(errors='replace')
for pattern in [r'^!',r'There were undefined',r'LaTeX Warning:.*undefined',r'multiply defined',r'Float too large',r'Label\(s\) may have changed']:assert not re.search(pattern,log,re.M),pattern
def table(s,key):return next(m.group() for m in re.finditer(r'\\begin\{(table|sidewaystable)\}.*?\\end\{\1\}',s,re.S) if r'\label{'+key+'}' in m.group())
def data(s,key):return table(s,key).split(r'\midrule',1)[1].split(r'\bottomrule',1)[0]
old=(before/'chapters/30-torque.tex').read_text(encoding='utf-8');new=(C/'chapters/30-torque.tex').read_text(encoding='utf-8')
for key in ['tab:FlxPinPreGof','tab:FlxPinGof','tab:ExtPinGof']:assert data(old,key)==data(new,key),key
old3=(before/'chapters/20-muscles.tex').read_text(encoding='utf-8');new3=(C/'chapters/20-muscles.tex').read_text(encoding='utf-8')
ob=data(old3,'tab:model_compare');nb=data(new3,'tab:model_compare');old_values=[]
for m in re.finditer(r'\\textbf\{[^}]*\}\s*&\s*\d+\s*&([^\n]+)',ob):old_values+=re.findall(r'\d+(?:\.\d+)?',m.group(1).split(r'\\')[0])
new_values=[]
for line in nb.splitlines():
 if '&' in line:new_values+=re.findall(r'\d+(?:\.\d+)?',line.split('&',3)[-1])
assert len(old_values)==36 and old_values==new_values,(old_values,new_values)
changes=[]
for p in (C/'chapters').glob('*.tex'):
 prev=before/'chapters'/p.name
 if prev.exists() and p.read_bytes()!=prev.read_bytes():changes.append(''.join(difflib.unified_diff(prev.read_text(encoding='utf-8').splitlines(True),p.read_text(encoding='utf-8').splitlines(True),fromfile='before/'+p.name,tofile='chapters/'+p.name)))
(S/'changes.diff').write_text('\n'.join(changes),encoding='utf-8')
for p in (before/'figs/20-muscles').rglob('*'):
 if p.is_file():assert p.read_bytes()==(C/'figs/20-muscles'/p.relative_to(before/'figs/20-muscles')).read_bytes(),p
doc=fitz.open(B/'main.pdf');offset=next(p for l,t,p in doc.get_toc() if l==1 and t=='Introduction')-1;pages={}
for chapter,keys in [('15-background',['fig:festo_bpa_cad','fig:steele_knee','fig:aarl_opensim_comparison','fig:masters_routing_and_arms']),('20-muscles',['fig:Fmax_fit','fig:F-P_contraction','fig:GenForceFit','tab:model_compare']),('30-torque',['fig:bktFrame','fig:xiProjection','tab:FlxPinPreGof','tab:FlxPinGof','tab:ExtPinGof','tab:AdoptedXi','fig:ExtensorRoute_redesigned']),('50-database',['ch:database']),('60-sns',['tab:sns_capability'])]:
 aux=(B/'chapters'/f'{chapter}.aux').read_text(errors='replace')
 for key in keys:
  m=re.search(r'\\newlabel\{'+re.escape(key)+r'\}\{\{([^}]+)\}\{(\d+)\}',aux);assert m,key
  n=int(m.group(2))+offset-1;name=key.split(':')[1];doc[n].get_pixmap(matrix=fitz.Matrix(1.4,1.4)).save(S/(name+'_page.png'));pages[key]={'number':m.group(1),'printed':m.group(2),'pdf':n+1}
shutil.copy2(B/'main.pdf',C/'Revision_C_user_sweep_20261023.pdf')
(S/'verification.json').write_text(json.dumps({'build':str(B),'pdf_pages':len(doc),'pages':pages,'force_plot_files_byte_identical':True,'table3_3_all_36_results_unchanged':True,'tables4_1_to4_3_data_bodies_unchanged':True},indent=2),encoding='utf-8')
print(json.dumps(pages,indent=2))
