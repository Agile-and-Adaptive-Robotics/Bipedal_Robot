"""Check source conservation, final build, page locations, and directive evidence."""
from pathlib import Path
import json,re,hashlib,xml.etree.ElementTree as ET
R=next(p for p in Path(__file__).resolve().parents if (p/'AGENTS.md').exists());D=R/'Documentation/Reports and Papers/Dissertation';N=Path(__file__).parent;A=D/'Archive/session1_sns_20261031';B=D/'Archive/build_session1_sns_20261031_v5'
read=lambda p:p.read_bytes().decode('utf-8').replace('\r\n','\n')
checks=[]
for name in ('60-sns.tex','95-AppendixD.tex'):
 p=D/'ProofFinal/chapters'/name;old=read(A/'ProofFinal/chapters'/name);new=read(p)
 labels=lambda s:set(re.findall(r'\\label\{([^}]+)\}',s))
 checks.append({'file':str(p),'removed_labels':sorted(labels(old)-labels(new)),'added_labels':sorted(labels(new)-labels(old)),'final_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'numeric_tokens_no_longer_present':sorted(set(re.findall(r'\d+(?:\.\d+)?',old))-set(re.findall(r'\d+(?:\.\d+)?',new)))})
 assert set(checks[-1]['removed_labels']) <= {'tab:spike_knee','fig:spike_knee'}
old=(A/'ProofFinal/thesis.bib').read_bytes();new=(D/'ProofFinal/thesis.bib').read_bytes();assert new.startswith(old)
log=read(B/'main.log');errors=re.findall(r'(?m)^!.*',log);warnings=re.findall(r'(?m)^LaTeX Warning:.*',log)
assert not errors and not warnings,(errors,warnings)
bbox=read(N/'REVISED_PDF_BBOX_20261031.xhtml');bbox=re.sub('[\x00-\x08\x0b\x0c\x0e-\x1f]','',bbox)
root=ET.fromstring(bbox);pages=[' '.join(w.text or '' for w in p.iter() if w.tag.endswith('}word')) for p in root.iter() if p.tag.endswith('}page')];locations={}
for key in ('6.2','6.4','6.5','D1','D6'):
 locations[key]=[i+1 for i,p in enumerate(pages) if re.search(r'Figure\s+'+re.escape(key)+r'\s*:',p)]
report={'machine':'easteregg2','build':str(B/'main.pdf'),'latex_errors':errors,'latex_warnings':warnings,'undefined_reference_or_citation':bool(re.search(r'(Reference|Citation).*undefined|undefined references',log)),'bibliography_old_bytes_preserved':True,'source_checks':checks,'figure_pdf_pages_1_based':locations,'page_count':len(pages)}
(N/'FINAL_CHECKS_20261031.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
