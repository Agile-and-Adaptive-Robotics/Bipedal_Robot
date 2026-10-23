from pathlib import Path
import re
C=Path(r'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022\Revision_C');B=Path((C/'review/latest_build.txt').read_text(encoding='utf-8'))
def page(file,label):
 s=(B/'chapters'/file).read_text(encoding='utf-8',errors='replace');return int(re.search(r'\\newlabel\{'+label+r'\}\{\{[^}]*\}\{(\d+)\}',s).group(1))
d=page('50-database.aux','ch:database');n=page('60-sns.aux','ch:sns');o=page('70-overall.aux','ch:overall')
for name in ['REPORT.md','README.md']:
 p=C/name;s=p.read_text(encoding='utf-8')
 for wrong,right in [('â€”','—'),('â€“','–'),('â†’','→')]:s=s.replace(wrong,right)
 if name=='REPORT.md':
  s=s.replace('# Revision C — final local sweep, 2026-10-23','# Revision C — uploaded and advisor review updated, 2026-10-23')
  s=re.sub(r'SADb occupies printed pages \*\*.*?\*\*; SNS \*\*.*?\*\*\.',f'SADb occupies printed pages **{d}–{n-1} ({n-d} pages)**; SNS **{n}–{o-1} ({o-n} pages)** in the local verification build.',s)
  s=s.replace('The new compact Chapter 4 tables and paired figure have been visually checked.','The restored Chapter 4 tables and paired figure have been visually checked.')
 else:s='> Uploaded to the existing Overleaf project on 2026-10-23. Final UI: Errors 0, Warnings 0; 40 straightforward comments resolved, 34 substantive threads open with replies, and no pending tracked edits. See review/overleaf_upload_20261023/REPORT.md.\n\n'+s
 p.write_text(s,encoding='utf-8')
print({'local_SADb':(d,n-1),'local_SNS':(n,o-1)})
