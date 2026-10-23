from pathlib import Path
import re,json
R=Path(r'D:\Github\Bipedal_Robot');C=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C';p=C/'review/final_sweep_20261023/app_payload_check.json'
raw=(R/'SADb_audit/app/sadb_app.html').read_text(encoding='utf-8');d=json.loads(re.search(r'<script type="application/json" id="sadb-data">(.*?)</script>',raw,re.S).group(1));papers=[r for r in d if r['rt']=='paper' and not r['ar']]
counts={k:sum(bool(r[k]) for r in papers) for k in ['n','an','ast','fb','af','rs','d']}
j=json.loads(p.read_text());j['retained_paper_field_coverage']=counts;p.write_text(json.dumps(j,indent=2));print(counts)
