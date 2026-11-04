"""Read-only comparison of saved closure and geometric strain fields."""
from pathlib import Path
import csv,json,numpy as np
here=Path(__file__).resolve().parents[1]/'packaging'
rs=json.loads((here/'build/saved_routes.json').read_text())
rows=[]
for r in rs:
 bounds=np.asarray(r['spec']['kin']['thetaRangesDeg']);s=np.asarray(r['sContraction']);strain=np.asarray(r['saved_strain'])
 i=np.abs(np.linspace(*bounds[0],s.shape[0])).argmin();j=np.abs(np.linspace(*bounds[1],s.shape[1])).argmin()
 rows.append(dict(site=r['site'],saved_sContraction_over_rest=s[i,j]/r['RestingL'],saved_geometric_strain=strain[i,j],negative_geometric_strain_cells=int((strain<0).sum()),cells=strain.size,note='Different saved fields; not a new equilibrium solve or a corrected result.'))
with (here/'output/saved_state_field_audit.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print(rows[0])
