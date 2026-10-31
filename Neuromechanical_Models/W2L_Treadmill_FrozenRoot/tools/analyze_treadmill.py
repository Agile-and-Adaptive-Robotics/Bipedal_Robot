import numpy as np,json
from pathlib import Path
p=Path('Neuromechanical_Models/W2L_Treadmill_FrozenRoot')
a=np.loadtxt(p/'Treadmill verification.txt',skiprows=1);a=a[(a[:,4]>.9)&(a[:,1]<=10)]
names=(p/'Treadmill verification.txt').read_text().splitlines()[0].split('\t')
out={names[i]:{'initial':float(a[0,i]),'final':float(a[-1,i]),'min':float(a[:,i].min()),'max':float(a[:,i].max())} for i in range(2,a.shape[1])}
for i,k in enumerate(names):
 if 'contacts' in k:
  hits=a[a[:,i]>0,1];out[k]['sample_count']=int(len(hits));out[k]['last_contact_time']=float(hits[-1]) if len(hits) else None
(p/'TEST_RESULTS.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k in ['belt_x','belt_y','root_y','slide_pos'] or 'contacts' in k},indent=2))

