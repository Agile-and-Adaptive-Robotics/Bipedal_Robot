from pathlib import Path
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np,json,hashlib
W=Path(__file__).resolve().parent;R=next(p for p in W.parents if p.name=='Bipedal_Robot')
A=W/'assets/identification';A.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':17,'axes.spines.top':False,'axes.spines.right':False,'mathtext.fontset':'stix'})
report=[]
for kind,fn,row,hold,tc,vc in [('flexor','minimizeFlxPin10_results.mat',76,[1,2],[6,7,8],[9,10,11]),('extensor','minimizeExt10mmX3_results.mat',31,[1,2,3],[8,9,10],[11,12,13])]:
    source=R/'Testing_Data/2022_02_Festo'/fn
    m=loadmat(source,simplify_cells=True,variable_names=['filtered_results','results_cv','population','lb','ub'])
    fr=m['filtered_results'];chosen=fr[row]
    fold=np.all(fr[:,hold]==chosen[hold],axis=1)
    points=fr[fold]
    fig,axes=plt.subplots(1,2,figsize=(11,4.1),layout='constrained')
    for ax,cols,name in zip(axes,[tc,vc],['Fitting tests','Validation tests']):
        ax.scatter(points[:,cols[0]],points[:,cols[2]],s=65,c='#4477AA',alpha=.8,label='Candidate')
        ax.scatter(chosen[cols[0]],chosen[cols[2]],s=180,marker='*',c='#AA3377',edgecolors='black',linewidths=.5,label='Adopted')
        ax.set(xlabel='Normalized RMSE',ylabel='Normalized maximum error',title=name)
        ax.grid(alpha=.16);ax.set_axisbelow(True)
    axes[1].legend(frameon=False,fontsize=14)
    fig.savefig(A/f'{kind}_pareto_projection.png',dpi=210,bbox_inches='tight');plt.close(fig)
    report.append({'case':kind,'source':str(source),'sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'adopted_row_1based':row+1,'held_out_test_indices':chosen[hold].tolist(),'same_fold_candidates':len(points),'plotted_columns_1based':[[x+1 for x in tc],[x+1 for x in vc]],'note':'Two-dimensional projections of the retained three-objective candidate set from the adopted fold. Validation was used in candidate filtering; no blind-validation claim.'})
equations={
 'route_cost':r'$J=10^5J_{\mathrm{worst}}+10^3J_{\mathrm{torque}}+\sum_a J_a+J_{\mathrm{geom}}+J_{\mathrm{len}}$',
 'shortfall':r'$s_j=\frac{\max(0,\,1.05|M^H_j|-|M^R_j|)}{\tau},\quad \tau=\max(1\,\mathrm{N\,m},\max_j|M^H_j|)$',
 'worst_mean':r'$J_{\mathrm{worst}}=\max_j s_j,\qquad J_{\mathrm{torque}}=\frac{1}{N}\sum_j s_j^2$',
 'ident_objective':r'$g_k(\boldsymbol{\chi})=\frac{1}{|\mathcal{T}|}\sum_{i\in\mathcal{T}}\frac{g_{k,i}(\boldsymbol{\chi})}{g^{(0)}_{k,i}},\quad k\in\{\mathrm{RMSE},\mathrm{FVU},\mathrm{Max}\}$',
}
for name,text in equations.items():
    fig=plt.figure(figsize=(12,1.3),dpi=240);fig.text(.5,.5,text,ha='center',va='center',fontsize=30,color='#202A25')
    fig.savefig(A/f'{name}.png',transparent=True,bbox_inches='tight',pad_inches=.04);plt.close(fig)
(W/'build/pareto_provenance.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
