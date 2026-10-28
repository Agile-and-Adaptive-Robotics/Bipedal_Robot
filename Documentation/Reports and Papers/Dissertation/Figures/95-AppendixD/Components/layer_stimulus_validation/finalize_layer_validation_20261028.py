from pathlib import Path
import json,csv,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
D=Path('D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation');A=D/'Archive/flow_and_validation_20261028';F=D/'Figures/95-AppendixD';P=D/'ProofFinal/chapters'
rows=json.loads((A/'stimulus_long_controls/stimulus_results.json').read_text());base=np.load(A/'stimulus_long_controls/baseline.npz')
plt.rcParams.update({'font.family':'Arial','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
fig,axs=plt.subplots(3,3,figsize=(8.4,7.3),sharex=True);checks=[]
for i,(name,label) in enumerate([('RG_plus50_2cycles','RG: +50 nA'),('PF_plus50_2cycles','PF: +50 nA'),('MNflex_minus500_2cycles','MN: −500 nA')]):
 row=next(r for r in rows if r['condition']==name);z=np.load(A/'stimulus_long_controls'/f'{name}.npz');t=z['t']-row['start_s'];mask=(t>=-2)&(t<=9);pre=t<0
 assert np.array_equal(z['V'][pre],base['V'][pre]);assert np.array_equal(z['q'][pre],base['q'][pre])
 def bursts(y):
  inds=np.flatnonzero((y[1:]>.5)&(y[:-1]<=.5))+1;return int(np.sum((z['t'][inds]>=row['start_s'])&(z['t'][inds]<row['start_s']+row['duration_s'])))
 check={'condition':name,'baseline_burst_onsets_during_stimulus':bursts(base['motor'][:,1]),'stimulated_burst_onsets_during_stimulus':bursts(z['motor'][:,1]),'pre_stimulus_voltage_and_mechanics_bit_exact':True};checks.append(check)
 for j,(key,col,scale,ylabel) in enumerate([('V',0,1,'RG extensor, mV'),('motor',1,1,'Flexor MN activation'),('q',0,180/np.pi,'Hip angle, deg')]):
  ax=axs[i,j];ax.plot(t[mask],base[key][mask,col]*scale,'--',color='.35',lw=1.2,label='Baseline');ax.plot(t[mask],z[key][mask,col]*scale,color='#0072B2',lw=1.1,label='Stimulus');ax.axvspan(0,row['duration_s'],color='#EECC66',alpha=.35);ax.set_ylabel(ylabel);ax.grid(alpha=.15)
 axs[i,0].set_title(label,loc='left',fontweight='bold')
for ax in axs[-1]:ax.set_xlabel('Time from stimulus onset, s')
fig.legend(*axs[0,0].get_legend_handles_labels(),loc='upper center',ncol=2,frameon=False);fig.tight_layout(rect=(0,0,1,.955));fig.savefig(F/'layer_stimulus_validation.pdf',bbox_inches='tight');fig.savefig(A/'layer_stimulus_validation.png',dpi=170);plt.close(fig)
(A/'whole_burst_checks.json').write_text(json.dumps(checks,indent=2));print(checks)
comp=F/'Components/layer_stimulus_validation';shutil.copy2(__file__,comp/Path(__file__).name)
allrows=[]
for series in ['stimulus_layers','stimulus_controls','stimulus_long_controls']:
 for row in json.loads((A/series/'stimulus_results.json').read_text()):allrows.append({'series':series,**row})
with (comp/'stimulus_results_all.csv').open('w',newline='') as f:
 fields=list(dict.fromkeys(k for r in allrows for k in r));w=csv.DictWriter(f,fields);w.writeheader();w.writerows(allrows)
for n in ['60-sns','95-AppendixD']:
 p=P/(n+'.tex');s=p.read_text(encoding='utf8')
 if n=='60-sns':
  s=s.replace('A one-second RG stimulus shifted both RG and motor timing by 0.307 cycles. A shorter PF stimulus and direct MN inhibition suppressed flexor output without a lasting phase shift. A longer PF stimulus shifted motor timing by 0.208 cycles while leaving the RG clock unchanged.','Stimuli lasting two baseline periods suppressed two expected flexor motor bursts at each layer. RG stimulation shifted both RG and motor timing by 0.211 cycles, whereas PF stimulation and direct MN inhibition preserved the original schedule. A separate one-second PF stimulus shifted motor timing by 0.208 cycles while leaving the RG clock unchanged.')
 else:
  s=s.replace('A subsequent stronger-current control series retained the same circuit and mechanical parameters','Subsequent stronger-current control series, including stimuli lasting two baseline periods (\\qty{3.876}{\\second}), retained the same circuit and mechanical parameters')
  s=s.replace('In particular, the longer PF pulse shifted the motor pattern without resetting the RG clock.','The two-cycle stimuli suppressed both expected flexor burst onsets at all three layers. Baseline and stimulated neural and mechanical trajectories were identical before stimulus onset. The one-second PF pulse shifted the motor pattern without resetting the RG clock, whereas the two-cycle PF pulse preserved both schedules. The duration dependence cautions against treating every motor phase shift as an RG reset.')
  s=s.replace('PF & $+50$ & 0.5 & 0.000 & 0.000 & Non-resetting deletion','PF & $+50$ & 0.5 & 0.000 & 0.000 & Brief suppression; schedule preserved')
  s=s.replace('MN & $-500$ & 0.5 & 0.000 & 0.000 & Non-resetting deletion','MN & $-500$ & 0.5 & 0.000 & 0.000 & Brief suppression; schedule preserved')
  s=s.replace('\\bottomrule\n\\end{tabularx}\n\\end{table}\n\n\\begin{figure}[p]\n\\centering\n\\includegraphics[width=\\textwidth,height=0.74',r'''RG & $+50$ & 3.876 & 0.211 & 0.211 & Resetting deletion \\
PF & $+50$ & 3.876 & 0.000 & 0.000 & Non-resetting deletion \\
MN & $-500$ & 3.876 & 0.000 & 0.000 & Non-resetting deletion \\
\bottomrule
\end{tabularx}
\end{table}

\begin{figure}[p]
\centering
\includegraphics[width=\textwidth,height=0.74''')
  s=s.replace('Rows show RG stimulation, PF stimulation, and direct flexor-MN inhibition; columns','All stimuli last two baseline periods. Rows show RG stimulation, PF stimulation, and direct flexor-MN inhibition; columns')
 p.write_text(s,encoding='utf8')
