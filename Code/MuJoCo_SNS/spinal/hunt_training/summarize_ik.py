"""Readable 3D IK/PF evidence with explicit source and constraint limits."""
import os,sys,json,xml.etree.ElementTree as ET
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
H=Path(__file__).resolve().parent;S=H.parent
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from motion_io import read_motion
R=H/'runs/ik_audit';audit=json.loads((R/'audit.json').read_text())
source=S.parents[2]/'Solid_Models/OpenSim/Gait2392_Robotbody'
tasks={n.attrib['name']:{e.tag:e.text for e in n} for n in ET.parse(source/'subject01_Setup_IK.xml').iter('IKCoordinateTask')}
requested=['hip_adduction_r','hip_rotation_r','lumbar_bending','lumbar_rotation','lumbar_extension','subtalar_angle_r','mtp_angle_r']
labels=['Hip x: abduction/adduction','Hip y: axial rotation','Lumbar x: lateral bending','Lumbar y: axial rotation','Lumbar z: flexion/extension','Subtalar: oblique inversion/eversion','MTP: toe flexion/extension']
plt.rcParams.update({'font.family':'Arial','font.size':10})
fig,axes=plt.subplots(4,2,figsize=(12,12),constrained_layout=True)
for record,color,style in zip(audit['full_motion_sources'],['#0077BB','#CC6677','#228833'],['-','--',':']):
    t,names,v,unit=read_motion(Path(record['path']))
    if 'subject01' in record['name']:start,end=.617,1.85
    else:start,end=float(t[0]),float(t[-1])
    phase=np.linspace(0.,1.,201,endpoint=False)
    for joint,title,ax in zip(requested,labels,axes.flat):
        if joint not in names:continue
        y=np.interp(start+phase*(end-start),t,v[:,names.index(joint)])
        if unit=='no':y=np.rad2deg(y)
        ax.plot(phase*100,y,color=color,ls=style,lw=1.6,label=record['name'])
        ax.set(title=title,xlabel='Source cycle (%)',ylabel='Source coordinate (degrees)')
axes.flat[-1].axis('off');handles,text=axes.flat[0].get_legend_handles_labels();axes.flat[-1].legend(handles,text,loc='center',fontsize=9)
fig.suptitle('Full 3D references: motion evidence, not identified neural PF populations',fontsize=13)
fig.savefig(H/'IK_COORDINATE_PATTERNS.png',dpi=160);plt.close(fig)
(H/'IK_COORDINATE_PATTERNS.alt.txt').write_text('Seven plots compare hip frontal/axial, three lumbar, subtalar, and MTP source coordinates. Experimental IK contains hip and lumbar variation but its foot coordinates are targeted to zero. Falisse prediction includes subtalar and MTP variation. The tutorial normal lumbar curves exactly counter pelvis rotations. Sources are separate examples, not a multi-subject ensemble.')
lines=['# Full 3D IK and PF-pattern audit','',
 'Machine: easteregg2. This is training evidence for the complete Gait2392 model, not a sagittal restriction.','',
 '## Scope and source limitations','',
 'The 43 cached references retain only sagittal hip/knee/ankle curves. All ten Ong references are planar predictions from a nine-DoF model with lumbar motion locked, so they cannot support independent out-of-plane/lumbar PF patterns. The original 32 running-trial full files are not present on this machine; their discarded coordinates cannot be reconstructed from caches. No multi-subject full-3D conclusion is claimed.','',
 'Three complete motion files were analyzed: the experimental subject01 IK, the illustrative Gait2392 normal example, and the recovered Falisse Case_40 prediction. The latter is not experimental IK. Its source commit and file hashes are in ignored runs/ik_sources/falisse_manifest.json.','',
 '| Coordinate | Experimental IK RoM | Normal example RoM | Falisse prediction RoM | Native muscles with >1 mm leverage |','|---|---:|---:|---:|---:|']
for joint,label in zip(requested,labels):
    values=[]
    for r in audit['full_motion_sources']:
        c=r['coordinates'][joint];values.append('missing' if c['status']=='missing' else f"{c['range_deg']:.3f} deg")
    lines.append('| '+label+' | '+' | '.join(values)+f" | {audit['muscle_authority'][joint]['muscles_over_1mm']} |")
lines+=['','## What this supports','',
 '- Retain distinct controllable hip frontal/axial and lumbar x/y/z channels. They can share the RG clock while having separate PF distribution, sensory gains, and motor recruitment. Distinct movement is not proof of an independent oscillator.','- The tutorial normal example sets all three lumbar coordinates exactly opposite the corresponding pelvis rotations (sum standard deviation zero). That is imposed trunk stabilization, not independent evidence of three new rhythm generators. The experimental lumbar curves are not that exact identity.','- Both subtalar and MTP coordinates in the experimental IK are explicitly targeted to their zero defaults by applied IKCoordinateTasks with weight 1000. The model coordinates are unlocked. Flat curves therefore reflect the IK setup, not absent human motion or absent neural control.','- The native Gait2392 model has four muscles per MTP and twelve per subtalar joint with nontrivial leverage over the sampled poses. These joints are mechanically controllable. Their periodic targets are not identified by the experimental zero-default solution.','- Falisse predicts approximately 37.5 degrees of MTP and 16.9 degrees of subtalar excursion. However, that model changed the toe axis and its toe actuation is passive. Its toe motion does not prove an active MTP PF pattern. The native model instead retains an oblique MTP axis and muscle authority.','- The native subtalar axis has a substantial anatomical-y component as well as x. Control its actual joint coordinate; do not manufacture a pure x hinge. The native MTP axis is also oblique, unlike the modified Falisse axis.','',
 'A sample-holdout sagittal regression was included as a descriptive dependency check only. It is not subject-level cross-validation or a neural model-selection test. Independent-subject activation/force data and perturbations remain necessary before claiming new biological primitives.','',
 '## Sign and training corrections','',
 'Actual source/OpenSim and converted/MuJoCo forward-kinematics perturbations require -1 for bilateral hip flexion/adduction and +1 for the remaining leg/lumbar coordinates. This replaces the old training data sign guesses only in the new campaign. All 17 joint coordinates remain free. Full root/frame and target retargeting assumptions are recorded in runs/targets/metadata.json.','',
 'The initial inverse-force target reconstruction uses native equality-aware moment arms, bounded native muscle activation, and force-velocity effects. It retains all 17 coordinates and uses an explicitly supported-in-place load reference. It is an engineering training target, not a new measurement of human muscle activation.','',
 '![Full coordinate patterns](IK_COORDINATE_PATTERNS.png)','',
 'Primary sources: [Falisse et al. 2022](https://doi.org/10.1371/journal.pone.0256311), [author model/results](https://github.com/antoinefalisse/predictsim_mtp), and [Ong et al. 2019 model scope](https://doi.org/10.1371/journal.pcbi.1006993). Local source: Solid_Models/OpenSim/Gait2392_Robotbody/subject01_Setup_IK.xml.']
(H/'IK_PF_AUDIT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print('Wrote IK_PF_AUDIT.md and IK_COORDINATE_PATTERNS.png')
