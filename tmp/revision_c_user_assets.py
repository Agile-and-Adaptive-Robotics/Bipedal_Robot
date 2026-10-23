from pathlib import Path
import fitz, numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Arc
C=Path(r'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022\Revision_C');F=C/'figs';S=C/'review/user_sweep_20261023'
# Native PDF relayout: preserve source image pixels and change only panel lettering.
d=fitz.open(F/'15-background/steeleknee.pdf');p=d[0]
for x,t in [(0,'(A)'),(200,'(B)')]:
 p.draw_rect(fitz.Rect(x,0,x+32,28),color=None,fill=(1,1,1));p.insert_text((x+3,20),t,fontsize=18,fontname='hebo')
d.save(S/'steeleknee_uppercase.pdf');d.close();(F/'15-background/steeleknee.pdf').write_bytes((S/'steeleknee_uppercase.pdf').read_bytes())
d=fitz.open(F/'30-torque/bktFrame.pdf');p=d[0];spans=[sp for b in p.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for sp in l['spans'] if sp['text'] in ['(A)','(B)','(C)','(D)']]
assert len(spans)==4
for sp in spans:p.add_redact_annot(fitz.Rect(sp['bbox']),fill=(1,1,1))
p.apply_redactions(images=0)
for sp in spans:p.insert_text(sp['origin'],sp['text'],fontname='hebo',fontsize=sp['size'])
d.save(S/'bktFrame_bold.pdf');d.close();(F/'30-torque/bktFrame.pdf').write_bytes((S/'bktFrame_bold.pdf').read_bytes())
# Split the existing pair into robot (right) and human (left) native PDF panels.
d=fitz.open(F/'15-background/masters_leg_model.pdf');r=d[0].rect
for name,clip in [('human_gait2392_paths.pdf',fitz.Rect(0,0,r.width/2,r.height)),('robot_pam_routes.pdf',fitz.Rect(r.width/2,0,r.width,r.height))]:
 out=fitz.open();p=out.new_page(width=clip.width,height=clip.height);p.show_pdf_page(p.rect,d,0,clip=clip);out.save(F/'15-background'/name)
plt.rcParams.update({'font.family':'Arial','font.size':11,'mathtext.fontset':'dejavusans','pdf.fonttype':42})
blue='#0072B2';orange='#D55E00';teal='#009E73';ink='#263442'
fig,axs=plt.subplots(1,2,figsize=(8.4,3.9))
u=np.array([.5,np.sqrt(3)/2]);q=np.array([2.,0]);h=q-q.dot(u)*u;n=np.array([u[1],-u[0]])
for i,ax in enumerate(axs):
 ax.set_aspect('equal');ax.set_xlim(-.55,3.6);ax.set_ylim(-1.85,3.15);ax.axis('off')
 ax.text(-.5,2.92,'(A)' if i==0 else '(B)',weight='bold',fontsize=13)
 ax.text(.15,2.92,'Force moment' if i==0 else 'Moment-arm projection',weight='bold')
 ax.plot([q[0]-1.6*u[0],q[0]+2.5*u[0]],[q[1]-1.6*u[1],q[1]+2.5*u[1]],color=blue,lw=1.7)
 ax.scatter([0,q[0]],[0,q[1]],s=21,c=ink)
 ax.text(-.32,-.1,'O');ax.text(2.02,-.28,'Q')
 ax.annotate('',q+1.6*u,q,arrowprops={'arrowstyle':'-|>','lw':2,'color':blue});ax.text(2.9,1.35,r'$\mathbf{F}=F\hat{\mathbf{u}}$',color=blue,fontsize=10)
 ax.annotate('',q,(0,0),arrowprops={'arrowstyle':'->','lw':1.6,'color':ink});ax.text(.7,.13,r'$\mathbf{d}$' if i==0 else r'$\mathbf{p}_{\mathrm{proj},xy}$',fontsize=11)
 ax.plot([0,h[0]],[0,h[1]],color=orange,lw=2);ax.text(.55,-.8,r'$r_k$',color=orange)
 ax.plot([h[0],q[0]],[h[1],q[1]],color=blue,lw=1.1,ls=':')
 corner=h-.12*n;ax.plot([corner[0],corner[0]+.12*u[0],h[0]+.12*u[0]],[corner[1],corner[1]+.12*u[1],h[1]+.12*u[1]],color=ink,lw=.8)
 ax.add_patch(Circle((-.27,.68),.1,fill=False,color=ink));ax.scatter([-.27],[.68],s=6,c=ink);ax.text(-.08,.62,r'$\hat{\mathbf{z}}_k$')
 if i==0:ax.text(-.35,-1.68,r'$M_z=\hat{\mathbf{z}}_k\cdot(\mathbf{d}\times\mathbf{F})=r_k F$',fontsize=10)
 else:
  ax.annotate('',1.2*n,(0,0),arrowprops={'arrowstyle':'-|>','lw':1.2,'color':teal});ax.text(.95,-.37,r'$\hat{\mathbf{n}}$',color=teal)
  ax.text(-.35,-1.68,r'$r_k=\mathbf{p}_{\mathrm{proj},xy}\cdot\hat{\mathbf{n}},\quad\hat{\mathbf{n}}=\frac{\mathbf{p}_f\times\hat{\mathbf{z}}_k}{\|\mathbf{p}_f\times\hat{\mathbf{z}}_k\|}$',fontsize=9)
fig.subplots_adjust(left=.03,right=.98,top=.98,bottom=.03,wspace=.20);fig.savefig(F/'15-background/moment_arm_methods.pdf');plt.close(fig)
fig,ax=plt.subplots(figsize=(7,2.1));ax.set_aspect('equal');ax.axis('off');ax.set_xlim(-2.6,2.6);ax.set_ylim(-.25,1.55)
ax.add_patch(Circle((0,0),.91,color='#E9EEF2',ec='#80909A',lw=1.1));angles=np.deg2rad(np.linspace(30,150,120));ax.plot(np.cos(angles),np.sin(angles),color=blue,lw=9,solid_capstyle='round')
for a in [30,150]:
 th=np.deg2rad(a);pt=np.array([np.cos(th),np.sin(th)]);direction=np.array([-np.sin(th),np.cos(th)])*(-1 if a==150 else 1);end=pt-direction*.8;ax.plot([pt[0],end[0]],[pt[1],end[1]],color=blue,lw=9,solid_capstyle='round');ax.plot([0,pt[0]],[0,pt[1]],ls='--',color=ink,lw=.8)
ax.annotate('',(0,1),(0,0),arrowprops={'arrowstyle':'<->','lw':1,'color':orange});ax.text(.08,.65,r'$R$',color=orange)
ax.add_patch(Arc((0,0),.75,.75,theta1=30,theta2=150,color=orange,lw=1.2));ax.text(-.75,1.25,'BPA following joint surface',fontsize=10)
ax.text(-2.55,.70,r'$\alpha=|\theta_{\mathrm{wrap}}-\theta_k|$',fontsize=11);ax.text(1.1,.70,r'$\Delta l=\chi_3 R\alpha\zeta^2$',fontsize=11);ax.text(1.1,.38,'Fitted usable-length loss',fontsize=9);ax.text(-.28,-.12,'Joint',fontsize=9)
fig.subplots_adjust(left=.01,right=.99,top=.99,bottom=.01);fig.savefig(F/'30-torque/bpa_wrapping_loss.pdf');plt.close(fig)
for rel in ['15-background/moment_arm_methods.pdf','30-torque/bpa_wrapping_loss.pdf','15-background/steeleknee.pdf','30-torque/bktFrame.pdf']:
 d=fitz.open(F/rel);d[0].get_pixmap(matrix=fitz.Matrix(1.6,1.6)).save(S/(Path(rel).stem+'_after.png'))
print('Saved native-vector method/wrapping diagrams, relabeled PDFs, and split model panels')

