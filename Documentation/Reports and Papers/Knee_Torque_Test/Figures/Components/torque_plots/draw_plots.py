"""Draw the paper's torque plots at their printed sizes from the MATLAB records."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import FixedLocator, FuncFormatter, MultipleLocator
from matplotlib.transforms import Bbox, TransformedBbox
import numpy as np
from PIL import Image
from scipy.io import loadmat


HERE = Path(__file__).resolve().parent
FIGURES = HERE.parents[1]
PAPERS = FIGURES.parents[1]
REPO = PAPERS.parents[1]
DISS = PAPERS / 'Dissertation'
VALIDATION = DISS / 'Notes/figure_names_and_overview_20261119/validation_figures'
INK = '#212121'
BLUE = '#0000FF'
GOLD = '#FFD700'
SALMON = '#FA8775'
REDUCED_PRESSURE = '#EA5F94'
MAGENTA = '#CD34B5'
GRAY = '#666666'


def configure(fonts=Path('C:/Windows/Fonts')):
    for file in ['arial.ttf', 'arialbd.ttf', 'arialbi.ttf']:
        font_manager.fontManager.addfont(fonts / file)
    plt.rcParams.update({'font.family':'Arial', 'font.weight':'bold', 'font.size':11,
                         'axes.labelweight':'bold', 'axes.edgecolor':INK, 'axes.linewidth':2,
                         'text.color':INK, 'xtick.color':INK, 'ytick.color':INK,
                         'axes.unicode_minus':True, 'pdf.fonttype':42, 'pdf.use14corefonts':False,
                         'path.simplify':False, 'mathtext.fontset':'custom',
                         'mathtext.rm':'Arial:weight=bold', 'mathtext.it':'Arial:style=italic:weight=bold',
                         'mathtext.bf':'Arial:weight=bold', 'mathtext.default':'it'})


def canvas(width,height):
    return plt.figure(figsize=(width/72,height/72),facecolor='white')


def at(fig, x, y, value, **style):
    w,h=fig.get_size_inches()*72
    return fig.text(x/w, y/h, value, va='center', **style)


def angle_label(fig,x,y):
    at(fig,x-13,y,'θ',fontsize=12.5,fontstyle='italic')
    at(fig,x-5.3,y-3.1,'k',fontsize=10.2)
    at(fig,x+.9,y,', °',fontsize=12.5)


def axes(fig, box, xlim, ylim, xticks, yticks, xminor=5, yminor=1):
    w,h=fig.get_size_inches()*72
    x,y,aw,ah=box
    ax=fig.add_axes([x/w,y/h,aw/w,ah/h],xlim=xlim,ylim=ylim)
    ax.spines[['right','top']].set_visible(False)
    ax.xaxis.set_major_locator(FixedLocator(xticks))
    ax.yaxis.set_major_locator(FixedLocator(yticks))
    ax.xaxis.set_minor_locator(MultipleLocator(xminor))
    ax.yaxis.set_minor_locator(MultipleLocator(yminor))
    ax.tick_params(direction='in',width=2,length=5.6,labelsize=10.5,pad=3)
    ax.tick_params(which='minor',direction='in',width=2,length=2.8)
    return ax


def dots(ax,x,y,color,title,area=60,alpha=.75,scale=.72,edge='none',edge_width=0):
    return ax.scatter(x,y,s=np.asarray(area)*scale**2,marker='o',c=color,edgecolors=edge,
                      linewidths=edge_width,alpha=alpha,label=title,zorder=4)


def curve(ax,x,y,color,title,width,pattern):
    return ax.plot(x,y,color=color,label=title,linewidth=width,linestyle=pattern,
                   solid_capstyle='round',dash_capstyle='butt',zorder=2)[0]


def legend(ax,entries,loc, size=10.4, anchor=None):
    names=[entry.get_label() for entry in entries]
    key=ax.legend(entries,names,loc=loc,bbox_to_anchor=anchor,prop={'family':'Arial','weight':'bold','size':size},
                  markerscale=1,scatterpoints=1,numpoints=1,handlelength=1.28,handletextpad=.34,
                  borderaxespad=.33,borderpad=.30,labelspacing=.17,fancybox=False,
                  framealpha=1,edgecolor=INK,facecolor='white')
    key.get_frame().set_linewidth(.9)
    key.set_zorder(8)
    return key


def multiple_panels(out,number,records):
    six=number==9
    w,h=500,440
    fig=canvas(w,h)
    xs=[54,208,362] if six else [55,293]
    aw=131 if six else 196
    bottoms=[260,43]
    ah=141 if six else 140
    titles=(['43.6 cm','43.2 cm + tendon','45.7 cm','40.0 cm','40.0 cm + tendon','46.5 cm']
            if six else ['48.5 cm','41.0 cm','45.7 cm','40.6 cm + tendon'])
    if six:
        continuation=loadmat(HERE/'extensor_extension.mat',simplify_cells=True)['extension']
        panel_d=records[3]
        endpoint,=panel_d['measured_angle'][(panel_d['measured_angle']>0) & (panel_d['measured']==0)]
    for k,(t,title) in enumerate(zip(records,titles)):
        row,col=divmod(k,len(xs))
        x,y=xs[col],bottoms[row]
        if six:
            ax=axes(fig,(x,y,aw,ah),(-135,endpoint),(-2,20),[-120,-90,-60,-30,0,endpoint],range(0,21,5))
            ax.xaxis.set_major_formatter(FuncFormatter(lambda v,_: f'{v:g}'.replace('-','−')))
            ax.get_xticklabels()[-1].set_ha('right')
        else:
            bottom=-55 if number==8 else -45
            ax=axes(fig,(x,y,aw,ah),(-120,15),(bottom,0),[-120,-90,-60,-30,0],range(-50 if number==8 else -40,1,10),yminor=2)
        at(fig,1 if col==0 else x-23,y+ah+22,f'({chr(65+k)})',fontsize=12.5,ha='left')
        at(fig,x+aw/2,y+ah+22,title,fontsize=11.3,ha='center')
        at(fig,x+aw/2,y+ah+8,'Training' if row==0 else 'Validation',fontsize=10.6,ha='center')
        measured=dots(ax,t['measured_angle'],t['measured'],BLUE,'Measured')
        if number!=8:
            hybrid=dots(ax,t['hybrid_angle'],t['hybrid'],GOLD,'Hybrid')
            measured.set_zorder(5)
            if six:
                # A marker centered on the measured endpoint needs its outer radius.
                right_edge=TransformedBbox(Bbox.from_extents(0,0,1+3.2/aw,1),ax.transAxes)
                for markers in (measured,hybrid):
                    markers.set_clip_box(right_edge)
        if number==7:
            hunt=curve(ax,t['angle'],t['hunt'],GRAY,'Hunt equation',2,'--')
            base=curve(ax,t['angle'],t['unoptimized'],SALMON,'Unoptimized',2.5,'-.')
            if k==3:
                legend(ax,[hybrid,measured,hunt,base],'lower left')
        else:
            angle,baseline,corrected=(t[name] for name in ('angle','unoptimized','optimized'))
            if six:
                tail=continuation[k]
                assert tail['id']==t['id'] and tail['angle'][0]>angle[-1] and tail['angle'][-1]==endpoint
                angle=np.r_[angle,tail['angle']]
                baseline=np.r_[baseline,tail['unoptimized']]
                corrected=np.r_[corrected,tail['optimized']]
            base=curve(ax,angle,baseline,SALMON,'Unoptimized',2,'--')
            fitted=curve(ax,angle,corrected,MAGENTA,'Optimized',2.5,'-')
            if six and k==2:
                legend(ax,[hybrid,measured,base,fitted],'upper right')
            if not six and k==3:
                legend(ax,[measured,base,fitted],'lower left')
    at(fig,30,214,'Torque, N m',ha='center',rotation=90,fontsize=12)
    angle_label(fig,264,12)
    fig.savefig(out,metadata={'Title':f'Pinned-knee torque, Figure {number}','CreationDate':None})
    plt.close(fig)


def biomimetic(out):
    r=json.loads((FIGURES/'Components/publication_20261028/Ext10mm_current_audit.json').read_text())
    fig=canvas(452,218)
    image=np.asarray(Image.open(FIGURES/'Figure components/KneeExt_10mm_passthru.png').convert('RGB'))
    photo=fig.add_axes([22/452,19/218,178/452,191/218])
    photo.imshow(image,interpolation='antialiased')
    photo.set_axis_off()
    ax=axes(fig,(266,36,176,164),(-125,15),(0,20),[-120,-90,-60,-30,0],range(0,21,5))
    base=curve(ax,r['angle'],r['baseline'],SALMON,'Improved BPA',2,'--')
    corrected=curve(ax,r['angle'],r['optimized'],MAGENTA,'Optimized model',2.5,'-')
    measured=dots(ax,r['measured_angle'],r['measured'],BLUE,'Measured')
    legend(ax,[base,corrected,measured],'upper right')
    at(fig,5,201,'(A)',fontsize=13)
    at(fig,201,201,'(B)',fontsize=13)
    at(fig,231,119,'Torque, N m',ha='center',rotation=90,fontsize=12)
    angle_label(fig,354,12)
    fig.savefig(out,metadata={'Title':'Biomimetic-knee extensor test','CreationDate':None})
    plt.close(fig)


def validation(out,role,prop):
    fig=canvas(425,300.9)
    ext=role=='extensor'
    ax=axes(fig,(50,39,354,249),(-125,15),(20,110) if ext else (-105,5),
            [-120,-90,-60,-30,0],range(20,111,10) if ext else range(-100,1,10),yminor=2)
    if ext:
        path=REPO/'Code/Matlab/Mesh_Optimization/Results/Vas_Pam_20mm_Result_cam_20261014_0416.mat'
        model=loadmat(path,variable_names=['ctxCamBest','predBest'],simplify_cells=True)
        ctx=model['ctxCamBest']
        design_x,design_y=ctx['phiD'],model['predBest']['TorqueZ']
        human_x,human_y=ctx['humanAngleD'],ctx['humanTorque']
        trials=np.genfromtxt(VALIDATION/'extensor_trials.csv',delimiter=',',names=True,dtype=None,encoding='utf8')
        x=trials['K_deg']
        measured_y,pred_y=trials['measured_torque_Nm_RECORD'],trials['pred_at_measured_P_Nm']
        pressures=trials['Pressure_kPa']
        ma=pa=60*pressures/620
        pred_alpha=1
    else:
        values={a['name']:a for a in json.loads((VALIDATION/'flexor_data.json').read_text())['series']}
        design_x,design_y=values['design']['x'],values['design']['y']
        human_x,human_y=values['human']['x'],values['human']['y']
        x=values['measured']['x']
        assert x==values['prediction']['x']
        measured_y,pred_y=values['measured']['y'],values['prediction']['y']
        pressure_record=DISS/'Figures/30-torque/Components/publication_20261028/flexor_pressure_publication.json'
        recorded=json.loads(pressure_record.read_text())
        np.testing.assert_array_equal(recorded['angle'],x)
        np.testing.assert_allclose(recorded['prediction'],pred_y,rtol=0,atol=1e-9)
        pressures=np.asarray(recorded['pressure_kPa'])
        ma=prop[0]['size']; pa=prop[1]['size']; pred_alpha=prop[1]['alpha']
        assert all(s['edge']=='none' for s in prop)
    design=curve(ax,design_x,design_y,MAGENTA,'Predicted, 620 kPa',2.5,'-')
    human=curve(ax,human_x,human_y,'black','Vastus medialis' if ext else 'Biceps femoris short head',3,':')
    prediction=dots(ax,x,pred_y,GOLD,'Prediction at measured pressure',pa,pred_alpha,scale=.85)
    angles=np.asarray(x)
    torques=np.asarray(measured_y)
    areas=np.broadcast_to(ma,angles.shape)
    entries=[design,human,prediction]
    reduced=pressures<620
    for mask,color,title in [(~reduced,BLUE,'Measured'),(reduced,REDUCED_PRESSURE,'Measured at reduced pressure')]:
        if np.any(mask):
            entries.append(dots(ax,angles[mask],torques[mask],color,title,areas[mask],scale=.85))
    legend(ax,entries,'upper right' if ext else 'lower left',10.5,
           (.986,.987) if ext else (.026,.043))
    at(fig,12,161,'Torque, N m',ha='center',rotation=90,fontsize=12)
    angle_label(fig,227,12)
    fig.savefig(out,metadata={'Title':f'Physical validation of the two-BPA {role}','CreationDate':None})
    plt.close(fig)


if __name__=='__main__':
    args=argparse.ArgumentParser(description=__doc__)
    args.add_argument('directory',type=Path)
    dest=args.parse_args().directory
    dest.mkdir(parents=True,exist_ok=True)
    configure()
    flex_data=loadmat(HERE/'plot_data.mat',simplify_cells=True)
    for n,name in [(7,'FlxPin_preopt_publication'),(8,'FlxPin_optimized_publication')]:
        multiple_panels(dest/(name+'.pdf'),n,flex_data['flex'])
    ext=loadmat(DISS/'Figures/30-torque/Components/extensor_figures/extensor_data.mat',simplify_cells=True)
    multiple_panels(dest/'ExtPin_complete_publication.pdf',9,ext['tests'])
    biomimetic(dest/'Ext10mm_52cm_current_publication.pdf')
    validation(dest/'Ext20_current_publication.pdf','extensor',flex_data['props'])
    validation(dest/'Redesigned_flexor_pressure_publication.pdf','flexor',flex_data['props'])
