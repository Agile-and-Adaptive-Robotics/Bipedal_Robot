"""Build the Gait2392 comparison figure from its cycle table.

Pass --before with the frozen-controller trajectory to regenerate the table.
An optional --after trajectory adds a fitted controller using the same reference.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
STEM = 'gait2392_ik_overlay'
JOINTS = ('hip', 'knee', 'ankle')


def reduce_cycles(path, tag, reference, kinematics, columns):
    with np.load(path, allow_pickle=False) as trace:
        n = np.count_nonzero(trace['t'] > 0) + 1
        t, q, neural = (trace[k][:n] for k in ('t', 'q', 'neuro'))
        # This runner sends each contact to both sensory channels. The archived
        # 'contact' array is their sum. Both channels carry the same normal force.
        force = trace['upward_N'][:n] if 'upward_N' in trace else trace['contact'][:n] * .5
    force_file = path.with_name('foot_ground.npz')
    if force_file.is_file():
        with np.load(force_file, allow_pickle=False) as measured:
            np.testing.assert_array_equal(measured['t'][:n], t)
            force = measured['upward_N'][:n]
    in_window = t >= 2.
    ik_time, _, _ = kinematics._read_mot(kinematics.IK_MOT)
    grf_time, grf_names, grf_data = kinematics._read_mot(kinematics.GRF_MOT)
    result = {}
    for col, side in enumerate(('r', 'l')):
        curves, periods, duty, _, used, count, onsets = kinematics.sim_side(t, q, neural, 2., side, force)
        if not used or count < 2:
            raise ValueError(f'Insufficient foot-loading cycles on side {side}')
        indices = (3, 4, 5) if side == 'r' else (8, 9, 10)
        _, _, cycles = kinematics._cycles_from_onsets(t[in_window], t[in_window], onsets, indices, q[in_window])
        occupancy = []
        for start, end in zip(onsets[:-1], onsets[1:]):
            keep = (t >= start) & (t < end)
            if .35 <= end-start <= 2.5 and keep.sum() >= 8:
                phase = 100*(t[keep]-start)/(end-start)
                occupancy.append(np.interp(kinematics.GRID, phase, (force[keep,col] > 50.).astype(float)))
        assert len(occupancy) == count
        columns[f'{side}_{tag}_stance_probability'] = np.mean(occupancy, axis=0)
        ref_name = 'ground_force_vy' if side == 'r' else '1_ground_force_vy'
        grf = grf_data[:, grf_names.index(ref_name)-1]
        a, b = kinematics._first_pair_in_ik(kinematics._loading_onsets(grf_time, grf), ik_time)
        columns[side+'_reference_stance'] = (np.interp(kinematics.GRID,100*(grf_time-a)/(b-a),grf)>50.).astype(float)
        errors = {}
        for joint in JOINTS:
            prefix = f'{side}_{joint}'
            columns[f'{prefix}_{tag}_mean_deg'] = curves[joint]
            columns[f'{prefix}_{tag}_sd_deg'] = np.std(cycles[joint],axis=0)
            columns[f'{prefix}_reference_deg'] = reference[side][joint]
            errors[joint] = float(np.sqrt(np.mean((curves[joint]-reference[side][joint])**2)))
        result[side] = dict(cycles=count, period_s=np.mean(periods), stance_fraction=duty, RMSE_deg=errors)
    return result


def plot(data, output):
    plt.rcParams.update({'font.family':'Arial','font.size':11.5,'pdf.fonttype':42,'axes.spines.top':False,
                         'axes.spines.right':False,'axes.linewidth':1.1,'axes.labelweight':'bold'})
    fig, panels=plt.subplots(4,2,figsize=(6.5,5.88),sharex=True,sharey='row',
                             gridspec_kw={'height_ratios':[1,1,1,.53]})
    has_fit='r_hip_fitted_mean_deg' in data.dtype.names
    x=data['cycle_percent']
    for col,side in enumerate(('r','l')):
        for row,joint in enumerate((*JOINTS,'stance')):
            ax=panels[row,col]
            if row<3:
                y=data[f'{side}_{joint}_model_mean_deg']
                sd=data[f'{side}_{joint}_model_sd_deg']
                ax.fill_between(x,y-sd,y+sd,color='#0072B2',alpha=.15,edgecolor='none')
                ref=data[f'{side}_{joint}_reference_deg']
                fit_key=f'{side}_{joint}_fitted_mean_deg'
                label=joint.capitalize()+' (deg)'
            else:
                y=data[f'{side}_model_stance_probability'];ref=data[f'{side}_reference_stance']
                fit_key=f'{side}_fitted_stance_probability';label='Cycles in\nstance (%)'
                y,ref=100*y,100*ref
                ax.set_ylim(-5,105);ax.set_yticks([0,50,100]);ax.set_xlabel('Gait cycle (%)')
            ax.plot(x,y,color='#0072B2',lw=1.55,label='Before tuning')
            if has_fit:
                yfit=data[fit_key]*(100 if row==3 else 1)
                ax.plot(x,yfit,color='#AA3377',ls='-.',lw=1.7,label='Stance fit')
            ax.plot(x,ref,color='black',ls='--',lw=1.55,label='Human recording')
            if col==0:ax.set_ylabel(label)
            if row==0:ax.set_title(('Right leg','Left leg')[col],fontweight='bold',fontsize=12)
            ax.set_xlim(0,100);ax.set_xticks([0,25,50,75,100]);ax.grid(alpha=.17)
            ax.tick_params(width=1.1,length=3)
    legend=fig.legend(*panels[0,0].get_legend_handles_labels(),loc='upper center',ncol=3 if has_fit else 2,
                      title='Legend',frameon=True,edgecolor='black',framealpha=1.,fancybox=False,
                      fontsize=10.5,title_fontsize=10.5,borderpad=.35)
    legend.get_title().set_fontweight('bold')
    fig.subplots_adjust(left=.15,right=.97,bottom=.085,top=.847,hspace=.20,wspace=.17)
    fig.savefig(output);plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--before',type=Path);p.add_argument('--after',type=Path)
    p.add_argument('--spinal',type=Path);args=p.parse_args()
    csv=HERE/(STEM+'.csv')
    if args.after and not args.before:p.error('--after requires --before')
    if args.before:
        default=next((root/'Code/MuJoCo_SNS/spinal' for root in HERE.parents if (root/'Code/MuJoCo_SNS/spinal/kine_ref.py').exists()),None)
        spinal=args.spinal or default
        if not spinal:p.error('Supply the spinal model directory with --spinal')
        sys.path.insert(0,str(spinal.resolve()));import kine_ref
        reference=kine_ref.load_reference()
        columns={'cycle_percent':kine_ref.GRID};metrics={};sources={}
        for tag,path in [('model',args.before),('fitted',args.after)]:
            if path is None:continue
            metrics[tag]=reduce_cycles(path,tag,reference,kine_ref,columns)
            sources[tag]={'trajectory_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
            recorded_load = path.with_name('foot_ground.npz')
            with np.load(path, allow_pickle=False) as trace:
                has_upward = 'upward_N' in trace
            sources[tag]['force_source'] = ('independent upward ground force' if has_upward or recorded_load.is_file()
                                            else 'half the archived sum of heel and toe channels')
            if recorded_load.is_file():
                sources[tag]['independent_force_sha256']=hashlib.sha256(recorded_load.read_bytes()).hexdigest()
        np.savetxt(csv,np.column_stack(list(columns.values())),delimiter=',',header=','.join(columns),comments='')
        (HERE/(STEM+'_metrics.json')).write_text(json.dumps(metrics,indent=2)+'\n',encoding='utf8')
        (HERE/(STEM+'_provenance.json')).write_text(json.dumps({'sources':sources,'normal_force_fraction_of_legacy_channel_sum':.5,'foot_load_threshold_N':50.,'start_s':2.,'cycle_interval_s':[.35,2.5],'reference_files':[kine_ref.IK_MOT.name,kine_ref.GRF_MOT.name]},indent=2)+'\n',encoding='utf8')
    plot(np.genfromtxt(csv,delimiter=',',names=True),HERE.parent/(STEM+'.pdf'))


if __name__=='__main__':
    main()
