"""Plot a short Stage 2 recording and phase-align its joint motion to walking IK."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

from kinematic_response import GRID, JOINTS, compare, reference
from native_pf_check import rises


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    body = json.loads((args.run / 'body.result.json').read_text())
    with np.load(args.run / 'capture/aligned_trace.npz') as z:
        t, tn = z['t_state'], z['t_neural']
        names = z['neural_names'].tolist()
        joint_names = z['joint_names'].tolist()
        addresses = [z['joint_qpos_addresses'][joint_names.index(joint)] for joint in JOINTS]
        angles = np.degrees(z['qpos'][:, addresses])
        volts = z['neural_voltage']
        bw = body['measured']['body_weight_N']
        load = z['foot_force_anatomical_N'][:, :, 1] / bw
    cadence = json.loads((args.run / 'CADENCE.json').read_text())
    start = max(7., cadence['full_native_since_neural_s'] or 0.)
    sample = (t >= start) & (t <= 14.)
    dt = np.median(np.diff(t))
    peaks = []
    for col in (1, 4):
        values = -angles[:, col]
        candidates, _ = find_peaks(values, distance=round(.48 / dt), prominence=12.)
        peaks.append(t[candidates[sample[candidates]]])
    edges = peaks[0]
    if len(edges) < 3 or len(peaks[1]) < 2:
        raise ValueError('Insufficient complete bilateral knee cycles for an overlay')
    rg_edges = rises(tn, volts[:, names.index('RG_E_r')])
    rg_edges = rg_edges[(rg_edges >= start) & (rg_edges <= 14.)]
    plot_edges = rg_edges[:4] if len(rg_edges) >= 3 else edges[:4]
    first, last = plot_edges[0], plot_edges[-1]
    fig, axes = plt.subplots(5, 2, figsize=(10, 8), sharex=True, sharey='row')
    for col, side in enumerate(('l', 'r')):
        joint_side = 1 if side == 'l' else 0
        for row, keys, labels in (
            (0, ('RG_E', 'RG_F'), ('Extensor', 'Flexor')),
            (1, ('PF_S2', 'PF_S3'), ('S2', 'S3')),
            (2, ('MN_vas_int', 'MN_bifemsh'), ('Vastus intermedius', 'Biceps femoris SH')),
        ):
            ax = axes[row, col]
            for key, label, color, style in zip(keys, labels, ('#D55E00', '#0072B2'), ('-', '--')):
                ax.plot(tn, volts[:, names.index(key + '_' + side)], color=color,
                        ls=style, lw=1.15, label=label)
            ax.legend(loc='upper right', fontsize=7)
        axes[3, col].plot(t, angles[:, 1 + 3 * joint_side], color='#333333', lw=1.2)
        axes[4, col].plot(t, load[:, joint_side], color='#009E73', lw=1.2)
        axes[4, col].set_ylim(-.02, .10)
        axes[0, col].set_title('Left' if col == 0 else 'Right')
        axes[4, col].set_xlabel('Time (s)')
        for ax in axes[:, col]:
            ax.set_xlim(first, last)
            for event in plot_edges[1:-1]:
                ax.axvline(event, color='.8', lw=.65, zorder=0)
    for ax, label in zip(axes[:, 0], ('RG (mV)', 'PF (mV)', 'MN (mV)',
                                     'Knee (°)', 'Vertical foot load (BW)')):
        ax.set_ylabel(label)
    fig.suptitle('Syn6 Stage 2 · Live afferents · Suspended air stepping')
    fig.tight_layout()
    fig.savefig(args.output / 'STAGE2_TRACES.png', dpi=170)
    plt.close(fig)

    cycles = np.array([np.column_stack([np.interp(a + GRID * (b - a), t, angles[:, i])
                                       for i in range(6)]) for a, b in zip(edges[:-1], edges[1:])])
    periods = [float(np.median(np.diff(events))) for events in peaks]
    burst_counts = {}
    for side in ('r', 'l'):
        events = rises(tn, volts[:, names.index('RG_E_' + side)])
        events = events[(events >= start) & (events <= 14.)]
        for channel in ('S2', 'S3'):
            voltage = volts[:, names.index(f'PF_{channel}_{side}')]
            onsets = tn[1:][(voltage[:-1] < 1.) & (voltage[1:] >= 1.)]
            burst_counts[f'{channel}_{side}'] = [
                int(np.count_nonzero((onsets >= a) & (onsets < b)))
                for a, b in zip(events[:-1], events[1:])]
    gait_refs = Path(__file__).resolve().parents[1] / 'gait_refs'
    refs = [reference(), reference(gait_refs / 'ong_speed_125.npz')]
    fig, axes = plt.subplots(4, 2, figsize=(10, 7.8), sharex=True, sharey='row',
                             gridspec_kw={'height_ratios': [1, 1, 1, .48]})
    for col, side in enumerate(('Left', 'Right')):
        joint_side = 1 if side == 'Left' else 0
        axes[0, col].set_title(side)
        for row, joint in enumerate(('Hip', 'Knee', 'Ankle')):
            ax = axes[row, col]
            i = 3 * joint_side + row
            for curve in cycles[:, :, i]:
                ax.plot(100 * GRID, curve, color='#0072B2', alpha=.18, lw=.8)
            ax.plot(100 * GRID, cycles[:, :, i].mean(axis=0), color='#0072B2', lw=1.8,
                    label='Air-stepping mean')
            for ref, label, color, style in zip(
                    refs, ('Subject01 walking', 'Ong et al. (2019) simulation, 1.25 m/s'),
                    ('#222222', '#D55E00'), ('--', ':')):
                ax.plot(100 * GRID, ref['cycle_deg'][:, i], color=color, ls=style, lw=1.5, label=label)
            if col == 0:
                ax.set_ylabel(f'{joint} (°)')
        phase_ax = axes[3, col]
        for j, (ref, color) in enumerate(zip(refs, ('#333333', '#D55E00'))):
            duty = ref['contact_duty'][joint_side]
            if duty is None:
                continue
            loading = ref['reference_left_loading_phase'] if joint_side else 0.
            onset = (loading - ref['right_knee_anchor_fraction']) % 1.
            y = 1 - j
            phase_ax.broken_barh([(0., 100.)], (y - .26, .52), facecolors='.88')
            spans = [(100 * onset, 100 * min(duty, 1. - onset))]
            if onset + duty > 1.:
                spans.append((0., 100 * (onset + duty - 1.)))
            phase_ax.broken_barh(spans, (y - .26, .52), facecolors=color)
        phase_ax.set_yticks([0, 1], ['Ong et al.', 'Subject01'])
        phase_ax.tick_params(labelleft=(col == 0), length=0)
        phase_ax.set_xlim(0., 100.)
        phase_ax.set_ylim(-.5, 1.5)
        phase_ax.set_xlabel('Cycle from right knee max flexion (%)')
    axes[0, 0].legend(fontsize=7)
    fig.suptitle('Suspended syn6 motion versus ground-walking IK · Shape comparison')
    fig.text(.5, .015, 'Benchmark bars  ·  Colored = stance  ·  Gray = swing  ·  Air trial has no stance',
             ha='center', fontsize=9)
    fig.tight_layout(rect=(0., .04, 1., 1.))
    fig.savefig(args.output / 'KINEMATIC_OVERLAY.png', dpi=170)
    plt.close(fig)
    results = dict(run=str(args.run), plotted_window_s=[float(first), float(last)],
                   cycle_window_s=[float(edges[0]), float(edges[-1])],
                   right_knee_cycles=len(cycles), simulated_median_period_s=periods,
                   pf_onset_counts_by_rg_cycle=burst_counts,
                   peak_vertical_foot_load_bw=load.max(axis=0).tolist(),
                   benchmarks=[dict(identity=ref['identity'], period_s=ref['period_s'],
                                    **compare(cycles, periods, ref)[0]) for ref in refs])
    (args.output / 'OVERLAY_METRICS.json').write_text(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
