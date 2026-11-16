"""Measure native PF timing against the locked syn6 primitive reference."""
import argparse
from copy import deepcopy
import hashlib
import json
import os
from pathlib import Path
import time

for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'

import numpy as np

import motor_adapter as motor
from controller_profiles_20261104.native_pf import motor_digest, motor_rows
from pf_timing_audit import primitive_target
from syn6_s5_timing import measure


def rises(t, voltage):
    events = []
    armed = False
    for i, value in enumerate(voltage):
        if value < 1.5:
            armed = True
        if armed and value >= 3.5:
            events.append(t[i])
            armed = False
    return np.asarray(events)


def assess(t, voltage, names, manifest, anchor):
    fsa_path = Path(__file__).resolve().parent.parent / 'fsa_results/fsa_backsolve.npz'
    phase = np.arange(200) / 200.
    period = manifest['target_period_s']
    results = {}
    with np.load(fsa_path) as fsa:
        for side in ('l', 'r'):
            edges = rises(t, voltage[:, names.index('RG_E_' + side)])
            pairs = list(zip(edges[:-1], edges[1:]))[-4:]
            periods = [float(end - start) for start, end in pairs]
            channels = {}
            for channel in range(1, 7):
                name = f'PF_S{channel}_{side}'
                expected = primitive_target(fsa, manifest, side,
                                            (phase + anchor[side]) * period, channel - 1)
                actual = np.clip(voltage[:, names.index(name)] / 5., 0., 1.)
                rows = []
                for start, end in pairs:
                    observed = np.interp(start + phase * (end - start), t, actual)
                    rows.append(measure(phase, 1., np.tile(observed[:, None], (1, 2)),
                                        np.tile(expected[:, None], (1, 2)), (0,))[0])
                    rows[-1]['side'] = side
                channels[name] = rows
            results[side] = dict(periods_s=periods, channels=channels,
                                 passed=sum(row['passed'] for rows in channels.values() for row in rows),
                                 checks=6 * len(pairs))
    return results


def plot_timing(folder, manifest, anchor, factors):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    phase = np.arange(200) / 200.
    fsa_path = Path(__file__).resolve().parent.parent / 'fsa_results/fsa_backsolve.npz'
    fig, axes = plt.subplots(6, len(factors), figsize=(3.7 * len(factors), 8),
                             sharex=True, sharey=True, squeeze=False)
    with np.load(fsa_path) as fsa:
        for column, factor in enumerate(sorted(factors)):
            with np.load(folder / f'phase_{factor}.npz') as data:
                t, volt, names = data['time'], data['volts'], data['names'].tolist()
            for side, color, style in (('l', '#0072B2', '--'), ('r', '#D55E00', '-')):
                events = rises(t, volt[:, names.index('RG_E_' + side)])
                start, end = events[-2:]
                for channel, ax in enumerate(axes[:, column], 1):
                    wanted = primitive_target(fsa, manifest, side,
                                              (phase + anchor[side]) * manifest['target_period_s'],
                                              channel - 1)
                    values = volt[:, names.index(f'PF_S{channel}_{side}')] / 5.
                    actual = np.clip(np.interp(start + phase * (end - start), t, values), 0., 1.)
                    ax.plot(100 * phase, wanted, color='0.65', ls=style, lw=1.,
                            label=f'{side.upper()} target')
                    ax.plot(100 * phase, actual, color=color, ls=style, lw=1.3,
                            label=f'{side.upper()} native')
                    ax.axhline(.2, color='0.5', ls=':', lw=.5)
                    if column == 0:
                        ax.set_ylabel(f'S{channel}')
                axes[0, column].set_title(f'{end - start:.3f} s per cycle')
    axes[0, len(factors) // 2].legend(loc='upper right', fontsize=7, ncols=2)
    fig.supxlabel('Phase from ipsilateral RG extensor onset (%)')
    fig.supylabel('Normalized PF output')
    fig.tight_layout()
    fig.savefig(folder / 'PF_TIMING.png', dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fit', type=Path, required=True)
    parser.add_argument('--design', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--body', type=Path)
    parser.add_argument('--after-full-handover', action='store_true')
    parser.add_argument('--calibration', type=Path, help='Reuse a frozen native SPEC.json')
    parser.add_argument('--plot', action='store_true')
    parser.add_argument('--period-factors', type=float, nargs='+', default=(1., .85, 1.15))
    args = parser.parse_args()
    if any(not np.isfinite(factor) or factor <= 0. for factor in args.period_factors):
        parser.error('Period factors must be positive and finite')
    parent = Path(json.loads((args.fit / 'MANIFEST.json').read_text())['fit'])
    raw = (parent / 'FINAL.json').read_bytes()
    check = motor.handle_request(dict(op='preflight', fit=str(parent),
                                      fit_sha256=hashlib.sha256(raw).hexdigest()))
    if not check['ready_for_optimization']:
        raise RuntimeError('Frozen readout reconstruction failed')
    baseline = motor.INSTANCE
    fixed = {**motor.FIXED, **json.loads((args.fit / 'FINAL.json').read_text())['parameters']}
    baseline.evaluate(fixed)
    digest = motor_digest(baseline.c, baseline.names)
    anchor = {}
    for side in ('l', 'r'):
        col = list(baseline.rg).index(baseline.net.idx['RG_E_' + side])
        events = rises(baseline.t, baseline.rg_trace[:, col])[-5:] / baseline.T
        anchor[side] = float(np.angle(np.exp(2j * np.pi * events).mean()) / (2 * np.pi) % 1)
    args.output.mkdir(parents=True, exist_ok=False)
    report = dict(reference_rg_extensor_phase=anchor, motor_sha256=digest,
                  target_period_s=baseline.T)
    if args.body:
        live_digest = (args.body / 'MOTOR_SHA256.txt').read_text().strip()
        if live_digest != digest:
            raise RuntimeError('Closed-loop MN inputs differ from the frozen replay')
        with np.load(args.body / 'capture' / 'aligned_trace.npz') as data:
            t = data['t_neural']
            keep = np.ones(t.shape, dtype=bool)
            if args.after_full_handover:
                cadence = json.loads((args.body / 'CADENCE.json').read_text())
                start = cadence['full_native_since_neural_s']
                if start is None:
                    raise ValueError('Trial did not reach full native handover')
                keep = t >= start
                report['full_native_since_neural_s'] = start
            report['closed_loop'] = assess(t[keep], data['neural_voltage'][keep],
                                           data['neural_names'].tolist(), baseline.manifest, anchor)
            report['capture_duration_s'] = float(np.ptp(t[keep]))
        if (args.body / 'STOP.json').exists():
            report['stop'] = json.loads((args.body / 'STOP.json').read_text())
    else:
        spec = (json.loads(args.calibration.read_text()) if args.calibration else
                dict(period_s=baseline.T, designs=json.loads(args.design.read_text())['designs']))
        (args.output / 'SPEC.json').write_text(json.dumps(spec, indent=2))
        manifest = deepcopy(baseline.manifest)
        manifest['controller_profile']['native_pf'] = spec
        replay = type(baseline)(manifest)
        replay.phase = 'D'
        if hasattr(replay, 'feedforward_only'):
            replay.feedforward_only = False
        old_weights = {(baseline.names[p], baseline.names[q]): baseline.base_g[i]
                       for i, (p, q) in enumerate(zip(baseline.c.pre, baseline.c.post))}
        for i, (p, q) in enumerate(zip(replay.c.pre, replay.c.post)):
            endpoints = replay.names[p], replay.names[q]
            if endpoints in old_weights:
                replay.base_g[i] = old_weights[endpoints]
        cols = np.r_[replay.rg, replay.pf]
        names = [replay.names[i] for i in cols]
        relay_ids = [i for i, name in enumerate(replay.names) if name.startswith('Timing')]
        relay_trace = np.empty((len(replay.t), len(relay_ids)))
        step = 0
        forward = replay.c.forward

        def capture(inputs):
            nonlocal step
            value = forward(inputs)
            relay_trace[step] = replay.core.V[relay_ids]
            step += 1
            return value

        replay.c.forward = capture
        original_inputs = replay.u.copy()
        membrane = replay.core.time_factor_membrane.copy()
        original_cells = np.array([not name.startswith('Timing') for name in replay.names])
        report['offline'] = {}
        report['cadence_protocol'] = ('Scale legacy-cell integration and intrinsic RG h-gate time '
                                      'constants together with replay afferent cadence. Keep native '
                                      'timing relay parameters and all MN synapses fixed.')
        for speed in args.period_factors:
            step = 0
            replay.core.time_factor_membrane[original_cells] = membrane[original_cells] / speed
            sample = (replay.t / speed) % baseline.T
            if speed == 1.:
                np.copyto(replay.u, original_inputs)
            else:
                for column in range(replay.u.shape[1]):
                    replay.u[:, column] = np.interp(sample, replay.t, original_inputs[:, column])
            params = dict(fixed, tau_scale=fixed.get('tau_scale', 1.) * speed)
            started = time.perf_counter()
            diagnostics = replay.evaluate(params)['diagnostics']
            if motor_digest(replay.c, replay.names) != digest:
                old = set(motor_rows(baseline.c, baseline.names))
                new = set(motor_rows(replay.c, replay.names))
                raise RuntimeError(f'Native build altered MN inputs: {sorted(old-new)[:3]} vs {sorted(new-old)[:3]}')
            elapsed = time.perf_counter() - started
            if step != len(replay.t):
                raise RuntimeError(f'Neural run stopped after {step} steps')
            volts = np.column_stack((replay.rg_trace, replay.pf_trace))
            report['offline'][str(speed)] = dict(
                checks=assess(replay.t, volts, names, manifest, anchor),
                elapsed_s=elapsed, diagnostics=diagnostics)
            np.savez_compressed(args.output / f'phase_{speed}.npz', time=replay.t,
                                names=names, volts=volts, relay_volts=relay_trace[:step],
                                relay_names=[replay.names[i] for i in relay_ids])
        if args.plot:
            plot_timing(args.output, manifest, anchor, args.period_factors)
    (args.output / 'TIMING.json').write_text(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
