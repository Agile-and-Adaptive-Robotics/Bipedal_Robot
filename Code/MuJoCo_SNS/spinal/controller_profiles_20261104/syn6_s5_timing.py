"""Experimental bilateral PF onset and shutoff relays on the frozen D replay."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import time

os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
from numba import njit
from scipy.optimize import differential_evolution

import motor_adapter as motor
from pf_timing_audit import features, primitive_target


@njit
def relay_clocks(rg, dt, nominal_period):
    """Schmitt-triggered phase timers, updated from measured RG transitions."""
    clocks = np.empty_like(rg)
    for channel in range(rg.shape[1]):
        armed = True
        last = -1
        period = nominal_period
        for i in range(len(rg)):
            voltage = rg[i, channel]
            if voltage < 1.5:
                armed = True
            if armed and voltage >= 3.5:
                if last >= 0:
                    period = (i - last) * dt
                last = i
                armed = False
            clocks[i, channel] = (i - last) * dt / period if last >= 0 else 2.0
    return clocks


@njit
def conductances(clocks, parameters, triggers):
    excitation = 0.0
    opening = 0.0
    for phase in range(len(triggers)):
        onset, offset, rise, fall, gain = parameters[phase * 5:phase * 5 + 5]
        delay = clocks[triggers[phase]]
        on = min(1.0, max(0.0, (delay - onset) / rise))
        off = min(1.0, max(0.0, (delay - offset) / fall))
        window = on * (1.0 - off)
        excitation += gain * window
        opening = max(opening, window)
    return excitation, 1.5 * (1.0 - opening)


@njit
def isolated_pf(clocks, parameters, triggers, membrane, leak, rest, bias):
    output = np.empty((len(clocks), 2))
    for side in range(2):
        voltage = rest[side]
        for i in range(len(clocks)):
            excite, shutoff = conductances(clocks[i, side * 2:side * 2 + 2], parameters, triggers)
            voltage += membrane[side] * (leak[side] * (rest[side] - voltage)
                                         + bias[side] + excite * (8.0 - voltage)
                                         - shutoff * voltage)
            output[i, side] = min(1.0, max(0.0, voltage / 5.0))
    return output


def measure(t, period, actual, wanted, cycles=range(7, 11)):
    phase = np.arange(200) / 200
    rows = []
    for side, label in enumerate(('l', 'r')):
        times = period * (np.asarray(list(cycles))[:, None] + phase)
        observed = np.interp(times, t, actual[:, side]).mean(axis=0)
        target = np.interp(times, t, wanted[:, side]).mean(axis=0)
        af, tf = features(observed), features(target)
        lag = (tf['fundamental_phase_deg'] - af['fundamental_phase_deg'] + 180) % 360 - 180
        errors = {}
        for edge in ('onset', 'offset'):
            key = edge + '_phase_percent'
            errors[edge] = None
            if len(af[key]) == len(tf[key]) and tf[key]:
                choices = [(np.roll(af[key], shift) - tf[key] + 50) % 100 - 50
                           for shift in range(len(tf[key]))]
                errors[edge] = min(choices, key=lambda delta: max(abs(delta))).tolist()
        passed = (all(value is not None and max(map(abs, value)) <= 5
                      for value in errors.values())
                  and abs(lag) <= 30 and abs(af['duty_fraction'] - tf['duty_fraction']) <= .15
                  and af['excursion'] >= .2)
        rows.append(dict(side=label, actual=af, target=tf, edge_errors_percent=errors,
                         fundamental_lag_deg=lag, max_error=float(np.max(abs(observed-target))),
                         passed=bool(passed)))
    return rows


def relay_design(t, period, clocks, target, channel):
    if channel == 5:
        return np.array([0, 1]), [(.36, .50), (.58, .73), (.03, .20), (.02, .14), (.5, 4),
                                  (.20, .33), (.34, .43), (.025, .14), (.02, .12), (.3, 3)]
    phase = np.arange(200) / 200
    reference = features(np.interp((phase + 6) * period, t, target[:, 1]))
    offsets = np.array(reference['offset_phase_percent']) / 100
    triggers, bounds = [], []
    for onset in np.array(reference['onset_phase_percent']) / 100:
        width = min((offsets - onset) % 1)
        ages = clocks[np.searchsorted(t, (6 + onset) * period), 2:4]
        candidates = [i for i, age in enumerate(ages) if .035 < age and age + width < .97]
        if not candidates:
            raise ValueError(f'S{channel} has no preceding RG event with a complete window')
        trigger = min(candidates, key=lambda i: ages[i])
        age = ages[trigger]
        triggers.append(trigger)
        bounds.extend([(max(0, age - .14), age + .025),
                       (age + width - .09, min(.98, age + width + .045)),
                       (.012, .18), (.012, .15), (.25, 5)])
    return np.array(triggers), bounds


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fit', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--verify', action='store_true')
    parser.add_argument('--channel', type=int, choices=range(1, 7), default=5)
    parser.add_argument('--bundle', type=Path, help='Apply all six saved channel calibrations together')
    args = parser.parse_args()
    final = json.loads((args.fit / 'FINAL.json').read_text())
    parent = Path(json.loads((args.fit / 'MANIFEST.json').read_text())['fit'])
    request = dict(op='preflight', fit=str(parent),
                   fit_sha256=hashlib.sha256((parent / 'FINAL.json').read_bytes()).hexdigest())
    if not motor.handle_request(request)['ready_for_optimization']:
        raise RuntimeError('Frozen readout reconstruction failed')
    replay = motor.INSTANCE
    fixed = {**motor.FIXED, **final['parameters']}
    original_diagnostics = replay.evaluate(fixed, time.time() + 120)['diagnostics']
    core = replay.c.c
    indices = [replay.names.index(f'PF_S{args.channel}_{side}') for side in ('l', 'r')]
    rg_indices = [replay.names.index('RG_' + phase + '_' + side)
                  for side in ('l', 'r') for phase in ('E', 'F')]
    rg_baseline = replay.rg_trace.copy()
    columns = [list(replay.pf).index(i) for i in indices]
    baseline = np.clip(replay.pf_trace[:, columns] / 5, 0, 1)
    mn_edges = np.isin(replay.c.post, replay.mn)
    frozen_weights = replay.c.gmax[mn_edges].copy()
    fsa_path = Path(__file__).resolve().parent.parent / 'fsa_results/fsa_backsolve.npz'
    with np.load(fsa_path) as fsa:
        target = np.column_stack([primitive_target(fsa, replay.manifest, side, replay.t, args.channel - 1)
                                  for side in ('l', 'r')])
    lagged_rg = np.vstack((replay.initial['V'][rg_indices], rg_baseline[:-1]))
    clocks = relay_clocks(lagged_rg, replay.dt, replay.T)
    membrane = core.time_factor_membrane[indices].copy()
    leak, rest, bias = (getattr(core, name)[indices].copy() for name in ('g_m', 'V_rest', 'i_b'))
    triggers, bounds = relay_design(replay.t, replay.T, clocks, target, args.channel)
    if args.bundle:
        if not args.verify:
            args.output.mkdir(parents=True, exist_ok=False)
        designs = {channel: json.loads((args.bundle / f'S{channel}' / 'RESULT.json').read_text())
                   for channel in range(1, 7)}
        if any(design['channel'] != channel or Path(design['source_fit']) != args.fit
               for channel, design in designs.items()):
            raise ValueError('Bundle contains a mismatched channel calibration')
    elif args.verify:
        parameters = np.array(json.loads((args.output / 'RESULT.json').read_text())['parameters'])
    else:
        args.output.mkdir(parents=True, exist_ok=False)
        training = (replay.t >= 2 * replay.T) & (replay.t < 7 * replay.T)

        def objective(values):
            prediction = isolated_pf(clocks, values, triggers, membrane, leak, rest, bias)
            error = prediction[training] - target[training]
            return np.mean(error ** 2) + .15 * np.mean((prediction[training] > .2) != (target[training] > .2))

        parameters = differential_evolution(objective, bounds, seed=17, popsize=7,
                                            maxiter=65, polish=False).x
    if not args.bundle:
        designs = {args.channel: dict(parameters=parameters.tolist(), triggers=triggers.tolist())}
    relay_rows = []
    for channel, design in designs.items():
        for side, label in enumerate(('l', 'r')):
            index = replay.names.index(f'PF_S{channel}_{label}')
            relay_rows.append((index, slice(2 * side, 2 * side + 2),
                               np.array(design['parameters']), np.array(design['triggers']),
                               float(core.i_b[index])))
    removed = np.isin(replay.c.post, [row[0] for row in relay_rows])
    assert not np.any(removed & mn_edges)
    forward = replay.c.forward
    elapsed = np.full(4, 2.0 * replay.T)
    periods = np.full(4, replay.T)
    armed = np.ones(4, dtype=bool)
    last = np.full(4, -1)
    step = 0

    def advance(inputs):
        nonlocal step
        elapsed[:] += replay.dt
        voltage = core.V[rg_indices]
        armed[voltage < 1.5] = True
        for event in np.flatnonzero(armed & (voltage >= 3.5)):
            if last[event] >= 0:
                periods[event] = (step - last[event]) * replay.dt
            last[event] = step
            elapsed[event] = 0
            armed[event] = False
        phase_age = elapsed / periods
        for index, channels, gains, sources, baseline_bias in relay_rows:
            excite, shutoff = conductances(phase_age[channels], gains, sources)
            core.i_b[index] = baseline_bias + excite * (8 - core.V[index]) - shutoff * core.V[index]
        replay.c.gmax[removed] = 0
        step += 1
        return forward(inputs)

    replay.c.forward = advance
    diagnostics = replay.evaluate(fixed, time.time() + 120)['diagnostics']
    assert np.array_equal(replay.c.gmax[mn_edges], frozen_weights)
    if args.bundle:
        checks = {}
        with np.load(fsa_path) as fsa:
            for channel in designs:
                cols = [list(replay.pf).index(replay.names.index(f'PF_S{channel}_{side}'))
                        for side in ('l', 'r')]
                actual = np.clip(replay.pf_trace[:, cols] / 5, 0, 1)
                target = np.column_stack([primitive_target(fsa, replay.manifest, side, replay.t, channel - 1)
                                          for side in ('l', 'r')])
                checks[str(channel)] = dict(
                    mean=measure(replay.t, replay.T, actual, target),
                    strides={str(cycle): measure(replay.t, replay.T, actual, target, (cycle,))
                             for cycle in range(7, 11)})
        result = dict(channels=checks, mn_weights_exact=True,
                      rg_trace_exact=bool(np.array_equal(rg_baseline, replay.rg_trace)),
                      baseline_mn_rmse=original_diagnostics['mn_command_rmse'],
                      mn_waveform_rmse=diagnostics['mn_command_rmse'],
                      scope='All 12 PF outputs driven simultaneously by hybrid RG-triggered relays in offline D replay')
        result['passed'] = all(row['passed'] for check in checks.values()
                               for rows in [check['mean'], *check['strides'].values()] for row in rows)
        result['designs'] = {str(channel): {key: design[key] for key in ('parameters', 'triggers')}
                             for channel, design in designs.items()}
        save_result(args, result)
        print(json.dumps({key: value for key, value in result.items() if key not in ('channels', 'designs')}))
        return
    actual = np.clip(replay.pf_trace[:, columns] / 5, 0, 1)
    rows = measure(replay.t, replay.T, actual, target)
    per_stride = {str(cycle): measure(replay.t, replay.T, actual, target, (cycle,))
                  for cycle in range(7, 11)}
    result = dict(channel=args.channel, triggers=triggers.tolist(), parameters=parameters.tolist(), sides=rows,
                  baseline=measure(replay.t, replay.T, baseline, target),
                  source_fit=str(args.fit), target_file=str(fsa_path), period_s=replay.T,
                  mn_weights_exact=True, rg_trace_exact=bool(np.array_equal(rg_baseline, replay.rg_trace)),
                  mn_waveform_rmse=diagnostics['mn_command_rmse'],
                  per_stride=per_stride,
                  passed=all(row['passed'] for row in rows)
                  and all(row['passed'] for stride in per_stride.values() for row in stride),
                  scope='Hybrid RG-triggered relay prototype coupled to one bilateral PF channel during offline D replay')
    save_result(args, result)
    if not args.verify:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 1, figsize=(9, 5), sharex=True, constrained_layout=True)
        mask = (replay.t >= 8 * replay.T) & (replay.t < 10 * replay.T)
        x = (replay.t[mask] / replay.T - 8) * 100
        for side, axis in enumerate(axes):
            for values, label, style in ((target, 'Primitive target', 'k--'),
                                         (baseline, 'Original PF', 'C1:'),
                                         (actual, 'Onset / shutoff relays', 'C0-')):
                axis.plot(x, values[mask, side], style, label=label)
            axis.axhline(.2, color='0.65', linewidth=.7)
            axis.set_ylabel(f'{("Left", "Right")[side]} S{args.channel} activity')
        axes[0].legend(ncol=3)
        axes[-1].set_xlabel('Reference cycle phase (%) across two strides')
        fig.savefig(args.output / f'S{args.channel}_TIMING.png', dpi=180)
    print(json.dumps(dict(passed=result['passed'], sides=rows, repeat=args.verify)))


def save_result(args, result):
    name = 'FRESH_REPEAT.json' if args.verify else 'RESULT.json'
    if args.verify:
        original = json.loads((args.output / 'RESULT.json').read_text())
        result['exact_repeat'] = result == original
        if not result['exact_repeat']:
            raise RuntimeError('Fresh replay differs')
    (args.output / name).write_text(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
