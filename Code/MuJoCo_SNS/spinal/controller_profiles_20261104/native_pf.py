"""Graded SNS PF relays with optional RG-driven clock-rate scheduling."""
import hashlib
import json
import math
import numpy as np

from sns_toolbox.connections import NonSpikingSynapse
from sns_toolbox.neurons import NonSpikingNeuron


CLOCK_TAU = 1.4
CLOCK_CEILING = 10.


def motor_rows(compiler, names):
    return sorted((names[pre], names[post], float(gain), float(lo), float(span),
                   float(reversal))
                  for pre, post, gain, lo, span, reversal in zip(
                      compiler.pre, compiler.post, compiler.gmax, compiler.elo,
                      compiler.erange, compiler.reversal)
                  if names[post].startswith('MN_'))


def motor_digest(compiler, names):
    return hashlib.sha256(json.dumps(motor_rows(compiler, names)).encode()).hexdigest()


def install(net, spec, connections):
    period = float(spec['period_s'])
    if not .5 <= period <= 3.:
        raise ValueError('Invalid native PF design period')
    designs = spec['designs']
    if set(designs) != {str(i) for i in range(1, 7)}:
        raise ValueError('Native syn6 requires six PF designs')

    def cell(name, tau, bias=0., initial=0.):
        net.idx[name] = len(net.net.populations)
        net.net.add_neuron(NonSpikingNeuron(membrane_capacitance=tau, bias=bias),
                           name=name, initial_value=initial)
        return name

    def edge(source, target, gain, reversal, lo, hi):
        synapse = NonSpikingSynapse(max_conductance=float(gain), reversal_potential=float(reversal),
                                   e_lo=float(lo), e_hi=float(hi))
        net.net.add_connection(synapse, source, target)

    removed = []
    for connection, source, target in list(connections(net)):
        if ((source.startswith('RG_') and target.startswith('PF_S'))
                or (source.startswith('PF_S') and target.startswith('PF_IN'))):
            if not spec.get('adaptive'):
                net.net.connections.remove(connection)
            removed.append((source, target))

    # A fast trigger receives excitation from the RG and delayed inhibition
    # from its follower. The resulting pulse resets a leaky ramp neuron.
    for side in ('l', 'r'):
        clocks = []
        for phase in ('E', 'F'):
            rg = f'RG_{phase}_{side}'
            follower = cell(f'TimingFollower_{phase}_{side}', .135)
            trigger = cell(f'TimingTrigger_{phase}_{side}', .009)
            clock = cell(f'TimingClock_{phase}_{side}', CLOCK_TAU, CLOCK_CEILING, CLOCK_CEILING)
            edge(rg, follower, 2., 8., 2.8, 3.8)
            edge(rg, trigger, 2., 8., 2.8, 3.8)
            edge(follower, trigger, 4., -4., .5, 3.5)
            edge(trigger, clock, 360., -.1, .8, 2.)
            clocks.append(clock)

        def threshold(phase, delay):
            reset_lag = spec.get('reset_lag_s', {}).get(str(phase), .045)
            return CLOCK_CEILING * -math.expm1(-max(.003, delay * period - reset_lag) / CLOCK_TAU)

        for channel, design in designs.items():
            quiet = cell(f'TimingQuiet_S{channel}_{side}', .032, 4., 4.)
            pf = f'PF_S{channel}_{side}'
            edge(quiet, pf, 2.8, 0., .5, 3.5)
            phases = design['triggers']
            values = design['parameters']
            if len(values) != 5 * len(phases) or any(phase not in (0, 1) for phase in phases):
                raise ValueError('Invalid native burst design')
            for burst, phase in enumerate(phases):
                onset, offset, rise, fall, gain = values[5 * burst:5 * burst + 5]
                if not 0 < onset < offset < 1 or min(rise, fall, gain) <= 0:
                    raise ValueError('Invalid native burst window')
                stem = f'S{channel}_{burst}_{side}'
                window = cell(f'TimingOn_{stem}', .012)
                shutoff = cell(f'TimingOff_{stem}', .01)
                lo = threshold(phase, onset)
                hi = max(lo + .025, threshold(phase, onset + rise))
                stop = threshold(phase, offset)
                edge(clocks[phase], window, 2., 8., lo, hi)
                edge(clocks[phase], shutoff, 2., 8., stop, threshold(phase, offset + fall))
                edge(shutoff, window, 6., -4., .15, 3.)
                edge(window, pf, gain, 8., .25, 5.3)
                edge(window, quiet, 3., -2., .35, 2.5)
        protection = spec.get('knee_premotor')
        if protection:
            relay = cell(f'TimingKneePremotor_{side}', .018)
            net.net.add_input(relay)
            net.inputs.append(f'KNEE_PREMOTOR_{side}')
            edge(relay, f'PF_S3_{side}', protection['flexor_gain'], 8., .25, 4.)
            edge(relay, f'TimingQuiet_S3_{side}', 3., -2., .25, 4.)
            edge(relay, f'PF_S2_{side}', protection['extensor_inhibition'], -3., .25, 4.)
            if spec.get('single_s3_phase'):
                hold = cell(f'TimingS3Hold_{side}', .018)
                net.net.add_input(hold)
                net.inputs.append(f'S3_HOLD_{side}')
                edge(hold, f'PF_S3_{side}', protection['flexor_gain'], 8., .25, 4.)
                edge(hold, f'TimingQuiet_S3_{side}', 3., -2., .25, 4.)
    net._n_true = len(net.net.populations)
    net.native_pf_removed = removed


class Cadence:
    """Use observed RG recurrence to set native clock rates and stage the handover."""

    def __init__(self, net, spec, dt):
        self.compiler = net.compiled
        self.dt = dt
        self.nominal = float(spec['period_s'])
        self.rg = np.array([net.idx[f'RG_{phase}_{side}']
                            for side in ('l', 'r') for phase in ('E', 'F')])
        self.clocks = np.array([net.idx[f'TimingClock_{phase}_{side}']
                                for side in ('l', 'r') for phase in ('E', 'F')])
        self.membrane = self.compiler.c.time_factor_membrane
        self.clock_base = self.membrane[self.clocks].copy()
        names = [row['name'] for row in net.net.populations]
        legacy = set(net.native_pf_removed)
        self.legacy_slots, self.relay_slots = [], []
        for slot, (pre, post) in enumerate(zip(self.compiler.pre, self.compiler.post)):
            source, target = names[pre], names[post]
            if (source, target) in legacy:
                self.legacy_slots.append(slot)
            if source.startswith('Timing') and target.startswith('PF_S'):
                self.relay_slots.append(slot)
        self.legacy_base = self.compiler.gmax[self.legacy_slots].copy()
        self.relay_base = self.compiler.gmax[self.relay_slots].copy()
        self.clock_rates = self.clock_base * self.nominal
        self.clock_work = np.empty_like(self.clock_rates)
        self.legacy_work = np.empty_like(self.legacy_base)
        self.relay_work = np.empty_like(self.relay_base)
        self.latch_s2 = bool(spec.get('latch_s2_until_rg', False))
        self.s3_phase = spec.get('single_s3_phase')
        self.s3_hold_ports = []
        self.s3_old_edges = []
        if self.s3_phase:
            start, end = self.s3_phase
            if not 0. <= start < end <= .9 or not spec.get('knee_premotor'):
                raise ValueError('Invalid single S3 phase window')
            self.s3_hold_ports = [net.inputs.index(f'S3_HOLD_{side}') for side in ('l', 'r')]
            for slot, (pre, post) in enumerate(zip(self.compiler.pre, self.compiler.post)):
                if names[pre].startswith(('TimingOn_S3_', 'TimingKneePremotor_')) and names[post].startswith(
                        ('PF_S3_', 'TimingQuiet_S3_')):
                    self.s3_old_edges.append(slot)
        self.s2_cells = [net.idx[f'PF_S2_{side}'] for side in ('l', 'r')]
        self.s2_windows = [
            [row for row, slot in enumerate(self.relay_slots)
             if self.compiler.post[slot] == cell
             and names[self.compiler.pre[slot]].startswith('TimingOn_')]
            for cell in self.s2_cells]
        self.startup = bool(spec.get('legacy_startup', True))
        release_s = float(spec.get('knee_premotor', {}).get('release_s', 0.))
        if not 0. <= release_s <= .3:
            raise ValueError('Premotor release must be between 0 and 300 ms')
        self.knee_release_decay = math.exp(-dt / release_s) if release_s else 0.
        self.knee_ports = ([net.inputs.index(f'KNEE_PREMOTOR_{side}') for side in ('l', 'r')]
                           if release_s or self.latch_s2 else [])
        self.knee_drive = np.zeros(len(self.knee_ports))
        self.shared_leg_clock = bool(spec.get('shared_leg_clock', False))
        self.period_alpha = float(spec.get('period_alpha', 1.))
        if not 0. < self.period_alpha <= 1.:
            raise ValueError('RG period update fraction must be in (0, 1]')
        self.flexion_phase = spec.get('reference_flexion_phase')
        if self.flexion_phase is not None and not .3 <= self.flexion_phase <= .8:
            raise ValueError('Invalid reference RG flexion phase')
        self.flexion_edges = []
        for clock in self.clocks[1::2]:
            slots = np.flatnonzero(self.compiler.pre == clock)
            self.flexion_edges.append((slots, self.compiler.elo[slots].copy(),
                                       self.compiler.erange[slots].copy()))
        self.handover_s = self.nominal * float(spec.get('handover_cycles', 1.))
        if not math.isfinite(self.handover_s) or self.handover_s <= 0.:
            raise ValueError('PF handover duration must be positive and finite')
        self.forward = self.compiler.forward
        self.compiler.forward = self.step
        self.reset()

    def reset(self):
        self.knee_drive.fill(0.)
        self.s2_latched = [False, False]
        self.armed = [True] * 4
        self.last = [None] * 4
        self.period = np.full(4, self.nominal)
        self.measured = [False] * 4
        self.time = 0.
        self.mix = 0. if self.startup else 1.
        self.full_native_since_s = None if self.startup else 0.
        for slots, lo, span in self.flexion_edges:
            self.compiler.elo[slots] = lo
            self.compiler.erange[slots] = span

    def align_flexion(self, i):
        extensor_onset = self.last[i - 1]
        if extensor_onset is None or not self.measured[i - 1]:
            return
        measured_phase = (self.time - extensor_onset) / self.period[i - 1]
        if not .2 <= measured_phase <= .9:
            return
        correction = float(np.clip(self.flexion_phase - measured_phase, -.18, .18))
        slots, lo, span = self.flexion_edges[i // 2]
        decay = math.exp(-correction * self.nominal / CLOCK_TAU)
        self.compiler.elo[slots] = np.maximum(.02, CLOCK_CEILING - (CLOCK_CEILING - lo) * decay)
        self.compiler.erange[slots] = span * decay

    def step(self, inputs):
        for side, port in enumerate(self.knee_ports):
            drive = max(float(inputs[port]), self.knee_release_decay * self.knee_drive[side])
            self.knee_drive[side] = drive
            inputs[port] = drive
        voltage = self.compiler.V
        for i, cell in enumerate(self.rg):
            if voltage[cell] < 1.5:
                self.armed[i] = True
            elif self.armed[i] and voltage[cell] >= 3.5:
                self.armed[i] = False
                if self.last[i] is not None:
                    interval = self.time - self.last[i]
                    if .55 * self.nominal <= interval <= 1.65 * self.nominal:
                        if self.measured[i] and self.period_alpha < 1.:
                            self.period[i] += self.period_alpha * (interval - self.period[i])
                        else:
                            self.period[i] = interval
                        self.measured[i] = True
                self.last[i] = self.time
                if i % 2 == 0:
                    self.s2_latched[i // 2] = False
                if i % 2 and self.flexion_phase is not None:
                    self.align_flexion(i)
        np.divide(self.clock_rates, self.period, out=self.clock_work)
        if self.shared_leg_clock:
            self.clock_work[1::2] = self.clock_work[::2]
        self.membrane[self.clocks] = self.clock_work
        if all(self.measured):
            self.mix = min(1., self.mix + self.dt / self.handover_s)
            if self.mix == 1. and self.full_native_since_s is None:
                self.full_native_since_s = self.time
        np.multiply(self.legacy_base, 1. - self.mix, out=self.legacy_work)
        np.multiply(self.relay_base, self.mix, out=self.relay_work)
        if self.latch_s2 and self.mix == 1.:
            for side, cell in enumerate(self.s2_cells):
                if self.knee_drive[side] >= 2. and voltage[cell] >= 1.:
                    self.s2_latched[side] = True
                if self.s2_latched[side]:
                    self.relay_work[self.s2_windows[side]] = 0.
        self.compiler.gmax[self.legacy_slots] = self.legacy_work
        self.compiler.gmax[self.relay_slots] = self.relay_work
        if self.s3_phase and self.mix == 1.:
            self.compiler.gmax[self.s3_old_edges] = 0.
            start, end = self.s3_phase
            for side, port in enumerate(self.s3_hold_ports):
                i = 2 * side
                age = (self.time - self.last[i]) / self.period[i]
                inputs[port] = 5. if start <= age < end else 0.
        self.time += self.dt
        return self.forward(inputs)
