"""Read-only actual Stage5 contact-PF graph/current audit; no neural/body steps."""
import os
os.environ['AARL_NET'] = 'w2lvar'
import argparse
import json
from pathlib import Path
import numpy as np
import integrated_hunt_design as I
from prepare_renshaw_integration import construct


def row(n, edge):
    c = n.compiled
    return dict(source=n.names[c.pre[edge]], target=n.names[c.post[edge]],
                g_uS=float(c.gmax[edge]), E_mV=float(c.reversal[edge]),
                low_mV=float(c.elo[edge]), span_mV=float(c.erange[edge]))


def run(condition, out):
    cfg = json.loads(Path(condition).read_text())
    source = (I.H/'physical_eval.py').read_text()
    results = []
    for scale in (.25, 1.):
        base, _ = construct(source, cfg, 4);b = base.compiled
        n, _ = construct(source, {**cfg, 'contact_pf_scale': scale}, 5);c = n.compiled
        ids = np.array([i for i, (pre, post) in enumerate(zip(c.pre, c.post))
                        if n.names[post].startswith('PF_') and n.names[pre].startswith(('HEEL_', 'TOEDF_'))])
        if not np.array_equal(np.flatnonzero(b.gmax != c.gmax), ids):
            raise AssertionError('Stage5 changed something besides its declared contact-PF slots')
        expected = n.original_g[ids]*scale
        unchanged = np.setdiff1d(np.arange(len(c.gmax)), ids)
        with I.unbound_electrical(n), I.unbound_electrical(base):
            for command in (0., .5, 1., 0.):
                n.apply_weights(command);base.apply_weights(command)
                if not np.array_equal(c.gmax[ids], expected) or not np.array_equal(c.gmax[unchanged], b.gmax[unchanged]):
                    raise AssertionError('Stage5 contact routes were overwritten or changed other paths')
            n.apply_weights(1.);base.apply_weights(1.)
        event = float(cfg.get('cross_contact_event_amplitude', 0.))
        heel_max = event if event > 0 else 1.5*float(I.C.P.G['heel_rge'])
        toe_max = 1.5*float(I.C.P.G['toe_rge'])
        ranges = {};paths = []
        for side in ('r', 'l'):
            for kind, bound in (('HEEL', heel_max), ('TOE', toe_max)):
                cell = n.idx[kind+'_'+side];port = n.ports[kind+'_c_'+side]
                if np.any(c.c.g_max_non[cell]) or np.any(c.c.g_ion[:, cell]):
                    raise AssertionError('Contact encoder is not passive isolated input')
                gain = c.c.input_connectivity[cell, port]
                high = (c.c.g_m[cell]*c.c.V_rest[cell]+c.c.i_b[cell]+gain*bound)/c.c.g_m[cell]
                ranges[kind+'_'+side] = dict(current_upper_nA=bound, voltage_upper_mV=float(high),
                    input_gain=float(gain), tau_s=float(c.c.c_m[cell]/c.c.g_m[cell]))
            relay = n.idx['TOEDF_'+side];toe = n.idx['TOE_'+side]
            incoming = np.flatnonzero((c.post == relay) & (c.gmax != 0))
            if len(incoming) != 1 or c.pre[incoming[0]] != toe:
                raise AssertionError('Unexpected TOEDF driving path')
            edge = int(incoming[0]);u = ranges['TOE_'+side]['voltage_upper_mV']
            g = c.gmax[edge]*np.clip((u-c.elo[edge])/c.erange[edge], 0, 1)
            v = (c.c.g_m[relay]*c.c.V_rest[relay]+c.c.i_b[relay]+g*c.reversal[edge])/(c.c.g_m[relay]+g)
            ranges['TOEDF_'+side] = dict(isolated_TOE_only_voltage_upper_mV=float(v))
        for edge in ids:
            item = row(n, edge);pre, post = c.pre[edge], c.post[edge]
            if item['source'].startswith('HEEL_'):
                high = ranges[item['source']]['voltage_upper_mV']
            else:
                high = ranges[item['source']]['isolated_TOE_only_voltage_upper_mV']
            g = c.gmax[edge]*np.clip((high-c.elo[edge])/c.erange[edge], 0, 1)
            item['source_high_mV'] = high
            item['current_at_target_2p5_mV_nA'] = float(g*(c.reversal[edge]-2.5))
            item['active_next_hop'] = [row(n, i) for i in np.flatnonzero((c.pre == post) & (c.gmax != 0))]
            if item['source'].startswith('HEEL_'):
                v = (c.c.g_m[post]*c.c.V_rest[post]+c.c.i_b[post]+g*c.reversal[edge])/(c.c.g_m[post]+g)
                item['isolated_heel_only_PF_IN_mV'] = float(v)
                for next_edge in item['active_next_hop']:
                    a = np.clip((v-next_edge['low_mV'])/next_edge['span_mV'], 0, 1)
                    next_edge['isolated_heel_only_current_at_target_2p5_nA'] = float(next_edge['g_uS']*a*(next_edge['E_mV']-2.5))
            paths.append(item)
        protected = [row(n, i) for i, (pre, post) in enumerate(zip(c.pre, c.post))
                     if n.names[pre].startswith('HEEL_') and n.names[post].startswith(('InE_', 'InF_', 'RG_'))]
        results.append({'contact_pf_scale': scale, 'restored_contact_PF_pairs': paths,
                        'source_encoder_bounds': ranges, 'unchanged_heel_IN_and_direct_RG': protected,
                        'all_other_gains_exact_to_stage4_under_commands': True})
    result = {'source_condition': str(condition), 'physical_eval_sha256': I.sha(I.H/'physical_eval.py'),
              'neural_or_body_steps': 0, 'stage5_definition': 'Ground contact-PF refinement plus measured support-reduction ladder; not renamed perturbation.',
              'conditions': results,
              'range_scope': 'Upper endpoints from the actual caller and passive/TOEDF equations, not recorded maxima. Heel event is bounded0..amplitude; positive-upward regional toe-load interpretation assumed. Raw regional formula has no explicit lower clip, so no universal nonnegative toe-current claim.',
              'current_scope': 'Engineering endpoint currents at recipient2.5mV and isolated input-only relay equilibria; other active PF/central currents and transients can change actual delivered currents.',
              'functional_note': 'TOEDF inhibitory destination is PF_FUNC_ANKLE_Z_N (dorsiflexor channel), replacing the old shared PF_KNEE-F. Loaded toe suppresses DF, rather than exciting it. Heel excites ipsilateral extensor-labelled PF interneurons, which inhibit the corresponding flexor PF half-centers.'}
    Path(out).write_text(json.dumps(result, indent=2, sort_keys=True), encoding='utf-8')
    print(json.dumps({'out': str(out), 'pairs_per_scale': [len(r['restored_contact_PF_pairs']) for r in results]}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser();parser.add_argument('condition');parser.add_argument('out')
    args = parser.parse_args();run(args.condition, args.out)
