"""Export separate PF motor pools and central sensory routes from the fitted graph."""
import argparse
import hashlib
from html import escape
import json
from pathlib import Path

import motor_adapter as motor
from controller_profiles_20261104.native_pf import motor_digest, motor_rows


def drawing(outgoing, feedback):
    """Nested drawings use a separate muscle page so dense wiring stays readable."""
    def node(kind, label, x, y, sub=None):
        result = dict(type=kind, label=label, x=x, y=y)
        if sub is not None:
            result['sub'] = sub
        return result

    def edge(source, target, tag, gain=1., sign='exc', points=()):
        return {'from': source, 'to': target, 'sign': sign, 'gain': gain, 'tag': tag,
                'pts': [dict(x=x, y=y) for x, y in points]}

    def spec(nodes, edges=()):
        return dict(nodes=nodes, edges=list(edges))

    def parameters(g, lo, hi, rev):
        return f'gmax {g:.6g} µS; Erev {rev:g} mV; activation {lo:g}–{hi:g} mV'

    central = {(row['source'], row['target']): row for row in feedback}
    by_source = {}
    for row in outgoing:
        if row[2] > 0 or row[:2] in central:
            by_source.setdefault(row[0], []).append(row)

    def destination_page(source, routes, expand_relays=True):
        nodes = [node('PORT-load', source, 110, 120)]
        edges = []
        for i, (_, target, g, lo, hi, rev) in enumerate(routes):
            y = 95 + 92 * i
            dest = target + ' → off-page'
            relay_routes = (by_source.get(target, ()) if expand_relays and
                            target.startswith('IaIN_') else ())
            nodes.append(node('SUB' if relay_routes else 'PORT-load', dest, 570, y,
                              destination_page(target, relay_routes, False) if relay_routes else None))
            route = central.get((source, target))
            tag = parameters(route['gmax_uS'] if route else g, lo, hi, rev)
            if route:
                tag += '; presynaptic modulation' + (' and load gate' if route['stance_gated'] else '')
            edges.append(edge(source, dest, tag, route['gmax_uS'] if route else g,
                              'inh' if rev < 0 else 'exc'))
        return spec(nodes, edges)

    def knee_comparison():
        nodes, edges = [], []
        for column, (side, channel, muscle_name) in enumerate(
                (('l', 2, 'rect_fem'), ('l', 3, 'bifemsh'),
                 ('r', 2, 'rect_fem'), ('r', 3, 'bifemsh'))):
            x = 245 + 640 * column + (100 if column >= 2 else 0)
            leg = 'Left' if side == 'l' else 'Right'
            act = f'{muscle_name}_{side}'
            pf = f'PF_S{channel}_{side}'
            motor = 'MN_' + act
            _, _, gain, lo, hi, rev = next(row for row in by_source[pf] if row[1] == motor)
            title = f'{leg} S{channel} · ' + ('knee extension' if channel == 2 else 'knee flexion')
            adapter = f'{act} · MN voltage to excitation'
            muscle = f'{act} · muscle'
            spindle = f'{act} · spindle'
            tendon = f'{act} · tendon organ'
            nodes.extend((node('PORT-load', title, x, 75),
                          node('MN', motor, x, 175),
                          node('PORT-load', adapter, x, 260),
                          node('MUSCLE', muscle, x, 345),
                          node('PORT-load', spindle, x - 150, 445),
                          node('PORT-load', tendon, x + 150, 445)))
            edges.extend((edge(title, motor, parameters(gain, lo, hi, rev), gain),
                          edge(motor, adapter, 'Graded MN voltage'),
                          edge(adapter, muscle, 'Normalized excitation 0–1'),
                          edge(muscle, spindle, 'Actuator length and velocity'),
                          edge(muscle, tendon, 'Positive tensile force')))
            for i, (kind, receptor, signal) in enumerate(
                    (('Ia', spindle, 'velocity'), ('II', spindle, 'length'),
                     ('Ib', tendon, 'tensile force'))):
                ax = x + 200 * (i - 1)
                encoder = f'{kind} {signal} adapter · {side}{channel}'
                aff = f'{kind}_{act}'
                routes = by_source.get(aff, ())
                central_targets = [r[1].replace('PF_', '').replace('RG_', 'RG-') for r in routes
                                   if (r[0], r[1]) in central]
                suffix = ', '.join(central_targets) if central_targets else 'local reflex'
                dest = f'{kind} off-page · {suffix} · {side}{channel}'
                nodes.extend((node('PORT-load', encoder, ax, 545),
                              node('SN-' + kind, aff, ax, 635),
                              node('SUB', dest, ax, 730, destination_page(aff, routes))))
                edges.extend((edge(receptor, encoder, signal),
                              edge(encoder, aff, 'Afferent input current'),
                              edge(aff, dest, 'Open for all exact central and local destinations')))
            aff = 'Ia_' + act
            relay = 'IaIN_' + act
            iain_input = next((row for row in by_source.get(aff, ()) if row[1] == relay), None)
            if iain_input:
                _, _, g, lo, hi, rev = iain_input
                targets = by_source.get(relay, ())
                off_page = f'{relay} → MN and IaIN targets'
                nodes.extend((node('IN-IaIN', relay, x - 200, 850),
                              node('SUB', off_page, x - 200, 945, destination_page(relay, targets, False))))
                edges.extend((edge(aff, relay, parameters(g, lo, hi, rev), g,
                                   points=((x - 285, 650), (x - 285, 805))),
                              edge(relay, off_page, 'Open for exact inhibitory targets', sign='inh')))
        return spec(nodes, edges)

    def muscle_page(pf, motor, g, lo, hi, rev):
        act = motor[3:]
        pf_port = pf + ' ← off-page'
        adapter = 'Excitation adapter · clip(V/EHI)'
        muscle = act + ' · MuJoCo muscle'
        spindle = 'Spindle · length and velocity'
        gto = 'Tendon organ · tensile force'
        nodes = [node('PORT-load', pf_port, 110, 95), node('MN', motor, 355, 95),
                 node('PORT-load', adapter, 605, 95), node('MUSCLE', muscle, 875, 95),
                 node('PORT-load', spindle, 250, 250), node('PORT-load', gto, 250, 580)]
        edges = [edge(pf_port, motor, parameters(g, lo, hi, rev), g),
                 edge(motor, adapter, 'Graded MN voltage'), edge(adapter, muscle, 'Muscle excitation 0–1'),
                 edge(muscle, spindle, 'Actuator length and velocity', points=((875, 174), (250, 174))),
                 edge(muscle, gto, 'Positive tensile force', points=((1010, 95), (1010, 516), (250, 516)))]
        for i, (kind, sensor, text) in enumerate(
                (('Ia', spindle, 'Velocity adapter · v / 0.6'),
                 ('II', spindle, 'Length adapter · normalized length'),
                 ('Ib', gto, 'Force adapter · tension / Fmax'))):
            y = (260, 415, 590)[i]
            aff = f'{kind}_{act}'
            nodes.extend((node('PORT-load', text, 540, y), node('SN-' + kind, aff, 820, y)))
            edges.extend((edge(sensor, text, 'Measured muscle signal'), edge(text, aff, 'Afferent input current')))
            routes = by_source.get(aff, ())
            groups = [('PF', [r for r in routes if r[1].startswith('PF_')]),
                      ('RG', [r for r in routes if r[1].startswith('RG_')]),
                      ('Local reflex', [r for r in routes if not r[1].startswith(('PF_', 'RG_'))])]
            for j, (layer, selected) in enumerate(groups):
                if selected:
                    label = f'{kind} → {layer} off-page ({len(selected)})'
                    nodes.append(node('SUB', label, 1110, y + 52 * j, destination_page(aff, selected)))
                    edges.append(edge(aff, label, 'Open for exact destinations and synapses'))
        exits = [row for row in by_source.get(motor, ()) if not row[1].startswith('MN_')]
        if exits:
            label = 'MN recurrent exits off-page'
            nodes.append(node('SUB', label, 355, 25, destination_page(motor, exits)))
            edges.append(edge(motor, label, 'Open for recurrent destinations'))
        return spec(nodes, edges)

    legs = []
    for side, leg, x in (('l', 'Left leg', 250), ('r', 'Right leg', 700)):
        channels = []
        for channel in range(1, 7):
            source = f'PF_S{channel}_{side}'
            motors = sorted((row for row in by_source[source] if row[1].startswith('MN_')),
                            key=lambda row: (-row[2], row[1]))
            members = [node('SUB', row[1][3:], 125 + 210 * (i % 5), 130 + 95 * (i // 5),
                            muscle_page(*row)) for i, row in enumerate(motors)]
            channels.append(node('SUB', f'S{channel} · {len(motors)} muscles',
                                 185 + 235 * ((channel - 1) % 3), 140 + 150 * ((channel - 1) // 3),
                                 spec(members)))
        legs.append(node('SUB', leg, x, 220, spec(channels)))
    comparison = knee_comparison()
    legs.append(node('SUB', 'Bilateral S2 / S3 knee comparison', 475, 80, comparison))
    return spec(legs) | {'_note': 'Editable fitted syn6 anatomy. Double-click the comparison or any subsystem.'}


def table(headers, rows):
    def cells(tag, values):
        return '<tr>' + ''.join(f'<{tag}>{escape(str(v))}</{tag}>' for v in values) + '</tr>'
    return ('<table>' + cells('th', headers) + ''.join(cells('td', row) for row in rows) + '</table>')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fit', type=Path, required=True)
    parser.add_argument('--body', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--connectome', type=Path, help='Write an editable nested Connectome Block Editor diagram')
    parser.add_argument('--editor', type=Path, help='Write a standalone editor with this diagram preloaded')
    args = parser.parse_args()
    parent = Path(json.loads((args.fit / 'MANIFEST.json').read_text())['fit'])
    raw = (parent / 'FINAL.json').read_bytes()
    preflight = motor.handle_request(dict(op='preflight', fit=str(parent),
                                          fit_sha256=hashlib.sha256(raw).hexdigest()))
    if not preflight['ready_for_optimization']:
        raise RuntimeError('Frozen readout reconstruction failed')
    replay = motor.INSTANCE
    fit = {**motor.FIXED, **json.loads((args.fit / 'FINAL.json').read_text())['parameters']}
    replay.evaluate(fit)
    compiled, names = replay.c, replay.names
    digest = motor_digest(compiled, names)
    if (args.body / 'MOTOR_SHA256.txt').read_text().strip() != digest:
        raise ValueError('Recorded run and export MN conductances differ')
    blocks = []
    outgoing = [(names[pre], names[post], float(g), float(lo), float(lo + span), float(rev))
                for pre, post, g, lo, span, rev in zip(compiled.pre, compiled.post, compiled.gmax,
                                                       compiled.elo, compiled.erange, compiled.reversal)]
    motors = motor_rows(compiled, names)
    for channel in range(1, 7):
        label = f'S{channel}'
        panels = []
        for side in ('r', 'l'):
            source = f'PF_{label}_{side}'
            rows = sorted((row for row in motors if row[0] == source and row[2] > 0),
                          key=lambda row: -row[2])
            panels.append(f'<section><h3>{source} → {len(rows)} motor pools</h3>' +
                          table(('MN', 'Muscle command destination', 'gmax (µS)', 'Presynaptic range (mV)', 'Erev (mV)'),
                                ((mn, mn[3:], f'{g:.5g}', f'{lo:g} to {lo + span:g}', f'{rev:g}')
                                 for _, mn, g, lo, span, rev in rows)) + '</section>')
        other = [(src, dst, f'{g:.5g}', f'{rev:g}') for src, dst, g, _, _, rev in outgoing
                 if src in (f'PF_{label}_r', f'PF_{label}_l') and not dst.startswith('MN_') and g > 0.]
        blocks.append(f'<details><summary>{label} motor-neuron projections</summary>'
                      + '<div class="panels">' + ''.join(panels) + '</div>'
                      + '<h3>Other compiled exits</h3>'
                      + '<p>Legacy PF → PF interneuron exits fade to zero during native handover. '
                        'Ia interneuron exits remain.</p>'
                      + table(('Source', 'Target', 'gmax (µS)', 'Erev (mV)'), other) + '</details>')
    feedback = [row for row in replay.net.central_revision['added'] if row['purpose'] == 'afferent']
    if args.connectome or args.editor:
        diagram = drawing(outgoing, feedback)
        payload = json.dumps(diagram, ensure_ascii=False, separators=(',', ':'))
        if args.connectome:
            args.connectome.write_text(payload, encoding='utf-8')
        if args.editor:
            source = Path(__file__).resolve().parents[1] / 'connectome_block_editor.html'
            editor = source.read_text(encoding='utf-8')
            blank = 'let models = [{ name: "Ben", data: { nodes: [], edges: [] } }];'
            if editor.count(blank) != 1:
                raise ValueError('Connectome editor boot structure changed')
            comparison = diagram['nodes'][-1]['sub']
            model = json.dumps([dict(name='S2 and S3 knee comparison', data=comparison),
                                dict(name='All PF and muscle subnetworks', data=diagram)],
                               ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
            editor = editor.replace(blank, f'let models = {model};')
            editor = editor.replace('"cme_models_v1"', '"cme_syn6_knee_comparison_20261105"')
            editor = editor.replace('loadIntoCanvas(cur(), true);\nstatus("ready',
                                    'loadIntoCanvas(cur(), true);\nfitView();\nstatus("ready')
            args.editor.write_text(editor, encoding='utf-8')
    for layer in ('PF', 'RG'):
        rows = [row for row in feedback if row['target_layer'] == layer]
        blocks.append(f'<details><summary>{layer} feedback · {len(rows)} bilateral routes</summary>'
                      + table(('Source', 'Target', 'Base gmax (µS)', 'Fit multiplier', 'Load gated'),
                              ((row['source'], row['target'], f"{row['gmax_uS']:.5g}",
                                f"{fit.get(row['afferent']+'_'+layer, 1.):.4g}",
                                'Yes' if row['stance_gated'] else 'No') for row in rows))
                      + '</details>')
    page = """<!doctype html><html lang="en"><meta charset="utf-8"><title>Syn6 wiring</title>
<style>body{font:16px system-ui;margin:2rem auto;max-width:1250px;padding:0 1rem;color:#182539}
details{border:1px solid #bcc8d7;border-radius:6px;margin:.8rem 0;padding:.7rem}
summary{cursor:pointer;font-weight:650}table{border-collapse:collapse;width:100%;font-size:.85rem}
td,th{text-align:left;padding:.3rem;border-bottom:1px solid #dbe2eb}
.panels{display:grid;grid-template-columns:1fr 1fr;gap:1rem}code{overflow-wrap:anywhere}
@media(max-width:850px){.panels{grid-template-columns:1fr}}</style>
<h1>Syn6 PF, motor pools, and central feedback</h1>
<p>Each PF cell projects a graded excitatory conductance to the listed motor pools.
An MN integrates PF, local reflex, and recurrent inputs. Its output is the normalized
muscle excitation <code>clip(V<sub>MN</sub>/E<sub>HI</sub>, 0, 1)</code> sent to the named MuJoCo actuator.
Multiple PF cells can feed the same MN. No force target exits an MN.</p>
<h2>Live Stage 2 feedback</h2>
<p>Muscle velocity → Ia · normalized length → II · tensile force → Ib.
Central chemical edges below are excitatory. The descending flexor command attenuates
presynaptic gain, and the Ib gain also follows the ipsilateral load proxy.
Heel routes are enabled in grounded trials, where heel events excite ipsilateral
InE and contralateral RG-F, with an ipsilateral PF inhibitory relay branch where
present. The Stage 2 air run disables heel events.</p>
<p>Knee angle and extension velocity → knee premotor relay. In the unrestricted
native experiment it excites S3, suppresses the S3 quiet shunt, and inhibits S2.
The single-window experiment leaves the S2 inhibition live and gives S3 its own
RG-timed relay after handover. This experimental path has no direct RG or MN exit.</p>
"""
    args.output.write_text(page + ''.join(blocks) + f'<p>MN SHA-256 <code>{digest}</code></p></html>',
                           encoding='utf-8')


if __name__ == '__main__':
    main()
