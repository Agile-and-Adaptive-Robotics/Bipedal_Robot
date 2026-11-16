"""Replay frozen serial pendulum trials and export their figures and wiring."""
import argparse
import json
from pathlib import Path

from serial_bench import CAPACITY, build, evaluate, model_xml, mujoco, np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont


CASES = (
    ('double', 2, 1, 'observer_rustuna_2link_20261106', 'last'),
    ('triple', 3, 1, 'shaped_rustuna_3link_20261106', 'last'),
    ('spatial4', 2, 2, 'spatial_observer_2link_20261106', 'diagonal'),
    ('spatial6', 3, 2, 'spatial_observer_3link_20261106', 'diagonal'),
)
COLORS = ('#0072B2', '#D55E00', '#009E73')


def plots(out, tag, model, log, capture, result):
    t = log[:, 0] - 1.5
    n = model.nv
    fig, panels = plt.subplots(4, n, figsize=(max(8, 2.55 * n), 7.8),
                               sharex=True, sharey='row', squeeze=False)
    for j in range(n):
        name = model.joint(j).name
        panels[0, j].set_title(f'Link {int(name[5]) + 1} · ' + ('pitch' if name.endswith('z') else 'roll'))
        for polarity, (color, style) in enumerate(zip(COLORS[:2], ('-', '--'))):
            label = ('Positive', 'Negative')[polarity]
            for row, key in ((0, 'premotor_mV'), (1, 'motor_mV'), (2, 'actuator_activation')):
                value = capture[key][:, j, polarity]
                if row == 2:
                    value = CAPACITY * value
                panels[row, j].plot(t, value, color=color, ls=style, lw=1.25, label=label)
        panels[3, j].plot(t, np.degrees(log[:, j + 1]), color='#222222', lw=1.35)
        panels[3, j].axhline(0, color='.6', lw=.6)
        for ax in panels[:, j]:
            ax.axvspan(0, .05, color='#F0E442', alpha=.4, zorder=0)
            ax.set_xlim(-.2, 1.35)
            ax.grid(alpha=.16)
        panels[3, j].set_xlabel('Time from push (s)')
    panels[0, 0].legend(loc='upper right', fontsize=8)
    for row, label in enumerate(('Premotor (mV)', 'Motor neuron (mV)',
                                 'Actuator torque\nmagnitude (N·m)', 'Joint angle (°)')):
        panels[row, 0].set_ylabel(label)
    axis = 'diagonal' if result['push_axis'] == 'diagonal' else 'distal pitch'
    fig.suptitle(f"{n} DOF · 20 N·m {axis} push · 32 ms task sensing delay", fontsize=14)
    fig.tight_layout()
    fig.savefig(out / f'{tag}_responses.png', dpi=155)
    plt.close(fig)


def animation(out, tag, model, log, result):
    data = mujoco.MjData(model)
    lengths = (.45, .42, .35)[:result['links']]
    # Orthographic view that exposes both horizontal directions in spatial cases.
    horizontal = np.array([.80, 0., -.60]) if result['dof'] > result['links'] else np.array([1., 0., 0.])
    scale = 300 / sum(lengths)
    font = ImageFont.truetype('arial.ttf', 20)
    small = ImageFont.truetype('arial.ttf', 15)
    frames = []
    for row in log[(log[:, 0] >= 1.1) & (log[:, 0] <= 3.3)][::6]:
        data.qpos[:] = row[1:1 + model.nq]
        mujoco.mj_forward(model, data)
        canvas = Image.new('RGB', (510, 435), '#f8fafc')
        draw = ImageDraw.Draw(canvas)
        origin = np.array([255., 375.])
        def xy(point):
            return tuple(origin + scale * np.array([point @ horizontal, -point[1]]))
        draw.line((255, 73, 255, 375), fill='#cbd5e1', width=2)
        draw.line((98, 375, 412, 375), fill='#334155', width=4)
        for j, length in enumerate(lengths):
            body = data.body(j + 1)
            start = xy(body.xpos)
            tip = xy(body.xpos + length * body.xmat.reshape(3, 3)[:, 1])
            draw.line([start, tip], fill=COLORS[j], width=14)
            x, y = start
            draw.ellipse((x - 8, y - 8, x + 8, y + 8), fill='#142536', outline='white', width=2)
        elapsed = row[0] - 1.5
        phase = ('PUSH' if 0 <= elapsed < .05 else
                 'RECOVERY' if elapsed >= .05 else 'BEFORE PUSH')
        draw.text((18, 14), f"{result['links']} links · {result['dof']} DOF · {phase}", fill='#142536', font=font)
        draw.text((18, 402), f"t = {elapsed:+.3f} s  |  20 N·m for 50 ms  |  32 ms task delay",
                  fill='#334155', font=small)
        frames.append(canvas)
    frames[0].save(out / f'{tag}_recovery.gif', save_all=True, append_images=frames[1:],
                   duration=30, loop=0, optimize=True)


def connectome(out, neural):
    def node(kind, label, x, y, sub=None):
        result = dict(type=kind, label=label, x=x, y=y)
        if sub is not None:
            result['sub'] = sub
        return result

    def edge(source, target, tag, gain=1., sign='exc'):
        return {'from': source, 'to': target, 'tag': tag, 'gain': gain, 'sign': sign}

    def axis_page():
        nodes, edges = [], []
        for x, side in ((420, 'POS'), (710, 'NEG')):
            inward = -1 if side == 'POS' else 1
            labels = [('PORT-load', 'Task command ' + side, x, 65),
                      ('HC-PF-E', 'PP_' + side, x, 170),
                      ('IN-PF', 'PP_IN_' + side, x + 80 * inward, 280),
                      ('MN', 'MN_' + side, x, 400),
                      ('PORT-load', 'Excitation adapter ' + side, x, 500),
                      ('MUSCLE', 'Filtered torque actuator ' + side, x, 600),
                      ('PORT-load', 'Virtual length and velocity ' + side, x + 275 * inward, 420),
                      ('SN-Ia', 'Ia_' + side, x + 270 * inward, 280),
                      ('SN-II', 'II_' + side, x + 270 * inward, 170),
                      ('SN-Ib', 'Ib_' + side, x + 270 * inward, 620),
                      ('IN-IbIN', 'IbRelay_' + side, x + 270 * inward, 520)]
            nodes.extend(node(*row) for row in labels)
            edges.extend((
                edge('Task command ' + side, 'PP_' + side, 'Sign split and inverse steady-state encoding'),
                edge('MN_' + side, 'Excitation adapter ' + side, 'Graded MN voltage'),
                edge('Excitation adapter ' + side, 'Filtered torque actuator ' + side,
                     '0.05 + 0.95 clip(V_MN / 5, 0, 1)'),
                edge('Filtered torque actuator ' + side, 'Ib_' + side, 'Activation × 1.5'),
                edge('Virtual length and velocity ' + side, 'Ia_' + side, 'Lengthening velocity proxy'),
                edge('Virtual length and velocity ' + side, 'II_' + side, 'Length proxy'),
            ))
        for row in build(parameters=neural).edges:
            g, rev = row['gmax_uS'], row['reversal_mV']
            edges.append(edge(row['source'], row['target'], f'{g:.6g} µS · Erev {rev:g} mV · 0–5 mV',
                              g, 'inh' if rev < 0 else 'exc'))
        return dict(nodes=nodes, edges=edges)

    tabs = []
    for tag, links, axes, _, _ in CASES:
        nodes = [node('PORT-load', 'Delayed joint angle and velocity', 225, 75),
                 node('PORT-load', 'Torque history and disturbance observer', 575, 75),
                 node('PORT-load', 'Coupled lag-aware task controller', 400, 180),
                 node('MUSCLE', 'MuJoCo pendulum and actuator states', 400, 590)]
        edges = [edge(nodes[i]['label'], nodes[2]['label'], 'State estimate') for i in (0, 1)]
        edges.extend((edge(nodes[3]['label'], nodes[0]['label'], '32 ms task sensor delay'),
                      edge(nodes[3]['label'], nodes[1]['label'], 'Actuator activation history')))
        for j in range(links * axes):
            name = f'Link {j // axes + 1} · ' + ('pitch' if j % axes == 0 else 'roll')
            nodes.append(node('SUB', name, 100 + 300 * (j % 3), 310 + 135 * (j // 3), axis_page()))
            edges.extend((edge(nodes[2]['label'], name, 'Signed axis command'),
                          edge(name, nodes[3]['label'], 'Two actuator excitations'),
                          edge(nodes[3]['label'], name, 'Local Ia, II, and Ib proxies')))
        tabs.append(dict(name=f'{links * axes} DOF · {tag}', data=dict(nodes=nodes, edges=edges)))
    (out / 'pendulum_connectome.json').write_text(json.dumps(tabs[-1]['data'], ensure_ascii=False), encoding='utf-8')
    editor_path = Path(__file__).resolve().parents[1] / 'connectome_block_editor.html'
    editor = editor_path.read_text(encoding='utf-8')
    blank = 'let models = [{ name: "Ben", data: { nodes: [], edges: [] } }];'
    if editor.count(blank) != 1:
        raise ValueError('Connectome editor boot structure changed')
    editor = editor.replace(blank, 'let models = ' + json.dumps(tabs, ensure_ascii=False) + ';')
    editor = editor.replace('"cme_models_v1"', '"cme_serial_pendulum_recovery"')
    editor = editor.replace('loadIntoCanvas(cur(), true);\nstatus("ready',
                            'loadIntoCanvas(cur(), true);\nfitView();\nstatus("ready')
    (out / 'pendulum_connectome_editor.html').write_text(editor, encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    root = args.record_root
    fit = json.loads((root / 'robust_delay_fit_v4/FINAL.json').read_text())['parameters']
    neural = {k: v for k, v in fit.items() if k not in ('kp', 'kd')}
    matched = json.loads((root / 'MATCHED_SERIAL_BENCHMARK_20261106.json').read_text())
    results = []
    for tag, links, axes, folder, push_axis in CASES:
        final = json.loads((root / folder / 'FINAL.json').read_text())
        best = final['best']
        expected = matched if axes == 1 else final['holdout']['conditions']
        reference = next(row for row in expected if row['push_Nm'] == 20. and
                         row['mirrored_sign'] == 1 and row.get('checkpoint', folder) == folder)
        model = mujoco.MjModel.from_xml_string(model_xml(links, axes))
        capture = {key: np.empty((1200, model.nv, 2)) for key in
                   ('premotor_mV', 'motor_mV', 'actuator_command', 'actuator_activation')}
        capture['push_Nm'] = np.empty((1200, model.nv))
        params = best['parameters']
        result, log = evaluate(model, np.asarray(best['gain']), 'SNS lag-aware', 1, 20., neural, .032,
                               params['delay_prediction'], 1.30, push_axis, True,
                               params['disturbance_gain'], capture=capture)
        if not result['passed'] or any(abs(result[key] - reference[key]) > 1e-9 for key in
                                      ('peak_angle_deg', 'final_rms_deg', 'settling_s')):
            raise RuntimeError(f'{tag} replay differs from its recorded checkpoint')
        capture = {key: val[:len(log)] for key, val in capture.items()}
        np.savez_compressed(args.output / f'{tag}_trace.npz', time_state=log, **capture)
        plots(args.output, tag, model, log, capture, result)
        animation(args.output, tag, model, log, result)
        results.append(dict(tag=tag, checkpoint=str(root / folder / 'FINAL.json'), **result))
        (args.output / 'REPLAY_RESULTS.json').write_text(json.dumps(results, indent=2))
        print(tag, result['settling_s'], flush=True)
    connectome(args.output, neural)


if __name__ == '__main__':
    main()
