"""Persistent MuJoCo evaluator for coupled pendulum SNS fits."""
import json
from pathlib import Path

from serial_bench import evaluate, feedback_gain, model_xml, mujoco


MODEL = None
NEURAL = None
AXES = 1
PREDICTOR = False
OBSERVER = False
SHAPED = False


def handle_request(request):
    global MODEL, NEURAL, AXES, PREDICTOR, OBSERVER, SHAPED
    operation = request['op']
    if operation == 'preflight':
        fit = json.loads(Path(request['single_fit']).read_text())['parameters']
        NEURAL = {key: value for key, value in fit.items() if key not in ('kp', 'kd')}
        AXES = int(request.get('axes', 1))
        PREDICTOR = bool(request.get('model_predictor', False))
        OBSERVER = bool(request.get('disturbance_observer', False))
        SHAPED = bool(request.get('joint_shaping', False))
        MODEL = mujoco.MjModel.from_xml_string(model_xml(int(request['links']), AXES))
        return dict(ok=True, ready_for_optimization=True, dof=MODEL.nv)
    if MODEL is None or operation not in ('evaluate', 'held_out'):
        raise ValueError('Serial pendulum preflight is required')
    parameters = request['params']
    gain = feedback_gain(MODEL, parameters['effort'], neural_tau=parameters['neural_tau_s'],
                         distal_effort=parameters.get('distal_effort', 1.),
                         velocity_weight=parameters.get('velocity_weight', 1.))
    if operation == 'evaluate':
        cases = [(1, 8., 0., 1., 'last'), (1, 13., .014, 1.1, 'last')]
        if AXES == 2:
            cases = [(1, 8., 0., 1., 'other'), (1, 13., .014, 1.1, 'diagonal')]
        if PREDICTOR:
            cases[-1] = (1, 16., .020, 1.25, cases[-1][-1])
        if OBSERVER:
            cases[-1] = (1, 19., .030, 1.32, cases[-1][-1])
        if SHAPED:
            cases[-1] = (1, 20., .032, 1.30, cases[-1][-1])
    else:
        directions = ('last', 'other', 'diagonal') if AXES == 2 else ('last',) * 3
        challenges = (((11., .008, .9), (17., .024, 1.18), (19., .030, 1.32))
                      if PREDICTOR else ((10., .006, .85), (16., .020, 1.25), (18., .026, 1.35)))
        if OBSERVER:
            challenges = ((12., .012, .95), (18., .027, 1.21), (20., .032, 1.30))
        if SHAPED:
            challenges = ((14., .016, 1.02), (18.5, .028, 1.15), (21., .033, 1.28))
        cases = [(sign, push, delay, scale, axis) for sign in (1, -1)
                 for (push, delay, scale), axis in zip(challenges, directions)]
    rows = []
    for sign, push, delay, scale, axis in cases:
        row, _ = evaluate(MODEL, gain, 'SNS lag-aware', sign, push, NEURAL, delay,
                          parameters['delay_prediction'], scale, axis, PREDICTOR,
                          parameters.get('disturbance_gain', 0.))
        rows.append(row)
    costs = [1000. * (not row['passed']) + 30. * (6. - row['simulated_s'])
             + .3 * row['peak_angle_deg'] + 2. * (row['final_rms_deg'] or 0.)
             + (row['settling_s'] if row['settling_s'] is not None else 4.)
             for row in rows]
    return dict(ok=True, value=float(sum(costs) / len(costs) + .25 * max(costs)),
                passed=sum(row['passed'] for row in rows), total=len(rows),
                conditions=rows, gain=gain.tolist())
