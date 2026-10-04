"""GOAL 2 cross-platform check: SNS-Toolbox (numpy backend) vs the Simulink
SpikingNeuron + HybridSpikingSynapse + NonSpikingNeuron circuit.

Same circuit as sns_units_test_spiking.m TEST 2:
  10 nA -> SpikingNeuron A (Cm 5 nF, Gm 1 uS, Vrest 0, Vth 8, Vreset 0)
  A -> SpikingSynapse (gmax 0.3 uS, ginc 0.1 uS, tau_syn 100 ms, Esyn 8 mV)
    -> NonSpikingNeuron B (Cm 200 nF, Gm 1 uS, Vrest 0)

Toolbox mapping (README_SNS_Simscape + library header): C[uF] = Cm[nF]/1000,
reversal_potential relative = block Esyn - Vrest_post (0 here), fixed threshold
= m 0, threshold_increment 0. The Simulink run (units_ref_spiking_sim.mat) is
loaded and spike counts + ISI stats + V_B stats are compared.
"""
import numpy as np
from sns_toolbox.networks import Network
from sns_toolbox.neurons import SpikingNeuron, NonSpikingNeuron
from sns_toolbox.connections import SpikingSynapse

dt = 1e-6
t_end = 0.5

nrn_a = SpikingNeuron(membrane_capacitance=5e-3, membrane_conductance=1.0,
                      resting_potential=0.0, bias=0.0,
                      threshold_time_constant=5.0,
                      threshold_initial_value=8.0,
                      threshold_proportionality_constant=0.0,
                      threshold_increment=0.0,
                      reset_potential=0.0)
nrn_b = NonSpikingNeuron(membrane_capacitance=0.2, membrane_conductance=1.0,
                         resting_potential=0.0, bias=0.0)
syn = SpikingSynapse(max_conductance=0.3, reversal_potential=8.0,
                     time_constant=0.1, transmission_delay=0,
                     conductance_increment=0.1)

net = Network(name='goal2 cross check')
net.add_neuron(nrn_a, name='A', color='darkorange')
net.add_neuron(nrn_b, name='B', color='aqua')
net.add_input(dest='A', name='I0')
net.add_connection(syn, 'A', 'B')
net.add_output('A', name='VA', spiking=False)
net.add_output('A', name='SA', spiking=True)
net.add_output('B', name='VB')

model = net.compile(backend='numpy', dt=dt, debug=False)
n = int(round(t_end / dt))
data = np.zeros([n + 1, net.get_num_outputs_actual()])
data[0, :] = model(np.array([0.0]))  # settle t=0 with zero input? use 10 like sim
for i in range(n):
    data[i, :] = model(np.array([10.0]))
data[n, :] = model(np.array([10.0]))

t = np.arange(n + 1) * dt
va = data[:, 0]
sa = -data[:, 1]          # spike outputs come out as -spikes (0 / -1)
vb = data[:, 2]
spikes = sa < 0
tspk_tb = t[spikes & ~np.r_[False, spikes[:-1]]]
isi_tb = np.diff(tspk_tb)

print(f"toolbox numpy backend: {len(tspk_tb)} spikes in {t_end} s, "
      f"ISI {isi_tb.mean()*1e3:.4f} +- {isi_tb.std()*1e3:.4f} ms, "
      f"V_B end {vb[-1]:.5f} mV, V_B mean(last 0.2 s) {vb[t > 0.3].mean():.5f} mV")

from scipy.io import loadmat
m = loadmat(r'D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\results\units_ref_spiking_sim.mat',
            squeeze_me=True)
tspk_sim = np.atleast_1d(m['tspk'])
vt = np.asarray(m['vb_t']).ravel()
vd = np.asarray(m['vb_d']).ravel()
vb_sim_end = vd[-1]
vb_sim_mean = np.interp(t[t > 0.3], vt, vd).mean()
isi_sim = np.diff(tspk_sim)

ok_count = len(tspk_sim) == len(tspk_tb)
d_isi = abs(isi_sim.mean() - isi_tb.mean()) * 1e3
d_vb_end = abs(vb_sim_end - vb[-1])
d_vb_mean = abs(vb_sim_mean - vb[t > 0.3].mean())
print(f"Simulink:              {len(tspk_sim)} spikes, "
      f"ISI {isi_sim.mean()*1e3:.4f} +- {isi_sim.std()*1e3:.4f} ms, "
      f"V_B end {vb_sim_end:.5f} mV, V_B mean(last 0.2 s) {vb_sim_mean:.5f} mV")
print(f"diffs: spike count {'MATCH' if ok_count else 'MISMATCH'}, "
      f"|mean ISI| {d_isi:.5f} ms, |V_B end| {d_vb_end:.5f} mV, "
      f"|V_B mean(>0.3s)| {d_vb_mean:.5f} mV")
passed = ok_count and d_isi < 0.05 and d_vb_end < 0.05 and d_vb_mean < 0.02
print("CROSS-CHECK", "PASS" if passed else "FAIL")
