# Decide empirically: are the Tutorial-8 pair's synapses dormant?
# Run (a) the full pair, (b) the same pair with gSyn = 0 (wired, zero gain),
# same kick; compare traces.
import numpy as np
from scipy.io import savemat
from sns_toolbox.neurons import NonSpikingNeuronWithPersistentSodiumChannel
from sns_toolbox.connections import NonSpikingSynapse
from sns_toolbox.networks import Network


def zinf(U, Kz, Sz, Ez):
    return 1 / (1 + Kz * np.exp(Sz * (Ez - U)))


def build(gsyn_scale=1.0):
    Er, R, S, Km, Kh, Em, Eh, Ena, Cm, Gm, tauHmax = (
        -60, 20, 0.05, 1, 0.5, -40, -60, 50, 5, 1, 300)
    delEna, delEm, delEh = Ena - Er, Em - Er, Eh - Er
    Gna = Gm * R / (zinf(R, Km, S, delEm) * zinf(R, Kh, -S, delEh) * (delEna - R))
    nrn = NonSpikingNeuronWithPersistentSodiumChannel(
        membrane_capacitance=Cm, membrane_conductance=Gm,
        g_ion=[Gna], e_ion=[delEna], k_m=[Km], slope_m=[S], e_m=[delEm],
        k_h=[Kh], slope_h=[-S], e_h=[delEh], tau_max_h=[tauHmax],
        name='HC', color='orange')
    delta = -0.01
    Ein = -100
    delEsyn = Ein - R
    gSyn = (-delta - delta * Gna * zinf(delta, Km, S, delEm)
            * zinf(delta, Kh, -S, delEh)
            + Gna * zinf(delta, Km, S, delEm) * zinf(delta, Kh, -S, delEh)
            * delEna) / (delta - delEsyn)
    syn = NonSpikingSynapse(max_conductance=gSyn * gsyn_scale,
                            reversal_potential=delEsyn)
    net = Network()
    net.add_neuron(nrn, name='HC0', color='blue')
    net.add_input('HC0')
    net.add_output('HC0')
    net.add_neuron(nrn, name='HC1', color='orange')
    net.add_output('HC1')
    net.add_connection(syn, 'HC0', 'HC1')
    net.add_connection(syn, 'HC1', 'HC0')
    return net


def run(net):
    dt, tMax = 1, 5000
    t = np.arange(0, tMax, dt)
    Iapp = np.zeros(t.size)
    Ipert = np.zeros(t.size)
    Ipert[1] = 1
    model = net.compile(backend='numpy', dt=dt)
    data = np.zeros([t.size, net.get_num_outputs_actual()])
    for i in range(t.size):
        data[i] = model([Iapp[i] + Ipert[i]])
    return t, data.T


t, full = run(build(1.0))
t2, zero = run(build(0.0))
d = np.max(np.abs(full - zero))
print("max |pair(gSyn) - pair(gSyn=0)| over 5 s: %.6g mV" % d)
print("full : HC0 t64=%.3f end=%.4f | HC1 end=%.4f" % (full[0][64], full[0][-1], full[1][-1]))
print("zero : HC0 t64=%.3f end=%.4f | HC1 end=%.4f" % (zero[0][64], zero[0][-1], zero[1][-1]))
savemat(r"D:\Github\Bipedal_Robot\Code\Matlab\SNS_Simscape\results\tutorial8_pair_probe.mat",
        {"t": t, "HC0_full": full[0], "HC1_full": full[1],
         "HC0_zero": zero[0], "HC1_zero": zero[1]})
