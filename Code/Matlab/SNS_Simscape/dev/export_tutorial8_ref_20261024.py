# Export the Tutorial-8 (persistent-Na half-center pair) numpy reference
# for the Simulink twin comparison. Mirrors tutorial8_run.py's cpg() with
# delta=-0.01; saves time + both traces to tutorial8_ref.mat.
import numpy as np
from scipy.io import savemat
from sns_toolbox.neurons import NonSpikingNeuronWithPersistentSodiumChannel
from sns_toolbox.connections import NonSpikingSynapse
from sns_toolbox.networks import Network


def zinf(U, Kz, Sz, Ez):
    return 1 / (1 + Kz * np.exp(Sz * (Ez - U)))


def run(delta=-0.01):
    Cm, Gm, Ena, Er, R = 5, 1, 50, -60, 20
    S = 0.05
    delEna = Ena - Er
    Km, Kh, Em, Eh = 1, 0.5, -40, -60
    delEm, delEh = Em - Er, Eh - Er
    tauHmax = 300
    Gna = Gm * R / (zinf(R, Km, S, delEm) * zinf(R, Kh, -S, delEh) * (delEna - R))

    neuron = NonSpikingNeuronWithPersistentSodiumChannel(
        membrane_capacitance=Cm, membrane_conductance=Gm,
        g_ion=[Gna], e_ion=[delEna],
        k_m=[Km], slope_m=[S], e_m=[delEm],
        k_h=[Kh], slope_h=[-S], e_h=[delEh], tau_max_h=[tauHmax],
        name='HC', color='orange')
    Ein = -100
    delEsyn = Ein - R
    gSyn = (-delta - delta * Gna * zinf(delta, Km, S, delEm)
            * zinf(delta, Kh, -S, delEh)
            + Gna * zinf(delta, Km, S, delEm) * zinf(delta, Kh, -S, delEh)
            * delEna) / (delta - delEsyn)
    syn = NonSpikingSynapse(max_conductance=gSyn, reversal_potential=delEsyn)

    net = Network()
    net.add_neuron(neuron, name='HC0', color='blue')
    net.add_input('HC0')
    net.add_output('HC0')
    net.add_neuron(neuron, name='HC1', color='orange')
    net.add_output('HC1')
    net.add_connection(syn, 'HC0', 'HC1')
    net.add_connection(syn, 'HC1', 'HC0')

    dt, tMax = 1, 5000
    t = np.arange(0, tMax, dt)
    Iapp = np.zeros(t.size)
    Ipert = np.zeros(t.size)
    Ipert[1] = 1
    model = net.compile(backend='numpy', dt=dt)
    data = np.zeros([t.size, net.get_num_outputs_actual()])
    for i in range(t.size):
        data[i] = model([Iapp[i] + Ipert[i]])
    return t, data.T, Gna, gSyn, delEsyn


t, traces, Gna, gSyn, delEsyn = run(-0.01)
print("Gna = %.6f uS, gSyn = %.6f uS, delEsyn = %.1f mV" % (Gna, gSyn, delEsyn))
print("HC0 V range [%.2f, %.2f] mV; final %.3f" %
      (traces[0].min(), traces[0].max(), traces[0][-1]))
savemat(r"D:\Github\Bipedal_Robot\Code\Matlab\SNS_Simscape\results\tutorial8_ref.mat",
        {"t": t, "HC0": traces[0], "HC1": traces[1],
         "Gna": Gna, "gSyn": gSyn, "delEsyn": delEsyn},
        do_compression=True)
print("saved results/tutorial8_ref.mat")
