# Empirically pin Tutorial-8 semantics: single HC neuron no-input rest,
# and the pair trajectory samples from the exported reference.
import numpy as np
from scipy.io import loadmat, savemat
from sns_toolbox.neurons import NonSpikingNeuronWithPersistentSodiumChannel
from sns_toolbox.networks import Network


def zinf(U, Kz, Sz, Ez):
    return 1 / (1 + Kz * np.exp(Sz * (Ez - U)))


Er, R, S, Km, Kh, Em, Eh, Ena, Cm, Gm, tauHmax = (
    -60, 20, 0.05, 1, 0.5, -40, -60, 50, 5, 1, 300)
delEna, delEm, delEh = Ena - Er, Em - Er, Eh - Er
Gna = Gm * R / (zinf(R, Km, S, delEm) * zinf(R, Kh, -S, delEh) * (delEna - R))
nrn = NonSpikingNeuronWithPersistentSodiumChannel(
    membrane_capacitance=Cm, membrane_conductance=Gm,
    g_ion=[Gna], e_ion=[delEna], k_m=[Km], slope_m=[S], e_m=[delEm],
    k_h=[Kh], slope_h=[-S], e_h=[delEh], tau_max_h=[tauHmax],
    name='HC', color='orange')

net = Network()
net.add_neuron(nrn, name='HC', color='blue')
net.add_input('HC')
net.add_output('HC')
model = net.compile(backend='numpy', dt=1)
u = 0.0
for i in range(5000):
    u = model([0.0])[0]
print("single HC, no input, from U=0: settles at U = %.4f" % u)
for u0 in (5.0, 20.0, 40.0):
    model2 = net.compile(backend='numpy', dt=1)
    uu = u0
    for i in range(3000):
        uu = model2([0.0])[0]
    print("  from U0=%.0f -> %.4f" % (u0, uu))

ref = loadmat(r"D:\Github\Bipedal_Robot\Code\Matlab\SNS_Simscape\results\tutorial8_ref.mat")
t = ref['t'].ravel(); h0 = ref['HC0'].ravel()
print("pair HC0: t=0: %.3f, t=100: %.3f, t=500: %.3f, t=1000: %.3f, "
      "t=2000: %.3f, t=4999: %.3f" % (
          h0[0], h0[100], h0[500], h0[1000], h0[2000], h0[-1]))
ipk = [i for i in range(1, 4999) if h0[i] > h0[i - 1] and h0[i] > h0[i + 1]
       and h0[i] > 15]
print("peaks>15mV at t =", [t[i] for i in ipk])
savemat(r"D:\Github\Bipedal_Robot\Code\Matlab\SNS_Simscape\results\tutorial8_probe.mat",
        {"u_rest_single": u})
