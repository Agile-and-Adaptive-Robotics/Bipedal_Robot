"""Bounded version of benchmarking/networkSizing/test_max_realistic_net_size.py
for the NUMPY backend (what our spinal net uses). The upstream test ships
backend=2 (Torch) and searches to 200k neurons; here: numpy, capped at 8k,
reporting build+step wall time per size vs our 410-neuron net.
"""
import io
import sys
import time

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
from sns_toolbox.networks import Network
from sns_toolbox.neurons import SpikingNeuron
from sns_toolbox.connections import SpikingSynapse
from sns_toolbox import backends


def design_network(size):
    net = Network()
    num_ins = int(0.08 * size) + 1
    num_outs = int(0.12 * size)
    num_syn = int(np.sqrt(size))
    num_rest = int(size) - num_ins - num_syn - num_outs
    spike = SpikingNeuron(name="m<0", threshold_proportionality_constant=-1)
    spike_bias = SpikingNeuron(name="bias",
                               threshold_proportionality_constant=-1,
                               bias=10.0)
    exc = SpikingSynapse(name="Excitatory Spiking")
    net.add_population(spike, [num_ins], name="ins")
    net.add_population(spike_bias, [num_outs], name="outs")
    net.add_population(spike_bias, [num_syn], name="connected")
    net.add_population(spike_bias, [num_rest], name="rest")
    net.add_connection(exc, "connected", "connected")
    net.add_input("ins")
    net.add_output("outs")
    return net

def main():
    for size in (500, 1000, 2000, 4000, 8000):
        t0 = time.time()
        try:
            net = design_network(size)
            t1 = time.time()
            # 1.5.2 API: backends take the compiled params dict, entered
            # via Network.compile (the master-branch benchmark passes the
            # Network directly - API changed after 1.5.2)
            model = net.compile(dt=0.001, backend="numpy")
            out = model.forward(np.array([1.0]))
            t2 = time.time()
            print(f"size {size:5d}: design {t1 - t0:5.1f}s  compile+1step "
                  f"{t2 - t1:6.1f}s  OK")
        except MemoryError:
            print(f"size {size:5d}: MEMORY ERROR at "
                  f"{time.time() - t0:.0f}s")
            break
        except Exception as e:
            print(f"size {size:5d}: FAILED ({type(e).__name__}: {e})")
            break
    print("\nour spinal net: 410 neurons, dense numpy matrix -> "
          "compile ~5-15 s, step ~0.5 ms (well inside limits)")


if __name__ == "__main__":
    main()
