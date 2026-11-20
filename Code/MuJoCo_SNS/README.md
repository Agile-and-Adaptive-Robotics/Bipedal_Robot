# Neuromechanical simulations

Install Python 3.10, MuJoCo 2.3.7, SNS-Toolbox 1.5.2, NumPy, SciPy, and Matplotlib in the same environment. Run the commands below from this directory.

| Command | Example |
| --- | --- |
| `python sns_bpa_demo.py` | Applies the pneumatic actuator force model to a simulated knee. |
| `python spinal/check_rhythm.py` | Tests the spinal circuit's rhythm in isolation. |
| `python spinal/runner.py --time 1` | Couples the spinal circuit to the included gait2392 human model for a one-second trial. |
| `python spinal/w2l_cpg/smoke_w2l.py` | Tests bilateral alternation in the walker circuit with synthetic contact input. |

The [desktop connectome editor](connectome_app/README.md) has its own dependencies. The [CAD-based robot](../../Neuromechanical_Models/AH001_MuJoCo/) also opens directly in MuJoCo. The default human and robot OpenSim models are in [Gait2392_Robotbody](../../Solid_Models/OpenSim/Gait2392_Robotbody/).
