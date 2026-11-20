# AH001 articulated robot model

Open `AH001.xml` in MuJoCo 2.3.7. The file contains the robot geometry exported from CAD, articulated feet, tendon routes, and reference-tension actuators. It uses a 2 ms physics step. The positive y direction is superior.

Keep the `meshes` directory beside the XML. To load the model in Python, run the following from this directory.

```python
import mujoco

model = mujoco.MjModel.from_xml_path("AH001.xml")
data = mujoco.MjData(model)
mujoco.mj_forward(model, data)
```
