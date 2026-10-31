from pathlib import Path
from lxml import etree as E
import json
p=Path('Neuromechanical_Models/W2L_Treadmill_FrozenRoot')
for gui in [False,True]:
 r=E.parse(str(p/('W2L_Treadmill_FrozenRoot.aproj' if gui else 'W2L_Treadmill_FrozenRoot_Standalone.asim')))
 for b in r.findall('.//RigidBody'):
  if b.findtext('Name') not in ['Root','femur_L','tibia_L','foot_L','femur_R','tibia_R','foot_R']:continue
  print('GUI',gui,b.findtext('Name'))
  print(E.tostring(b,pretty_print=True).decode().split('<ChildBodies>')[0][-9000:])
