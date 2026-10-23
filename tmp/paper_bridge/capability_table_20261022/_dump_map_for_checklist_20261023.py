import json
m = json.load(open(r'D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Routing27_v2\actuator_map_27.json'))
print('rom_table_deg:')
for k, v in m['rom_table_deg'].items():
    print(' ', k, v)
print('actuator count:', len(m['actuators']))
for a in m['actuators']:
    print(a['id'], a['master_name'], '| primary:', a['primary_dof'],
          '| secondaries:', a['secondary_dofs'])
