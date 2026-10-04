"""List stock gait2392 right-side muscle MIFs to ground the 27-actuator map.

Read-only. Prints muscle name and max_isometric_force for every right-side
Thelen2003Muscle in the stock model, plus candidate group sums used by the
routing-campaign map decisions.
"""
import xml.etree.ElementTree as ET

MODEL = (r'D:\Github\Bipedal_Robot\Solid_Models\OpenSim'
         r'\Gait2392_Robotbody\gait2392_simbody.osim')

tree = ET.parse(MODEL)
root = tree.getroot()

mif = {}
for mu in root.iter('Thelen2003Muscle'):
    name = mu.get('name')
    f = mu.find('max_isometric_force')
    if name and f is not None and name.endswith('_r'):
        mif[name] = float(f.text)

for name in sorted(mif):
    print(f'{name:14s} {mif[name]:8.1f}')

print()
groups = {
    'glut_max1+2+3': ['glut_max1_r', 'glut_max2_r', 'glut_max3_r'],
    'glut_med1+2+3': ['glut_med1_r', 'glut_med2_r', 'glut_med3_r'],
    'glut_min1+2+3': ['glut_min1_r', 'glut_min2_r', 'glut_min3_r'],
    'ALL gluteals (9)': ['glut_max%d_r' % i for i in (1, 2, 3)]
    + ['glut_med%d_r' % i for i in (1, 2, 3)]
    + ['glut_min%d_r' % i for i in (1, 2, 3)],
    'glut_max1+med1+min1': ['glut_max1_r', 'glut_med1_r', 'glut_min1_r'],
    'add_mag1+2+3': ['add_mag1_r', 'add_mag2_r', 'add_mag3_r'],
    'vasti (3)': ['vas_med_r', 'vas_int_r', 'vas_lat_r'],
    'iliacus+psoas': ['iliacus_r', 'psoas_r'],
    'hamstrings (4)': ['bifemlh_r', 'semiten_r', 'semimem_r', 'grac_r'],
}
print('--- candidate group sums ---')
for label, members in groups.items():
    total = sum(mif[m] for m in members if m in mif)
    print(f'{label:22s} {total:8.1f}')
