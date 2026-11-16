"""Bounded read-only mechanical-source audit. No mj_forward/mj_step or simulation."""
from pathlib import Path
import ast,csv,hashlib,json,xml.etree.ElementTree as ET
import numpy as np
import mujoco

HERE=Path(__file__).resolve().parent
ROOT=next(p for p in HERE.parents if p.name=='Bipedal_Robot')
SP=ROOT/'Code/MuJoCo_SNS/spinal'
OS=ROOT/'Solid_Models/OpenSim/Gait2392_Robotbody'
MODEL=OS/'mjc/gait2392_simbody/gait2392_simbody_cvt3.xml'
sources=[SP/'_muscle_force_compare.py',SP/'muscle_force_compare.csv',SP/'runner.py',SP/'bsolve_ik.py',SP/'bsolve_out.npz',OS/'gait2392_simbody.osim',OS/'subject01_simbody.osim',OS/'ResultsBSolve/zz_bsolve_StaticOptimization_force.sto',OS/'ResultsBSolve/zz_bsolve_StaticOptimization_activation.sto',MODEL]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
rows=list(csv.reader((SP/'muscle_force_compare.csv').open()))
stats=[dict(muscle=r[0],R2=float(r[1]),corr2=float(r[2]),lag=int(r[3]),mean_ratio=float(r[4]),reference_peak_N=float(r[5]),mj_peak_N=float(r[6]),R2_act1=float(r[7]))for r in rows[1:]]
passing=[r['muscle']for r in stats if .87<=r['R2']<=.98 and .9<=r['mean_ratio']<=1.1]
z=np.load(SP/'bsolve_out.npz',allow_pickle=False);t=z['t']
def sto(p):
    lines=p.read_text(errors='replace').splitlines();i=next(i for i,s in enumerate(lines)if s.strip().lower()=='endheader')
    names=lines[i+1].split();data=np.array([[float(x)for x in s.split()]for s in lines[i+2:]if s.strip()])
    return names,data
names,force=sto(OS/'ResultsBSolve/zz_bsolve_StaticOptimization_force.sto')
for r in stats:
    current_peak=float(np.interp(t,force[:,0],force[:,names.index(r['muscle'])]).max())
    r['current_reference_interpolated_peak_N']=current_peak
    r['saved_vs_current_reference_peak_abs_difference_N']=abs(current_peak-r['reference_peak_N'])

# Execute only the current pure XML transformer isolated from the runner AST.
# No imports of the runner, no writes to its model directory, no model dynamics.
tree=ast.parse((SP/'runner.py').read_text())
function=next(n for n in tree.body if isinstance(n,ast.FunctionDef)and n.name=='patch_xml')
module=ast.Module(body=[function],type_ignores=[])
namespace={'np':np,'mujoco':mujoco,'PRUNE_MUSCLES':{'quad_fem_r','quad_fem_l','gem_r','gem_l','peri_r','peri_l'}}
exec(compile(module,str(SP/'runner.py'),'exec'),namespace)
raw_text=MODEL.read_text();raw=ET.fromstring(raw_text)
model=mujoco.MjModel.from_xml_path(str(MODEL));data=mujoco.MjData(model)
patched_text=namespace['patch_xml'](raw_text,model,data)
patched_air_text=namespace['patch_xml'](raw_text,model,data,no_ground=True)
patched=ET.fromstring(patched_text);air=ET.fromstring(patched_air_text)
raw_joints={x.attrib['name']:x for x in raw.find('worldbody').iter('joint')if 'name'in x.attrib}
pat_joints={x.attrib['name']:x for x in patched.find('worldbody').iter('joint')if 'name'in x.attrib}
raw_act={x.attrib['name']:x for x in raw.findall('./actuator/general')}
pat_act={x.attrib['name']:x for x in patched.findall('./actuator/general')}
osim=ET.parse(OS/'gait2392_simbody.osim')
osim_fmax={x.attrib['name']:float(x.findtext('max_isometric_force'))for x in osim.iter()if x.tag in ['Thelen2003Muscle','Millard2012EquilibriumMuscle']}
repairs={'hip_axes':{},'force_scale_entries':{},'disabled_actuators':{},'driver_limits':{},'patella_routes':{},'pathpoint_followers':[]}
for n in ['hip_flexion_r','hip_flexion_l','hip_adduction_r','hip_adduction_l']:
    repairs['hip_axes'][n]={'OLD':raw_joints[n].attrib['axis'],'NEW':pat_joints[n].attrib['axis']}
for n in pat_act:
    base=n.rsplit('_',1)[0]
    if base in ['ercspn','intobl','extobl','ext_hal']:
        repairs['force_scale_entries'][n]={'OLD_active':float(raw_act[n].attrib['gainprm'].split()[2]),'NEW_active':float(pat_act[n].attrib['gainprm'].split()[2]),'passive_unchanged':float(raw_act[n].attrib['biasprm'].split()[2])==float(pat_act[n].attrib['biasprm'].split()[2]),'source_OpenSim_Fmax_N':osim_fmax[n]}
    if n in namespace['PRUNE_MUSCLES']:
        repairs['disabled_actuators'][n]={'active':float(pat_act[n].attrib['gainprm'].split()[2]),'passive':float(pat_act[n].attrib['biasprm'].split()[2])}
for side in ['r','l']:
    for base in ['knee_angle','hip_flexion','ankle_angle','subtalar_angle','mtp_angle','hip_adduction','hip_rotation']:
        n=base+'_'+side;repairs['driver_limits'][n]={'OLD':raw_joints[n].attrib['limited'],'NEW':pat_joints[n].attrib['limited'],'range_unchanged':raw_joints[n].attrib.get('range')==pat_joints[n].attrib.get('range')}
    name='rect_fem_'+side+'_tendon'
    repairs['patella_routes'][name]={'OLD':[x.attrib['site']for x in raw.find('./tendon/spatial[@name="'+name+'"]')],'NEW':[x.attrib['site']for x in patched.find('./tendon/spatial[@name="'+name+'"]')]}
for eq in raw.findall('./equality/joint'):
    n=eq.attrib['joint1']
    repairs['pathpoint_followers'].append({'joint':n,'OLD_range':raw_joints[n].attrib.get('range'),'NEW_limited':pat_joints[n].attrib.get('limited'),'NEW_armature':pat_joints[n].attrib.get('armature')})
repairs['contact_pairs']={'OLD':len(raw.findall('./contact/pair')),'NEW_ground':len(patched.findall('./contact/pair')),'NEW_air':len(air.findall('./contact/pair'))}
repairs['pathpoint_collision_flags']={'raw_virtual_flags':sum(x.attrib.get('contype')=='2'and x.attrib.get('conaffinity')=='2'for x in raw.findall('.//geom')),'patched_virtual_flags_still2':sum(x.attrib.get('contype')=='2'and x.attrib.get('conaffinity')=='2'for x in patched.findall('.//geom'))}
out={'scope':'CSV summary/primary script/OpenSim reference/pose-clock audit, plus current static compilation and pure XML-transform inspection; no model forward calculation or integration.',
    'csv_header_columns':len(rows[0]),'csv_data_columns':sorted(set(len(r)for r in rows[1:])),
    'csv_header':rows[0],'actual_row_order':['muscle','R2_matched','R2shape_matched','lag_frames','mj_over_os_mean','OpenSim_peak_N','MuJoCo_matched_peak_N','R2_act1'],
    'force_stats':stats,'passing_eight':passing,'frame_count':len(t),'first_last_time_s':[float(t[0]),float(t[-1])],'median_frame_dt_s':float(np.median(np.diff(t))),
    'saved_npz_arrays':{k:list(z[k].shape)for k in z.files},'current_unpatched_compiled_actuator_count':int(model.nu),
    'model_repairs':repairs,'source_sha256':{str(p):sha(p)for p in sources if p.is_file()}}
(HERE/'HISTORICAL_MECHANICAL_RESULTS.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k]for k in ['csv_header_columns','csv_data_columns','passing_eight','frame_count','first_last_time_s','median_frame_dt_s','current_unpatched_compiled_actuator_count']},indent=2))
print('reference_peak_max_abs_difference_N',max(x['saved_vs_current_reference_peak_abs_difference_N']for x in stats))
print('repair_counts',len(repairs['force_scale_entries']),len(repairs['disabled_actuators']),len(repairs['driver_limits']),len(repairs['pathpoint_followers']),repairs['contact_pairs'])
