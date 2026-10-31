from pathlib import Path
import json,hashlib
from lxml import etree as E

base=Path(__file__).resolve().parents[2]
for folder,name,title in [('W2L_VirtualWalker_20261031','W2L_VirtualWalker','Original controller'),('W2L_VirtualWalker_Newer_20261031','W2L_VirtualWalker_Newer','II/contact/bilateral-RG controller')]:
 w=base/folder;b=json.loads((w/'BUILD.json').read_text());v=json.loads((w/'VERIFICATION.json').read_text());test=json.loads((w/'TEST_RESULTS.json').read_text())
 audit={}
 for source in b['source']:
  src=Path(source);out=w/(name+('.aproj' if src.suffix=='.aproj' else '_Standalone.asim'));rs=E.parse(str(src));ro=E.parse(str(out))
  muscles=lambda r:[E.tostring(x) for x in r.iter('RigidBody') if x.findtext('Type','').startswith('LinearHill')]
  drawings=lambda r:[x.text for x in r.iter('DiagramXml')]
  for text in drawings(ro):
   if text and text.strip():E.fromstring(text.encode())
  audit[str(src)]=dict(neural=v['preservation'][str(src)],muscle_and_stretch_receptor_payloads_unchanged=muscles(rs)==muscles(ro),neural_diagram_text_unchanged=drawings(rs)==drawings(ro),neural_diagrams_parsed=len(drawings(ro)),output_sha256=hashlib.sha256(out.read_bytes()).hexdigest())
  assert audit[str(src)]['muscle_and_stretch_receptor_payloads_unchanged'] and audit[str(src)]['neural_diagram_text_unchanged']
 (w/'PRESERVATION.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
 poses=b.get('poses_deg',{})
 gap_text='; '.join(f"{n}: {(y-b['floor_top_y'])*1000:.3f} mm" for n,y in b['foot_sole_y'].items())
 text=f'''W2L passive hip-wing recording copy: {title}

Open {name}.aproj in AnimatLab. This is a separate copy; the earlier recording projects and source models are unchanged. Play runs 10.1 seconds at normal speed (MatchPhysicsStep). Reset between recordings. The corresponding {name}.asim and {name}_Standalone.asim are supplied.

The inherited corrected pose is {json.dumps(poses)}. Tibial/femoral face alignment, hinge-reference changes, and preserved anatomical joint limits are documented by the geometry agent in the source geometry folder. Hip and knee positions are asymmetric as requested. Ankle poses remain inside the source model's original ranges rather than forcing an incompatible level-foot posture.

The Root pelvis is NOT frozen and has no added joint. Two cyan boxes are rigid compound children of Root, extending laterally along z. Each is 120 mm anterior-posterior by 25 mm high by 300 mm lateral, centered at root coordinates (0,-40,+/-250) mm, with declared mass 0.001 g each. Their coordinates and GUI matrices are synchronized.

Two static, invisible, collision-enabled rails run along world x. Each rail is 4 m long, 40 mm thick, and 180 mm wide, centered at z=+/-320 mm. Their upper surface is y={b['rail_top_y']:.9f} m. The existing moving floor is preserved and level, with top y={b['floor_top_y']:.6f} m. AnimatLab axes here are x anterior/travel, y superior, and z lateral.

Initial physical box sole clearances: {gap_text}. The right leg is intentionally in swing. Wing undersides start {1000*b['wing_to_rail_gap']:.3f} mm above the rails, 12 mm more than the lowest stance sole clearance. The lower foot contact-sensor surface can be slightly below the physical foot box. Native first foot-sensor contact times (seconds, null means no contact during this test): {json.dumps(v['first_foot_sensor_contacts_s'])}. Wing geometry reaches the rails at approximately {v['first_approximate_wing_contact_s']} seconds, reconstructed from body positions and a 0.5 mm threshold, not a measured contact force.

Only the new rail/wing material has sliding friction 0.05 and zero restitution. This is passive physical support, not a hidden balance constraint. Pitch, roll, yaw, and translation remain mechanically possible. Preliminary friction trials in the original-copy build folder used the earlier hinge-reference setup and are not comparisons of the final model.

Native verification: exit {test['exit_code']}; 5001 data samples over 0-10 seconds, run completed at 10.101 seconds. Root Y stayed {v['root_y_range_m'][0]:.6f}-{v['root_y_range_m'][1]:.6f} m and ended at {v['end_root_y_m']:.6f} m. Net root X displacement was {v['root_x_displacement_m']:+.6f} m. Motor-neuron and joint activity were recorded. These are supported controller experiments; they do not establish unsupported walking, balanced gait, or right-leg weight bearing.

See TEST_RESULTS.json, VERIFICATION.json, BUILD.json, and PRESERVATION.json for signals, clearances, source hashes, neural/muscle preservation, and drawing-XML checks. Raw signals are in verification.npz and Virtual walker verification.txt. The original-copy folder includes verification_comparison.png comparing Root height, hip/knee joint changes, and hip motor-neuron voltages in both controllers. New chart angle values are changes relative to the initialized pose, not absolute anatomical angles.

GUI-open verification was not performed, to avoid disrupting the user's open AnimatLab application. Project XML, embedded drawing XML, linked resources, and native standalone runtime were checked. This pass does not create a screen recording.

Rebuild from the repository root with D:\\Anaconda\\envs\\myo\\python.exe Neuromechanical_Models/W2L_VirtualWalker_20261031/tools/build_virtual_walker.py (add --newer for the newer controller). Tests and plots are generated by run_test.py and plot_and_verify.py. Generated projects, traces, and graphics are ignored by Git.
'''
 (w/'README.txt').write_text(text,encoding='utf-8')
 print(folder,'preservation checked')
