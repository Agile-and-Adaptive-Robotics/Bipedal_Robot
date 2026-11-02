import path from 'node:path';
export async function apply(h){
 const {d,rec,begin,newSlide,im,tx,line,addAfter,changed,a,GREEN,INK}=h;
 const f=name=>a(path.join('mechanics',name));
 begin(4,'Actuator force and moment arm determine joint torque',
  'Sutton and colleagues model a leg joint as a pendulum with muscle actuation, elastic forces, and damping. This is the paper’s generic leg-joint schematic, not a detailed human knee. The highlighted force and moment arm explain the conversion from actuator force to torque. This defense next asks whether compliant BPAs can provide enough force and contraction to reproduce biological joint movement. Source: Sutton, Szczecinski, Quinn, and Chiel (2023), PNAS Nexus 2, pgad298, Figure 1B. https://doi.org/10.1093/pnasnexus/pgad298',
  'Sutton et al. (2023), PNAS Nexus 2:pgad298, Fig. 1B');
 im(4,f('Sutton2023_fig1B_swing.png'),55,111,490,353);
 tx(4,'Joint torque',590,160,320,42,28,true,GREEN);
 tx(4,'M = F r⊥',580,224,340,74,44,true,INK);
 tx(4,'Force × moment arm',566,326,370,45,28,false,INK);
 for(const n of [8,31]){
  begin(n,'The robot’s available torque covers only part of the human range',
   'Historical human and robot biceps-femoris-long-head torque calculations are reproduced from the original geometry and force scripts and original input MAT files. Each retains its own original joint-angle domain. Both panels use identical torque limits and color limits; differences are not produced by separate scaling. These are the earlier baseline model calculations, not the later corrected redesign or a new experimental result. The coverage gap motivates improved BPA force prediction and joint design. See mechanics source manifest for generating scripts, input hashes, and shared scales.',
   'Bolen (2019), original BFLH calculations · shared torque and color scales');
  tx(n,'Human model',110,110,300,33,26,true,GREEN,'center');tx(n,'Robot model',560,110,300,33,26,true,GREEN,'center');
  im(n,f('historical_bflh_common_scales.png'),23,147,914,326);
 }
 begin(29,'OpenSim muscle paths determine the human knee-torque target',
  'The actual Gait2392 knee and biceps-femoris-short-head path are shown at a specified knee coordinate. Force acts along the muscle path. OpenSim computes the moment arm for that generalized knee coordinate; multiplying by muscle force gives the muscle contribution to knee torque. The knee includes coupled translations, so a fixed hinge-center approximation must not be substituted silently for its generalized moment arm. The overlay is generated from the actual model geometry and state, with source and numeric checks recorded in the mechanics manifest.',
  'Gait2392 OpenSim model · biceps femoris short head · Dissertation Ch. 4');
 im(29,f('opensim_knee_BFSH_overlay_landscape.png'),45,105,870,367);
 begin(33,'Improved BPA force prediction does not resolve the torque mismatch',
  'This single example is the 48 cm file-label pinned-flexor test, with measured unloaded rest length 48.5 cm. The plotted torque and the next slide’s measured BPA length come from the same saved test. The rigid geometry is kept, while the BPA force characterization is improved; the residual discrepancy motivates loaded-length corrections. Measured blue points are drawn above hybrid predictions. The saved source identifies the separate load-cell and fish-scale measurements and their pressure readings. No correction factors were refitted.',
  'Pinned-flexor test, l₀ = 485 mm · saved 48 cm test record');
 im(33,f('pinned48_preoptimization_torque.png'),30,107,900,366);
 begin(34,'The measured BPA is shorter than the rigid-geometry prediction',
  'Same pinned-flexor test as the preceding torque plot. Blue points are the measured inflated BPA lengths in the source InflatedLength array, with the same measured knee angles. The dashed curve is the geometric prediction from the source MuscleLength array. All lengths are converted from meters to millimeters. The mismatch motivates a signed constant offset, effective compliance, and bending corrections; the following slides explain these terms before their identification. This is not a comparison across different actuators or tests.',
  'Same pinned-flexor test · measured and predicted active BPA length');
 im(34,f('pinned48_measured_vs_predicted_length.png'),30,107,900,366);
 // Shorten only the editable illustrative Δl overlay; source geometry is untouched.
 const arc=rec.get(36).elements.filter(e=>e.kind==='shape'&&e.position?.left>840&&e.position?.top>330&&e.position?.top<400);
 if(arc.length)d.delete(arc.map(e=>e.id));
 const cx=818.75,cy=330.32,r=63,start=.12,end=.12+.2*(1.05-.12);let prev=null;
 for(let k=0;k<=10;k++){const p=start+(end-start)*k/10,point=[cx+r*Math.cos(p),cy+r*Math.sin(p)];if(prev)line(36,...prev,...point,false,'#CD34B5',4);prev=point;}
 const dl=rec.get(36).elements.find(e=>e.kind==='textbox'&&e.text==='Δl');if(dl)d.resolve(dl.id).position={left:884,top:333,width:60,height:37};
 d.slides.getItem(35).speakerNotes.text+=' The magenta illustrative length-loss segment is shortened to one-fifth of the previous annotation; it is not a new measured bending loss.';changed.push(36);
 begin(45,'The single-BPA flexor falls short of the human torque target',
  'The original dissertation/torque-paper biomimetic-flexor panel compares a 20 mm BPA with the selected maximum human BFSH torque. The original rigid-body prediction and corrected prediction are shown along with measured torque; reduced-pressure points are separately marked. The corrected model still does not make the earlier single-actuator route meet the target. This is the reason to change the route optimization. Use the original panel and legend; do not describe the fitted model as a physical optimization.',
  'Dissertation Ch. 4 · single-BPA biomimetic flexor validation');
 im(45,f('biomimetic_flexor_publication.png'),90,105,780,369,{left:.5,top:0,right:0,bottom:0});
 begin(48,'The redesigned flexor meets the human target at the tested angles',
  'The torque-paper publication figure compares the human BFSH reference, redesigned prediction at 620 kPa, measured torque, and prediction at the measured pressure. The measured blue points are opaque and drawn above the yellow hybrid predictions. The physical validation contains 15 measured points from -119 to +5 degrees; measured margins over the selected human magnitude are +0.42 N m at -105 degrees to +17.0 N m at -14 degrees. These discrete isometric tests do not establish dynamic gait or complete-leg behavior. Pressure differences explain why the measured points are not expected to lie on the 620 kPa curve.',
  'Torque paper, redesigned 20 mm flexor · measured pressure used for hybrid prediction');
 im(48,f('redesigned_flexor_publication.png'),75,107,810,366);
 const cam=newSlide('cam_force_pressure','The cam changes actuator loading and the pressure needed for torque',
  'Both curves use the finalized saved direct-route and cam designs with the same adopted design factors. The source verifies that the saved cam baseline equals the saved direct-route design. The left panel compares total BPA force at a fixed 620 kPa supply; this is not the force required to produce the same torque. The right panel compares pressures back-calculated to match the same human torque target, from the saved October 28 pressure table. Cam pressure is lower through much of the middle range but not at every angle; do not claim a uniform reduction. This slide is a model comparison, followed by physical cam-extensor validation. No route optimization was rerun.',
  'Saved direct and cam designs · force at 620 kPa; pressure for the same human target');
 im(cam,f('cam_force_and_pressure_comparison.png'),25,111,910,356);addAfter(50,cam);
 begin(51,'The cam-assisted extensor meets the selected human torque target',
  'Use the finalized dissertation/torque-paper physical validation, preserving the source model curves and measured data. Blue measurements are drawn after yellow predictions at measured pressure, so the measured result remains visible. The human reference is the selected vastus medialis muscle, not the complete human extensor group. The model capacity and measured-pressure predictions are distinct. These results are isometric measurements at the tested angles, not a walking demonstration.',
  'Dissertation Ch. 4 · cam-assisted extensor physical validation');
 im(51,f('cam_extensor_physical_publication.png'),45,108,870,366);
 const reverse=newSlide('reverse_pulley','Design factors support the design of a more biomimetic knee',
  'Original dissertation reverse-pulley knee schematic. Two parallel BPAs drive a moving sheave. In the ideal two-part cable arrangement, tendon travel is twice BPA shortening and tendon tension is half the combined BPA force. The corrected loaded BPA force and length are inputs to this design calculation. Friction, packaging, clearance, and load-sharing still require physical checks. This is a mechanism-design example and transition to the whole-leg coverage comparison, not a claim that this particular knee has been dynamically validated.',
  'Dissertation reverse-pulley knee mechanism · ideal cable relation');
 im(reverse,f('reverse_pulley_biomimetic.png'),25,107,500,366);
 tx(reverse,'Two parallel BPAs',565,150,345,42,29,true,GREEN);
 tx(reverse,'2× tendon travel',565,232,345,45,30,true,INK);
 tx(reverse,'½ combined BPA force',565,307,365,70,29,false,INK);
 // Place the mechanism immediately before the former whole-leg comparison.
 addAfter(51,reverse);
 const remove=rec.get(52).elements.filter(e=>e.position?.top>=100&&e.position?.top<478&&e.kind!=='chart').map(e=>e.id);if(remove.length)d.delete(remove);
 const chart=rec.get(52).elements.find(e=>e.kind==='chart');if(!chart)throw Error('Native whole-leg chart missing');d.resolve(chart.id).position={left:90,top:132,width:780,height:296};
 tx(52,'Angular-domain coverage, averaged across 27 sites',95,438,770,35,25,false,INK,'center');changed.push(52);
}
