import fs from 'node:fs/promises';
import path from 'node:path';
export async function apply(h){
 const {d,rec,begin,newSlide,im,tx,line,block,a,old,fig,changed,GREEN,INK,W}=h;
 const snap=await d.inspect({select:'$.slides[*]',include:['id','index','elements.id','elements.kind','elements.position','elements.text']});
 rec.clear();for(const r of snap.records)rec.set(r.index+1,r);
 await fs.writeFile('build/pass2_elements.json',JSON.stringify(snap.records));
 const textOf=e=>typeof e.text==='string'?e.text:JSON.stringify(e.text??'');
 function replace(n,rx,t,pos){const e=rec.get(n).elements.find(e=>e.kind==='textbox'&&rx.test(textOf(e)));if(!e)throw Error('Missing text '+n+' '+rx);const obj=d.resolve(e.id);obj.text=t;if(pos)obj.position=pos;changed.push(n);return obj;}
 function note(n,t){d.slides.getItem(n-1).speakerNotes.text += '\n'+t;changed.push(n);}
 function clearWhere(n,p){const ids=rec.get(n).elements.filter(p).map(e=>e.id);if(ids.length)d.delete(ids);changed.push(n);}
 function white(n,x,y,w,h){d.slides.getItem(n-1).shapes.add({geometry:'rect',position:{left:x,top:y,width:w,height:h},fill:'#FFFFFF',line:{fill:'none',width:0}});}
 // Background follows the causal problem before the goals and roadmap.
 begin(3,'BPAs provide lightweight, compliant tensile actuation','A braided pneumatic actuator expands radially and shortens when pressurized. Its low actuator mass, large tensile force, and compliance make it useful for muscle-like attachment arrangements. The Puppy photograph is an example of pneumatic actuation in a biomimetic robot, not an energy-efficiency experiment. Air supply and control hardware must be included in any system-level efficiency claim. Source: submitted dissertation background and Bolen, Movement and Control of Biomimetic Humanoid Robots, slide 9. Transition: the force is useful only if limb geometry produces the required joint torque.','Dissertation Ch. 2 · AARL Puppy robot from Bolen’s earlier talk');
 im(3,a('bpas.png'),35,150,390,220);im(3,path.join(W,'build/supervisor_sources/old_s09_p0.jpg'),465,130,455,315);
 tx(3,'Pressure → shortening and tension',35,398,390,64,26,true,GREEN,'center');
 begin(4,'Muscle force and its moment arm determine joint torque','The force vector F and perpendicular moment arm rk determine the torque Mz about the joint axis. This is the actual source schematic from the dissertation background. The physical photo is context; the moving four-bar knee center is not identified with a visible bolt. Transition: feedback determines when and how strongly the actuators pull.','Dissertation Ch. 2, force-moment construction');
 im(4,a('moment_arm_methods.png'),40,115,520,350,{left:0,top:0,right:.5,bottom:0});
 const knee=im(4,a('knee_cutout_approved.png'),562,230,350,166);knee.rotation=270;
 begin(6,'The robot reproduces selected human anatomy with pneumatic actuation','The human and robot models share the intended musculoskeletal organization, while the robot uses selected BPA paths and mechanical joints. This is an anatomical mapping, not a claim of identical force or movement. Transition: matching anatomy does not guarantee matching contraction range or joint torque.','Dissertation Ch. 2, human muscle paths and robot BPA paths');
 im(6,a('human_gait2392_paths.png'),65,118,290,310);im(6,a('robot_pam_routes.png'),410,118,290,310);im(6,fig('aarl_bpa_biped.png'),745,210,170,215);
 tx(6,'Human muscle paths',35,438,340,35,25,true,GREEN,'center');tx(6,'Robot BPA paths',380,438,340,35,25,true,GREEN,'center');tx(6,'Physical robot',738,438,190,35,22,true,GREEN,'center');
 begin(7,'Limited BPA contraction constrains the robot’s joint range','Historical model comparison from the older Movement and Control talk. The blue region is the OpenSim reference and the red region is the robot model. BPA contraction limits and mechanical interference constrain achievable hip angles; the plot is not a new physical range measurement and contraction is not the only limiting factor. Transition: joint torque can also fall short within the attainable range.','Bolen, Movement and Control of Biomimetic Humanoid Robots · historical model comparison');
 im(7,path.join(W,'build/supervisor_sources/old_s14_p0.png'),240,106,480,365);
 begin(8,'The earlier robot model did not reproduce human knee torque','Historical biceps-femoris-long-head comparison in the submitted background. These are model predictions over hip and knee angles, not current physical measurements. The source uses different vertical and color ranges, which are retained and identified; compare the numerical torque values, not apparent surface height alone. Transition: these mismatches motivate both better robot design and controlled tests of biological principles.','Dissertation Ch. 2 · Bolen (2019), biceps femoris long head');
 tx(8,'Human model',50,108,405,35,26,true,GREEN,'center');tx(8,'Earlier robot model',505,108,405,35,26,true,GREEN,'center');
 im(8,fig('masters_torque_comparison.png'),25,151,440,270,{left:.045,top:.332,right:.50,bottom:.334});
 im(8,fig('masters_torque_comparison.png'),505,151,430,270,{left:.535,top:.332,right:.015,bottom:.334});
 tx(8,'Original torque and color scales differ',120,440,720,35,23,false,INK,'center');
 begin(9,'Biology guides robot design; robot experiments test movement hypotheses','The first goal is to use biological principles to improve robots for unstructured environments. The reciprocal goal is to use increasingly biomimetic robots for controlled interventions that would be difficult or impossible in humans. For example, alter one actuator property or sensory pathway, measure the resulting movement, and generate a biological prediction to test. The return arrow represents a research process; this dissertation does not claim a new human experiment. Transition: the same understanding can inform devices that assist people.','Dissertation Ch. 1, research motivation and objectives');
 im(9,a('human_gait2392_paths.png'),70,171,220,240);im(9,fig('aarl_bpa_biped.png'),670,171,220,240);
 tx(9,'Biological evidence',35,113,305,39,27,true,GREEN,'center');tx(9,'Biomimetic robot',615,113,310,39,27,true,GREEN,'center');
 line(9,300,212,650,212);tx(9,'Mechanics and neural principles',315,155,325,55,24,false,INK,'center');
 line(9,650,365,300,365);tx(9,'Controlled interventions →\ntestable biological predictions',306,290,343,65,24,false,INK,'center');
 tx(9,'Improve autonomous movement in unstructured environments',100,442,760,35,24,true,GREEN,'center');
 begin(10,'Biological mechanics and control also inform prosthetic-limb design','A third, longer-term goal is to improve prosthetic and orthotic devices. Hugh Herr is shown wearing BiOM T2 powered ankle prostheses, an external example of biologically inspired mechanics and control. This dissertation does not claim a new clinical device. MIT News, April 17, 2014: https://news.mit.edu/2014/hugh-herr-bionic-ankle-emulates-nature-0417 . Photo: Bryce Vickmark / MIT News, CC BY-NC-ND, reproduced without alteration. Transition: this dissertation addresses the underlying actuator, joint, and modeling-tool gaps.','BiOM T2 ankle prostheses · Photo: Bryce Vickmark / MIT News (2014)');
 im(10,a('biom_prostheses.jpg'),115,112,730,350);
 begin(11,'Three gaps limit prediction and testing of a biomimetic biped','The first gap is force prediction across initial BPA lengths, pressure, and contraction. The second is predicting joint torque when the hardware bends and attachment systems deform. The third is testing neural controllers with complex mechanical models in inspectable, reusable environments. These define the main talk.','Main-talk research questions');
 im(11,a('bpas.png'),35,157,265,175);im(11,a('knee_cutout_approved.png'),345,157,265,175);im(11,fig('W2LNEW_TREADMILL_STILL.png'),655,145,265,190,{left:0,top:.12,right:0,bottom:0});
 ['Force across lengths,\npressures, and contraction','Joint torque with\ncompliance and bending','Tests of coupled\nneural and body models'].forEach((t,i)=>tx(11,t,30+310*i,373,285,85,25,true,GREEN,'center'));
 // Source14 is omitted from the order: its joint schematic now appears at4.
 const lc=replace(13,/Contracted\s+at 620 kPa/,'Free contraction\nat 620 kPa');lc.text.style={fontFamily:'Arial',fontSize:23,color:GREEN,bold:true};
 // Closing the same 450 mm example shown before the proposed model.
 begin(27,'The new model improves force predictions for the 20 mm BPAs','The lower two panels of the published comparison show the 300 mm and 450 mm initial lengths. The 450 mm case closes the discrepancy introduced earlier. Proposed-model predictions are purple; existing model series retain the original plot colors here. These are comparisons against the source measurements, not a new test.','Bolen et al. (2026), Actuators 15(5):230');
 im(27,old('s027_img0.png'),30,115,900,350,{left:0,top:.50,right:0,bottom:.02});
 note(32,'This historical biarticular example is biceps femoris long head (bifemlh); it is distinct from the short-head redesign target introduced earlier.');
 tx(18,'Current BPA\nlength, lₘ',670,290,245,80,27,true,GREEN);line(18,655,329,585,329);changed.push(18);
 replace(19,/Hold pressure at/,'Hold the BPA at lₘ = l₀');replace(19,/Vary initial BPA/,'Pressurize to 620 kPa');replace(19,/Measure maximum/,'Measure F₆₂₀ across initial lengths');
 note(19,'Maximum force is measured with the actuator held at its unloaded rest length, before allowing shortening. l_m=l_0 for this reference measurement.');
 for(const n of [35,37])note(n,'The displayed length equation is the no-artificial-tendon case; subtract artificial-tendon length when one is present.');
 // Place the requested photo with the actual pinned test, not a different apparatus.
 clearWhere(40,e=>e.kind==='image'&&e.position?.left>790&&e.position?.top>90);
 const p=im(40,'D:/GitHub/Bipedal_Robot/Pictures/PictureAnalysis/20220307_222908.jpg',826,170,105,165);p.rotation=180;
 tx(40,'χ₁, χ₂ fixed; χ₀, χ₃ refitted',210,455,570,28,22,false,INK,'center');
 note(42,'The plotted biomimetic-knee tests use the original dissertation apparatus inset. The earlier route optimization appeared adequate under rigid-body assumptions, but physical tests failed to meet the human target, so the route optimization was modified.');changed.push(42);
 // Model display is large enough to see motion from the auditorium.
 begin(127,'The redesigned flexor route can be checked across knee motion','Original complete-cycle repository GIF from Knee_Flexor_20mm.gif, followed by physical validation. The native PowerPoint crop enlarges the existing animation without changing any frame. It illustrates the selected route and wrapping; quantitative validation is the following measured torque plot. The cached animation does not have a generation manifest establishing byte-identical measured-test geometry.','Repository model animation · redesigned 20 mm knee flexor');
 im(127,'D:/GitHub/Bipedal_Robot/Knee_Flexor_20mm.gif',110,112,430,355,{left:465/1074,top:125/646,right:(1074-685)/1074,bottom:(646-530)/646});
 tx(127,'BPA 1',600,185,270,42,28,true,'#CD34B5');tx(127,'BPA 2',600,255,270,42,28,true,'#9D02D7');tx(127,'Human BFSH',600,325,285,42,28,true,'#FA8775');
 begin(128,'The direct extensor route met the model target but overloaded hardware','This is the original direct-route extensor animation, not the later cam. Its generator loads Vas_Pam_20mm_Result.mat. Physical testing of the direct route failed printed and machined hardware. The cam then increased moment arm and reduced the required actuator force. This native crop preserves every GIF frame.','Repository model animation · direct 20 mm knee extensor');
 im(128,'D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Results/Knee_Extensor_20mm.gif',130,112,430,350,{left:520/1074,top:190/646,right:(1074-675)/1074,bottom:(646-445)/646});
 tx(128,'Direct BPA route',605,210,270,66,28,true,GREEN);tx(128,'Next design change:\nincrease the moment arm',585,330,320,90,27,false,INK);
 replace(126,/Rigid-body route predictions/,'The single-BPA flexor falls short of the human torque target');
 white(126,397,281,200,44);tx(126,'Rigid-body model',399,282,196,22,20,true,INK);tx(126,'Corrected model',399,304,196,22,20,true,INK);
 note(126,'The single20mm flexor remains below the selected human target at high flexion, including the reported−103° to−78° region. Legend wording distinguishes correction-factor fitting from the next route-redesign step.');
 begin(43,'Route optimization includes compliance and bending','The search varies attachment positions, BPA resting length, and artificial-tendon length. At every joint angle, fixed identified design factors and the inner displacement solve determine corrected force and torque. The target is the selected human torque magnitude plus5percent in the required direction. The worst-angle shortfall is penalized; pressure, contraction, bone clearance, wrapping, and tendon feasibility constrain the design. surrogateopt finds candidates and patternsearch refines them. The animation illustrates the geometric variables; measured validation follows.','Dissertation Ch. 4, constrained knee-route redesign');
 tx(43,'Choose attachment positions\nand actuator–tendon lengths',35,130,515,68,27,true,GREEN,'center');im(43,a('route_variables.png'),40,220,510,44);
 line(43,290,279,290,317);tx(43,'Predict torque at every joint angle',35,331,515,65,26,true,GREEN,'center');im(43,a('route_target.png'),35,399,520,42);
 im(43,'D:/GitHub/Bipedal_Robot/Knee_Flexor_20mm.gif',595,113,305,320,{left:465/1074,top:125/646,right:(1074-685)/1074,bottom:(646-530)/646});
 tx(43,'Check pressure, contraction, clearance, and wrapping',65,451,830,27,23,false,INK,'center');
 begin(125,'Force and attachment displacement must be solved together','Original source detail of one bracket. The3-D bracket displacement is projected along the BPA force. At both ends these projections, plus tendon extension, contribute to series displacement delta. At each pose and fixed chi, fzero solves the scalar force balance. The F620(l0) multiplier is the maximum force reference at620kPa; normalized FBPA contains pressure dependence. This avoids double pressure scaling. The result updates attachment positions, routed length, force direction, moment arm, and torque. Parameter identification is a separate outer problem on the next slide.','Dissertation Ch. 4, displacement projection and force balance');
 im(125,a('xiProjection.png'),35,112,405,280,{left:.03,top:.52,right:.34,bottom:.02});
 tx(125,'Force displaces the attachments',475,144,440,70,28,true,GREEN);tx(125,'Displacement changes\nBPA length and force',475,265,440,80,28,true,GREEN);
 tx(125,'One bracket detail',75,386,335,30,22,false,INK,'center');im(125,a('force_balance.png'),65,428,830,44);
 // Distinguish the inner force balance from identification and route redesign.
 newSlide(129,'The design factors are fitted before testing new configurations','At each candidate parameter set chi, the inner scalar root solves displacement separately at every pose. gamultiobj identifies chi from torque residuals. Flexor: fit 48.5 and 41.0 cm; validation 45.7 and 40.6 cm with 3.5 cm tendon. Extensor: retain the same chi1/chi2 pair, fit chi0/chi3 on 43,43+tendon,46 cm; validate on40,40+tendon,47 cm. Validation results also informed Pareto acceptance, so these are not blind external tests. The adopted parameter set is then fixed for route redesign. Source: Chapter4, AppendixC, evaluator fold records.','Dissertation Ch. 4 and Appendix C · parameter identification');
 tx(129,'Fit the parameters',55,125,390,40,28,true,GREEN,'center');tx(129,'Evaluate other tests',515,125,390,40,28,true,GREEN,'center');
 block(129,'Flexor: χ₀, χ₁, χ₂\n2 fitting tests',50,194,360,100);block(129,'2 flexor\nvalidation tests',550,194,360,100);line(129,410,244,550,244);
 block(129,'Extensor: χ₀, χ₃\n3 fitting tests',50,328,360,100);block(129,'3 extensor\nvalidation tests',550,328,360,100);line(129,410,378,550,378);
 tx(129,'Shared χ₁, χ₂ retained for the extensor',140,447,680,28,24,true,GREEN,'center');
 // Scientific condition labels and cross-model provenance.
 begin(53,'The revised walker separates left and right rhythms and adds feedback','The newer bilateral walker has separate left/right rhythm generators, commissural connections, type-II length feedback, and heel/toe contact feedback. The figures are original source crops: bilateral rhythm generators and commissural/contact pathways. They are partial views of a model, not a complete biological circuit. Length feedback reports limb configuration; contact feedback reports interaction with the surface.','Dissertation Ch. 6 · newer bilateral W2L controller');
 im(53,old('s053_img0.png'),35,130,270,285,{left:.065,top:.025,right:.80,bottom:.55});
 im(53,old('s053_img0.png'),350,190,570,175,{left:.24,top:.47,right:.23,bottom:.30});
 tx(53,'Left / right rhythm',35,423,275,35,24,true,GREEN,'center');tx(53,'Contact and commissural pathways',335,389,590,35,24,true,GREEN,'center');
 tx(53,'Type-II length feedback is also retained',345,439,565,30,22,false,INK,'center');
 note(56,'Separate recovered-W2L experiment, not the preceding W2L108 demonstration. Controller, drive, and support are fixed across feedback conditions. Contact-phase numbers are medians over 4.5–30 s; the displayed traces are a12–15 s excerpt. A phase of0.5 denotes alternation.');
 tx(56,'Recovered W2L · same controller, drive, and support',100,104,760,28,21,false,INK,'center');
 replace(56,/^Feedback intact$/,'Feedback intact',{left:145,top:138,width:330,height:35});replace(56,/^Feedback disabled$/,'Feedback disabled',{left:540,top:138,width:330,height:35});
 replace(56,/Blue solid:/,'Blue solid: left     Orange dashed: right',{left:240,top:175,width:600,height:30});
 replace(56,/Contact phase: 0.498/,'Median phase: 0.498 cycle');replace(56,/Contact phase: 0.223/,'Median phase: 0.223 cycle');
 tx(56,'Median over 4.5–30 s · 0.5 = alternation',250,469,660,22,17,false,INK,'center');
 replace(77,/0\.71.*left/,'Duty: 0.71 right · 0.65 left · reference 0.61');
 note(77,'Duty-side labels corrected against overlay_metrics.json: right0.710142857; left0.648571429. Underlying figure/data unchanged.');
 note(78,'Each human and BPA force curve is normalized to its own peak; this comparison concerns curve shape, not absolute force agreement.');
 for(const rx of [/Vastus lateralis:/,/Semimembranosus:/]){const e=rec.get(78).elements.find(e=>rx.test(textOf(e)));d.resolve(e.id).text.style={fontFamily:'Arial',fontSize:24,color:GREEN,bold:true,alignment:'center'};}
 tx(78,'Each curve is divided by its own peak',150,115,700,32,23,false,INK,'center');
 replace(100,/Table 6\.1.*contribution 7/,'Deng et al. (2019) · Dissertation Table 6.1, literature-controller tools');
 tx(101,'Blue: stimulus · Dashed: baseline',210,445,660,30,22,false,INK,'center');changed.push(101);
 replace(103,/Source circuit/,'Imported rat-model organization');replace(103,/Specified stimulus test/,'Separate reduced-circuit test');
 for(const e of rec.get(103).elements.filter(e=>e.kind==='textbox'&&e.position?.top===438))d.resolve(e.id).text.style={fontFamily:'Arial',fontSize:23,color:GREEN,bold:true,alignment:'center'};
 note(103,'The left diagram organizes the imported public biarticular rat model. The right result is a separate reduced14-neuron Deng-derived deletion experiment. This is not a matched replication of the complete imported model.');
 note(44,'Coverage values are equal-site means of trapezoid-weighted angular-domain coverage, not unweighted fractions of sampled points.');
 replace(44,/Coverage of sampled/,'Angular-domain coverage\n(mean across 27 sites)');
 replace(44,/Dissertation Ch\. 4, Isometric/,'27-site model · 620 kPa · actuator sizing and pulley ratios · inherited routes');
 d.slides.getItem(22).speakerNotes.text='The free-contraction measurements show substantial variance; the source does not establish a conclusive predictive relationship with initial length. The maximum contraction must be measured for each actuator. Source: published Chapter3, maximum contraction and discussion.';changed.push(23);
 note(24,'Published model development used537 measured points, with80percent for fitting and20percent for testing. Counts are kept in notes to respect the requested visual focus. The selected-length comparisons shown later are not claimed to be a separate held-out-actuator study.');
 note(63,'Reduced14-neuron Deng-derived neural circuit. The displayed RG current pulse is+50nA for3.876s, two baseline periods. This is a neural-response experiment.');
 tx(63,'Reduced Deng circuit · RG pulse: +50 nA for 3.876 s',75,111,825,30,22,false,INK,'center');changed.push(63);
 await (await import('./sadb_questions_fix.mjs')).apply(h);
}
