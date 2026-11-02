import fs from 'node:fs/promises';
import fss from 'node:fs';
import path from 'node:path';
import {FileBlob,PresentationFile} from 'file:///C:/Users/Ben%20Bolen/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/@oai/artifact-tool/dist/artifact_tool.mjs';
const W=path.resolve('.'), prior=path.resolve('../C3NS_Revision_20261102_v19');
const d=await PresentationFile.importPptx(await FileBlob.load(path.join(prior,'output/Bolen_PhD_Defense_RoseGreen_20261102_v19.pptx')));
const snap=await d.inspect({select:'$.slides[*]',include:['id','index','elements.id','elements.kind','elements.position','elements.text']});
const rec=new Map(snap.records.map(x=>[x.index+1,x]));
const GREEN='#174D24',INK='#202A25';
const changed=[];
const a=n=>path.join(W,'assets',n), old=n=>path.join(prior,'assets',n);
const fig=n=>path.resolve('../C3NS_Revision_20261031_v18/submitted_figs',n);
function tx(n,text,x,y,w,h,size=25,bold=false,color=INK,alignment='left'){
 const s=d.slides.getItem(n-1).shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
 s.text=text;s.text.insets=0;s.text.style={fontFamily:'Arial',fontSize:size,bold,color,alignment};return s;
}
function im(n,file,x,y,w,h,crop){
 const blob=fss.readFileSync(file),ext=path.extname(file).toLowerCase();
 if(crop&&(ext==='.png'||ext==='.gif')){
  const iw=(ext==='.png'?blob.readUInt32BE(16):blob.readUInt16LE(6))*(1-crop.left-crop.right),ih=(ext==='.png'?blob.readUInt32BE(20):blob.readUInt16LE(8))*(1-crop.top-crop.bottom),sc=Math.min(w/iw,h/ih);
  const nw=iw*sc,nh=ih*sc;x+=(w-nw)/2;y+=(h-nh)/2;w=nw;h=nh;
 }
 return d.slides.getItem(n-1).images.add({blob,contentType:(ext==='.jpg'||ext==='.jpeg')?'image/jpeg':ext==='.gif'?'image/gif':'image/png',fit:'contain',alt:path.basename(file),position:{left:x,top:y,width:w,height:h},...(crop?{crop}:{})});
}
function line(n,x1,y1,x2,y2,arrow=true){
 return d.slides.getItem(n-1).shapes.add({geometry:'connector',kind:'straight',position:{left:Math.min(x1,x2),top:Math.min(y1,y2),width:Math.abs(x2-x1),height:Math.abs(y2-y1),horizontalFlip:x2<x1,verticalFlip:y2<y1},line:{fill:GREEN,width:2.5},...(arrow?{tail:{type:'triangle',width:'med',length:'med'}}:{})});
}
function begin(n,title,notes,cite){
 const ids=rec.get(n).elements.filter(e=>e.position&&e.position.top<478).map(e=>e.id);if(ids.length)d.delete(ids);
 tx(n,title,18,19,850,82,title.length>85?30:32,true,GREEN);
 d.slides.getItem(n-1).speakerNotes.text=notes;
 const source=rec.get(n).elements.find(e=>e.position?.top>=478&&e.position?.top<510&&e.position?.width>650&&e.kind==='textbox');
 if(cite&&source)d.resolve(source.id).text=cite;
 changed.push(n);
}
function block(n,text,x,y,w,h){d.slides.getItem(n-1).shapes.add({geometry:'rect',position:{left:x,top:y,width:w,height:h},fill:'#F3F8F3',line:{fill:GREEN,width:2}});tx(n,text,x+12,y+16,w-24,h-20,27,true,GREEN,'center');}
// Preserve the original GIF objects, timings, and media relationships on slide 2.
const title2=rec.get(2).elements.find(e=>e.kind==='textbox'&&e.position?.top<80&&e.position?.width>650);
d.resolve(title2.id).text='Search-and-rescue robots must move reliably over uneven terrain';
d.slides.getItem(1).speakerNotes.text='Wheels are an effective choice for warehouse robots operating on smooth, controlled floors. Search-and-rescue robots face a different problem: rubble, obstacles, uncertain contact, and balance. Biology offers useful mechanical and neural principles for that application. These clips motivate the research; they do not establish an energy-efficiency advantage for this robot. Sources: Amazon Robotics; DARPA Robotics Challenge, 2015.';
changed.push(2);
begin(5,'Neural control depends on feedback from the moving body','This is a functional control diagram. Descending drive enters the neural controller, whose motor output activates the actuators. Body motion and contact change muscle length, muscle force, and foot-contact signals; sensory pathways return those signals to the controller. No single negative-feedback summing junction is implied: the neural pathways include excitatory and inhibitory connections. Source: submitted dissertation, motivation and neuromechanical models.','Dissertation Chapters 1 and 6');
block(5,'Neural\ncontroller',110,205,240,100);block(5,'Actuators\nand body',610,205,240,100);block(5,'Sensory\npathways',385,372,230,92);
tx(5,'Descending drive',100,113,270,36,25,false,INK);line(5,230,155,230,205);
line(5,350,255,610,255);tx(5,'Motor output',375,211,210,36,25,false,INK,'center');
line(5,730,305,730,418,false);line(5,730,418,615,418);
tx(5,'Muscle length, force,\nand foot contact',745,326,200,78,22,false,INK);
line(5,385,418,230,418,false);line(5,230,418,230,305);
tx(5,'Sensory\nfeedback',35,337,175,62,23,false,INK,'center');
begin(6,'Biomimetic robots let us test hypotheses about human movement','A second research goal is to use biomimetic robots to study biomechanics and neuroscience. We can change actuator mechanics or sensory feedback and observe the resulting movement under controlled conditions. The human model supplies anatomical reference geometry; the robot provides an experimental system.','Dissertation Chapters 1 and 2');
im(6,fig('Gait2392_human_reference_oblique.png'),80,115,350,310);im(6,fig('aarl_bpa_biped.png'),560,115,300,310);
tx(6,'Human musculoskeletal model',35,437,440,36,24,true,GREEN,'center');tx(6,'Biomimetic robot',520,437,380,36,24,true,GREEN,'center');
begin(7,'Biological mechanics and control also inform prosthetic-limb design','A third, longer-term goal is to improve prosthetic and orthotic devices. This photograph shows Hugh Herr wearing two BiOM T2 powered ankle prostheses, an established example of biologically inspired mechanics and control. It is external motivating work, not a device developed or clinically tested in this dissertation. Source: MIT News, Rob Matheson, April 17, 2014, https://news.mit.edu/2014/hugh-herr-bionic-ankle-emulates-nature-0417 . Photograph: Bryce Vickmark. MIT News makes this photograph available under CC BY-NC-ND; reproduced without alteration.','BiOM T2 ankle prostheses · Photo: Bryce Vickmark / MIT News (2014)');
im(7,a('biom_prostheses.jpg'),115,112,730,350);
begin(11,'Three gaps guide the design and testing of a biomimetic biped','The main talk addresses three questions: predicting BPA force across initial lengths; predicting and validating joint torque while accounting for compliance and bending; and testing neural controllers with mechanical models and sensory feedback.','Main-talk roadmap');
im(11,a('bpas.png'),35,160,265,175);im(11,a('knee_cutout_approved.png'),345,160,265,175);im(11,fig('W2LNEW_TREADMILL_STILL.png'),650,140,270,215);
['Predict BPA force','Predict and validate\njoint torque','Test neuromechanical\nmodels'].forEach((t,i)=>tx(11,t,30+310*i,383,285,78,26,true,GREEN,'center'));
function mechanicsMap(n,title,question){
 begin(n,title,question+' This section concerns the actuator and joint mechanics.','Dissertation Chapters 3 and 4');
 im(n,a('bpas.png'),35,151,330,200);im(n,a('knee_cutout_approved.png'),455,169,445,191);
 tx(n,'BPA force',35,115,330,40,28,n===12,GREEN,'center');tx(n,'Joint torque',490,115,355,40,28,n===28,GREEN,'center');
 tx(n,question,75,400,810,68,30,true,GREEN,'center');
}
mechanicsMap(12,'Braided Pneumatic Actuator Modeling','How much force can a selected BPA produce?');
mechanicsMap(28,'Isometric Knee Torque','Does the selected BPA produce the required joint torque?');
begin(14,'Joint torque depends on force and its perpendicular moment arm','The dissertation background gives this fixed-axis geometric construction. F is the force, r_k is its perpendicular moment arm, and M_z is the moment about the out-of-page joint axis. The photograph supplies physical context. The biomimetic four-bar knee has a moving instantaneous center, so its visible gold bolt is not labeled as a fixed joint axis.','Dissertation Ch. 2, force moment and moment-arm projection');
im(14,a('moment_arm_methods.png'),25,109,535,350,{left:0,top:0,right:.5,bottom:0});
const knee=im(14,a('knee_cutout_approved.png'),555,220,350,166);knee.rotation=270;
begin(17,'Existing models mispredict force for a 450 mm BPA','This is the 20 mm diameter, 450 mm initial-length comparison from the published force paper. Only the existing Sarosi and M-B predictions are shown here, grouped as Others, to establish the modeling gap. The gray reference is equality between measured and predicted force. The proposed model is introduced after the measurements and appears in later validation plots. Underlying points are retained from the published comparison.','Bolen et al. (2026), Actuators 15(5):230 · 20 mm BPA, initial length 450 mm');
im(17,a('slide17_comparison.png'),35,105,890,365);
begin(19,'We first measure maximum force at each initial BPA length','Before modeling the full pressure-contraction relation, we examine maximum isometric force at 620 kPa as a function of uninflated resting length. The test frame fixes the actuator length while the load cell measures force. The next slide presents the 10 mm actuators; the 20 mm results follow separately. Source: submitted Chapter 3, Maximum Force at 620 kPa.','Bolen et al. (2026), Actuators 15(5):230');
im(19,old('s018_img0.png'),35,110,390,360);
tx(19,'Hold pressure at 620 kPa',475,170,440,65,29,true,GREEN);
tx(19,'Vary initial BPA length, l₀',475,272,440,65,29,true,GREEN);
tx(19,'Measure maximum force, F₆₂₀',475,374,450,65,29,true,GREEN);
for(const [n,title,file,eq,note] of [
 [20,'Shorter 10 mm BPAs produce less maximum force','s020_img0.png','max_force_b.png','The general arctangent form is F620(l0)=b1 atan[b2(l0-b3)]. The coefficients determine the force scale, length scale, and length offset. This slide shows the 10 mm results; the next slide shows the 20 mm results.'],
 [22,'Pressure and initial length jointly determine maximum BPA force','s022_img0.png','max_force_a.png','The pressure-dependent model uses a1 and a2 with P in both the prefactor and the arctangent argument: Fmax(l0,P)=a1 P atan[a2 P(l0-0.0075 m)]. The 7.5 mm offset is fixed in the published equation, not introduced here as a newly fitted coefficient.']]){
 begin(n,title,note+' Source: Bolen et al. (2026), submitted Chapter 3.','Bolen et al. (2026), Actuators 15(5):230');
 im(n,old(file),150,104,660,320);im(n,a(eq),n===20?200:135,435,n===20?560:690,34);
}
begin(30,'Muscle geometry determines the joint’s moment arm','Original dissertation background figure, unchanged. Panel A shows the force-moment construction; panel B shows the signed moment-arm projection used in the torque calculation. The joint axis points out of the page. The force is in the illustrated plane; out-of-plane forces require the knee-axis projection stated in Chapter 4.','Dissertation Ch. 2, equivalent geometric constructions of joint torque');
im(30,a('moment_arm_methods.png'),30,106,900,365);
begin(35,'A signed offset corrects the predicted BPA length','Chi0 is a signed constant correction to predicted current BPA length, not to measured resting length l0 and not an elastic deflection. lm=lLMT-2lfitting-chi0. lLMT is the routed origin-to-insertion length and lfitting is an end-fitting length. The adopted offsets are +3.94 mm for the flexor and -6.38 mm for the extensor. No geometric length discrepancy is dimensioned on the stock photograph. Source: Chapter 4, equation Lm and adopted parameter table.','Dissertation Ch. 4, length offset and system stiffness');
im(35,a('bpas.png'),130,115,700,130,{left:0,top:.5,right:0,bottom:0});
im(35,a('length_offset.png'),115,290,730,65);
tx(35,'χ₀ accounts for the unloaded length discrepancy',85,399,790,60,28,true,GREEN,'center');
// Retain the existing bending illustration and physical attachment callout, but remove
// every unsupported photo dimension from the chi0 column.
const leftIds=rec.get(37).elements.filter(e=>e.position&&e.position.top>=110&&e.position.top<478&&e.position.left<325).map(e=>e.id);d.delete(leftIds);
tx(37,'χ₀: length offset',30,120,275,40,25,true,GREEN);im(37,a('bpas.png'),30,185,275,64,{left:0,top:.5,right:0,bottom:0});
im(37,a('length_offset.png'),25,304,288,45);tx(37,'Constant offset',35,400,270,40,24,true,GREEN,'center');
d.slides.getItem(36).speakerNotes.text='Chi0 corrects the unloaded length discrepancy and is signed. The source equation defines the correction; the stock BPA photograph is not a dimensional drawing. Chi1 and chi2 describe effective system stiffness, including brackets, fixtures, artificial tendons, and loading hardware. The middle callout identifies the attachment system, not a coordinate axis. Chi3 scales length loss due to bending along the wrapped R alpha arc. Source: dissertation Chapter 4 and Appendix C.';
changed.push(37);
begin(39,'The design factors improve all three measures of pinned-flexor fit','Native table reproduces the dissertation pinned-flexor goodness-of-fit values before and after correction, with all three reported measures. Fitting tests: 48.5 cm and 41.0 cm. Validation tests: 45.7 cm and 40.6 cm with a 3.5 cm artificial tendon. Means are the source table means, not newly computed pooled residuals. RMSE and maximum absolute residual are in N m; FVU is dimensionless and lower is better. Source: Chapter 4, Table FlxPinGof.','Dissertation Ch. 4, pinned-knee flexor goodness of fit');
tx(39,'Before correction',229,118,348,35,26,true,GREEN,'center');tx(39,'With design factors',577,118,348,35,26,true,GREEN,'center');
const gof=[
 ['l₀; lₜₑₙ (cm)','RMSE','FVU','Max.','RMSE','FVU','Max.'],
 ['Fit: 48.5; —','7.617','2.1278','25.531','2.165','0.1718','5.467'],
 ['Fit: 41.0; —','9.694','2.1149','21.076','1.084','0.0277','1.766'],
 ['Fit mean','8.656','2.1213','23.303','1.624','0.0998','3.617'],
 ['Val: 45.7; —','8.803','1.9380','26.264','1.368','0.0468','3.614'],
 ['Val: 40.6; 3.5','6.373','0.5935','16.859','1.611','0.0427','5.891'],
 ['Validation mean','7.588','1.2657','21.562','1.490','0.0447','4.752']
];
const table=d.slides.getItem(38).tables.add({rows:gof.length,columns:7,left:35,top:164,width:890,height:270,columnWidths:[194,116,116,116,116,116,116],values:gof});
table.cells.block({row:0,column:0,rowCount:7,columnCount:7}).assign({textStyle:{fontFamily:'Arial',fontSize:23,color:INK,alignment:'center'},margins:{left:3,right:3,top:5,bottom:5}});
table.borders.assign({fill:'#D4DFD6',width:.7});
for(let r=0;r<7;r++)for(let c=0;c<7;c++){
 table.getCell(r,c).fill=r===0?GREEN:(r===3||r===6)?'#E7F1E8':'#FFFFFF';
 table.getCell(r,c).text.style={fontFamily:'Arial',fontSize:23,color:r===0?'#FFFFFF':INK,alignment:'center',bold:r===0||r===3||r===6};
}
tx(39,'RMSE and Max.: N·m     FVU: dimensionless',130,447,700,30,21,false,INK,'center');
begin(36,'Bracket stiffness is expressed in the local attachment frame','The dissertation CAD schematics show the local bracket-frame axes for the pinned flexor and the biomimetic flexor with a 20 mm BPA. The fitted stiffness values represent the combined attachment system. The next method slide explains the coordinate rotations before force is projected into these directions. Source: Chapter 4, Figure bktFrame, panels A and D.','Dissertation Ch. 4, bracket-frame orientations');
im(36,a('bktFrame.png'),40,120,410,320,{left:0,top:0,right:.5,bottom:.5});
im(36,a('bktFrame.png'),515,120,405,320,{left:.5,top:.5,right:0,bottom:0});
tx(36,'Pinned-knee flexor',35,442,420,35,24,true,GREEN,'center');tx(36,'Biomimetic-knee flexor',500,442,430,35,24,true,GREEN,'center');
begin(43,'Route optimization must satisfy the torque target with loaded geometry','The revised search varies origin and insertion positions, BPA resting length, and artificial-tendon length. Every candidate uses the identified compliance and bending terms before predicting torque across knee angle. The worst-angle shortfall is weighted most strongly so a good mean cannot compensate for failure at one angle. The required envelope is the selected human muscle magnitude plus 5 percent, in the required torque direction. Pressure, contraction, bone clearance, wrapping, and tendon feasibility are checked. Surrogate optimization supplies candidates and pattern search refines them. Source: Chapter 4, Constrained Knee-Route Redesign.','Dissertation Ch. 4, constrained knee-route redesign');
tx(43,'Choose attachment points and muscle–tendon lengths',65,119,830,46,28,true,GREEN,'center');
im(43,a('route_variables.png'),155,176,650,57);
tx(43,'Evaluate corrected torque throughout the range of motion',60,268,840,50,27,true,GREEN,'center');
im(43,a('route_target.png'),145,333,670,58);
tx(43,'Pressure · contraction · clearance · wrapping',80,427,800,38,26,false,INK,'center');
function newSlide(n,title,notes,cite){
 d.slides.add({width:960,height:540});rec.set(n,{elements:[]});begin(n,title,notes,cite);
 tx(n,cite,34,481,850,28,10,false,'#52615A');tx(n,String(n),918,518,30,20,13,false,'#FFFFFF','center');
}
newSlide(124,'Two rotations align the bracket frame with the attachment point','Original dissertation Figure xiFrameTransform, split into its two panels without changing the geometry. First rotate about the fixed space-frame z axis toward the sagittal projection. Then pitch about the current bracket y axis to point the bracket x axis toward the attachment. The adopted transform is Rsb=Rz(theta_z)Ry(theta_y)^T. Frame orientation and the configuration-specific stiffness-axis assignments must be kept distinct.','Dissertation Ch. 4, bracket-frame transformation');
im(124,a('xiFrameTransform.png'),25,152,440,290,{left:0,top:.15,right:.08,bottom:.49});
im(124,a('xiFrameTransform.png'),490,152,440,290,{left:0,top:.63,right:.08,bottom:.025});
tx(124,'1. Rotate toward the sagittal projection',35,112,420,66,24,true,GREEN,'center');
tx(124,'2. Tilt toward the attachment point',500,112,420,66,24,true,GREEN,'center');
newSlide(125,'Force and attachment displacement must be solved together','The source figure distinguishes three-dimensional bracket displacement from the scalar displacement along the BPA path. The two projected bracket deflections and the tendon extension act in series. BPA force depends on contraction, which itself changes with delta. The evaluator uses fzero to solve Fmax FBPA(epsilon*(delta),P*) - keq delta = 0 at each knee angle, then updates the attachment positions, routed length, force direction, moment arm, and torque. This inner force-balance solve is distinct from parameter identification and route optimization.','Dissertation Ch. 4, series compliance and force-balance equation');
im(125,a('xiProjection.png'),40,112,330,330);
tx(125,'Force changes attachment displacement',420,150,495,75,29,true,GREEN);
tx(125,'Displacement changes BPA length and force',420,270,495,85,29,true,GREEN);
im(125,a('force_balance.png'),350,394,570,67);
newSlide(126,'Rigid-body route predictions overestimate the physical knee’s torque','Enlarged right-hand panel from the original slide 42: 20 mm biomimetic-knee flexor. The route looked adequate under rigid-body assumptions, but the physical measurements did not meet the selected human torque target over the full range. The measured shortfall motivated modifying the route optimization to include the identified design factors and feasibility requirements. The source legend uses Unoptimized/Optimized for the torque-model corrections, not for the later redesigned route. Reduced-pressure measurements are identified separately.','Dissertation Ch. 4, biomimetic-knee flexor validation');
im(126,old('s042_img3.png'),105,111,750,360,{left:.49,top:0,right:0,bottom:0});
newSlide(127,'The redesigned flexor route can be checked across knee motion','Original repository animation of the two-BPA 20 mm knee-flexor route. This is a model animation, followed by the physical validation results on the next slide. The animation makes attachment locations, wrapping, and the predicted torque and pressure variation visible. Source: Code/Matlab/Mesh_Optimization/Results/Knee_Flexor_20mm.gif.','Repository model animation · redesigned 20 mm knee flexor');
im(127,'D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Results/Knee_Flexor_20mm.gif',45,110,870,360);
newSlide(128,'The direct extensor route met the model target but overloaded hardware','Original direct-route extensor animation, not the later cam route. The generator loads Vas_Pam_20mm_Result.mat and uses the direct-route distal-ring geometry. The direct route met the predicted torque envelope at actuator forces that broke printed and machined parts. The following slide explains why the cam was added to increase moment arm and reduce the required BPA force; the physical cam results follow.','Repository model animation · direct 20 mm knee extensor');
im(128,'D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Results/Knee_Extensor_20mm.gif',45,110,870,360);
await (await import('./pass2_content.mjs')).apply({d,rec,begin,newSlide,im,tx,line,block,a,old,fig,changed,GREEN,INK,W});
await (await import('./pass3_content.mjs')).apply({d,rec,begin,im,tx,a,old,changed,GREEN,INK});
const order=[...Array.from({length:37},(_,i)=>i+1),124,125,129,38,39,40,41,42,126,43,127,45,128,46,47,44,48,...Array.from({length:75},(_,i)=>i+49)].filter(n=>![14,50].includes(n));
await fs.writeFile('build/slide_order.json',JSON.stringify(order));
await fs.writeFile('build/changed_slides.json',JSON.stringify([...new Set(changed)].sort((a,b)=>a-b)));
await (await PresentationFile.exportPptx(d)).save(path.join(W,'build/artifact_draft.pptx'));
console.log('Authored slides '+changed.join(', '));
