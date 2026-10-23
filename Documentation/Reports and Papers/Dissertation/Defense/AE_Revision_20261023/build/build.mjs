import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath,pathToFileURL} from 'node:url';
import {Presentation,PresentationFile} from '@oai/artifact-tool';
const dir=path.dirname(fileURLToPath(import.meta.url)), root=path.dirname(dir);
const manifest=JSON.parse(await fs.readFile(path.join(dir,'asset_manifest.json'),'utf8'));
const p=Presentation.create({slideSize:{width:1280,height:720}});
p.theme.defaultFont='Calibri';
const specs=[];
const C={title:'#AF7B51',text:'#222222',muted:'#595959',teal:'#00796B'};
function text(s,t,x,y,w,h,size=26,color=C.text,bold=false){
 const a=s.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
 a.text=t;a.text.style={typeface:'Calibri',fontSize:size,bold,color};a.text.insets=0;return a;
}
async function image(s,key,x,y,w,h){
 const source=manifest[key]?.asset??path.join(root,'assets',key+'.png');
 s.images.add({blob:new Uint8Array(await fs.readFile(source)),contentType:'image/png',alt:key.replaceAll('_',' '),fit:'contain',position:{left:x,top:y,width:w,height:h}});
}
function table(s,values,y=190,h=330){
 const t=s.tables.add({rows:values.length,columns:values[0].length,left:50,top:y,width:1180,height:h,values});
 for(let r=0;r<values.length;r++)for(let c=0;c<values[0].length;c++){
  const cell=t.getCell(r,c);cell.text.style={typeface:'Calibri',fontSize:25,color:r===0?C.teal:C.text,bold:r===0};cell.fill='#FFFFFF';
 }
 return t;
}
function base(title,section,cite,notes,backup=false){
 const s=p.slides.add();s.background.fill='#FFFFFF';
 const n=p.slides.items.length;
 text(s,section.toUpperCase(),30,15,1200,24,16,C.teal,true);
 text(s,title,30,49,1220,110,43,C.title,true);
 text(s,cite,35,672,1160,35,18,C.muted);
 text(s,String(n).padStart(2,'0'),1210,675,40,25,16,C.muted);
 s.speakerNotes.text=notes+'\n\nSources: '+cite;
 if(backup)s.hidden=true;
 specs.push({slide:n,title,section,cite,notes,backup,assets:[]});
 return s;
}
async function fig(title,section,cite,notes,key,callout='',backup=false){
 const s=base(title,section,cite,notes,backup);
 await image(s,key,45,175,1190,callout?420:480);
 if(callout)text(s,callout,65,614,1150,43,26,C.teal,true);
 specs.at(-1).assets=[key];return s;
}
async function pair(title,section,cite,notes,a,b,labels=[]){
 const s=base(title,section,cite,notes);
 await image(s,a,40,175,585,440);await image(s,b,655,175,585,440);
 labels.forEach((v,i)=>text(s,v,50+i*615,621,565,40,25,C.teal,true));
 specs.at(-1).assets=[a,b];return s;
}
let s=base('Movement and Control of\nBiomimetic Humanoid Robots','Dissertation defense','Portland State University · Agile & Adaptive Robotics Laboratory',
 'Opening, 30 seconds. This dissertation develops and tests components needed to study a biomimetic robot under neural control. It establishes actuator and joint results and selected simulation capabilities. It does not demonstrate harness-free walking. Advisor: M. Anthony Hunt.');
text(s,'Ben P. Bolen',65,260,650,60,38,C.text,true);text(s,'Mechanical Engineering PhD defense\nAdvisor: M. Anthony Hunt, Ph.D.\nNovember 6, 2026',65,340,700,160,29);await image(s,'robot',850,170,330,465);specs.at(-1).assets=['robot'];

s=base('Four gaps connect actuator mechanics to neural-control experiments','The argument','Revision C · Ch. 1, research gaps and objectives',
 'Allow 1 minute. Read the argument in chapter order: predict force across actuator lengths; predict the torque delivered through compliant routes; recover source-grounded sensory evidence; and connect neural circuits to mechanical models with representative verification. These are project-specific gaps, not claims that no related tools exist elsewhere. A research robot is the long-term goal. Each contribution needs its own evidence.');
table(s,[['Knowledge gap','Contribution tested here'],['Force across actuator lengths','Resting-length-aware force model'],['Torque through compliant routes','Identified corrections and redesigned-route tests'],['Scattered sensory evidence','Linked, searchable SADb resource'],['Untested model connections','Interfaces with matched circuit checks']],175,420);

await pair('This dissertation adds models and tests to an existing robot program','Scope and attribution','Revision C · Ch. 1; Ch. 6, Model Conversion and Actuator Coupling',
 '1 minute. Credit the AARL body plan, Morrow’s earlier robot and the Li walker mechanics, Steele’s biomimetic knee, and Deng’s published controller. OpenSim supplies human anatomy and targets; MyoConverter supplies initial conversion; MuJoCo supplies physics; SNS-Toolbox supplies neuron and synapse formulations. My contributions include actuator characterization with collaborators, correction identification, route-design and testing work with collaborators, SADb curation and application, model repairs and audits, and graphical circuit/interface implementations. The 92-actuator model is a modeling scope, not a physical robot completed in this dissertation. Ask me to explain the specific change and evidence rather than implying I invented these packages.','robot','hardware',['Existing anatomy and mechanisms','New route modeling and tests']);

await fig('Isometric tests measured force across length, pressure, and contraction','Actuator characterization','Bolen, Elzein, Pang, and Hunt (2026), Actuators · Revision C, Ch. 3',
 '1 minute. Introduce a braided pneumatic actuator (BPA). Characterization used 10 mm and 20 mm Festo actuators, up to 620 kPa, with 537 force-pressure-contraction observations: 321 for 10 mm and 216 for 20 mm. Resting lengths: 112–518 mm for 10 mm, 300–509 mm for 20 mm. These are isometric measurements. Credit collaborators as coauthors. Valve dynamics and hysteresis are outside this static model.','jig','Braided pneumatic actuators (BPAs): 537 observations');
await pair('Maximum isometric force decreases as resting length becomes shorter','Actuator result','Bolen et al. (2026), Actuators · Revision C, Fig. Fmax_fit',
 '1 minute. These existing MATLAB plots retain the experimental data and model curves. Explain the x-axis as resting length and the y-axis as force at 620 kPa. Filled points are experimental data, the solid line is the new fit, and the dashed manufacturer curve differs most for short actuators. Avoid a universal novelty claim beyond the published scope.','max10','max20');

s=base('At 112 mm, the new model reduces force error from 42.1% to 4.5%','Actuator result','Revision C · Ch. 3, Effects of Resting Length; 10 mm BPA at 620 kPa',
 '1 minute. This is one specific measured case, not an average error or an independent holdout claim. Measured force 350.9 N, new model 335.3 N, manufacturer tool 498.6 N. Errors quoted verbatim from Chapter 3. The chart is a native editable chart reconstructed from the reported numeric values, not digitized from an image.');
s.charts.add('bar',{position:{left:120,top:185,width:1030,height:420},categories:['Measured','New model','Manufacturer tool'],series:[{name:'Force (N)',values:[350.9,335.3,498.6],fill:'#2E5A88'}],hasLegend:false,dataLabels:{showValue:true,position:'outEnd',textStyle:{fontSize:25}},yAxis:{title:'Force (N)',textStyle:{fontSize:23}},xAxis:{textStyle:{fontSize:23},title:'10 mm diameter, 112 mm resting length, 620 kPa'}});
text(s,'350.9 N measured',65,620,1150,40,27,C.teal,true);

await pair('Normalized force curves describe both tested actuator diameters','Actuator model','Bolen et al. (2026), Actuators · Revision C, Fig. GenForceFit',
 '1 minute. Define dimensionless force, pressure, and contraction verbally. Each diameter has its own fitted coefficients. The complete 537-point set is not shown in these representative isoclines. The source used an 80/20 training/validation split. A random point split does not establish generalization to every new actuator or dynamic condition. The force model also requires each actuator’s measured maximum contraction.','surface10','surface20');
await fig('The new model has the lowest error in the reported model comparisons','Actuator evidence','Bolen et al. (2026), Actuators · Revision C, Fig. Comparison and model-comparison table',
 '1 minute. Explain the measured-versus-predicted axes and diagonal. Four matched diameter/length cases compare this model with Sárosi and Martens–Boblan predictions. Scope the assertion to these reported comparisons. Do not imply exhaustive superiority over every BPA model. Complete statistics and coefficients remain in the dissertation. This ends the actuator section: a force prediction still needs a transmission model before it predicts joint torque.','comparison');

await fig('Pinned and biomimetic knees separate identification from transfer tests','Joint torque: experimental design','Revision C · Ch. 4, test stands; four-bar knee credited to Steele et al.',
 '1 minute. The pinned knee supports identifying transmission effects. The biomimetic four-bar knee has a moving instantaneous center, so load-cell forces must be transformed to that center. Separate force characterization, correction identification, route optimization, and physical tests throughout this section. These apparatus views are existing figure components.','knee_jig');
s=await fig('Measured actuator length exposes errors in the rigid-route prediction','Joint torque: diagnosis','Revision C · Ch. 4, Fig. FlxPin_preopt; hybrid-method description',
 '1 minute. The hybrid method substitutes measured actuator length into the characterized force calculation while retaining the joint mechanics. Its closer agreement with measured torque implicates transmission/length assumptions. This is diagnostic evidence, not proof that the force model has zero error. Show where the yellow hybrid points approach the blue measured points while the rigid prediction diverges.','hybrid');
text(s,'Torque (N·m)',45,170,220,35,24,C.text,true);
await fig('The corrected model accounts for length offset, compliance, and wrapping','Joint torque: model contribution','Revision C · Ch. 4, Improved Torque Model; Appendix C',
 '1 minute. Chi0 is a constant length offset. Chi1 and chi2 encode effective system stiffness, lumping brackets, fixtures, and the cable winch; they are not independently measured bracket material properties. Chi3 represents lost active length during wrapping. One stiffness pair serves the adopted flexor and extensor configurations. The flexor identification has a flat stiffness valley, so the fitted pair is not unique. The externally imposed redesign overshoot criterion helped choose the adopted pair.','bracket','Effective system stiffness includes the fixtures and transmission');

s=base('Identification reserves tests before assessing corrected predictions','Joint torque: identification','Revision C · Ch. 4, Tables FlxPinGof and ExtPinGof; Appendix C',
 '1 minute. gamultiobj provides the Pareto fronts for identification. Report adopted-member fold assignments: flexor training 48.5 and 41.0 cm, validation 45.7 and 40.6 cm with tendon. The unused 47.9 cm test is a separate post hoc evaluation with encoder correction. Extensor uses three training and three validation tests; three early tests with bracket shifts are excluded. Row 77 of the flexor front and row 32 of the pair-locked extensor front are the picks of record. Route redesign later uses surrogateopt plus patternsearch, a different optimization task.');
table(s,[['Adopted fold','Training','Held out','Not used'],['Pinned flexor','2 tests','2 tests','1 test'],['Pinned extensor','3 tests','3 tests','3 tests']],205,270);
text(s,'Identify corrections: gamultiobj',65,525,1120,40,29,C.teal,true);
text(s,'Redesign routes: surrogateopt, then patternsearch',65,580,1120,45,29,C.teal,true);

s=await fig('The corrected flexor model reproduces the held-out torque trend','Joint torque: representative evidence','Revision C · Ch. 4, Fig. FlxPin_optimized; Table FlxPinGof',
 '1 minute. Both panels are held-out validation tests of the current adopted fold. The adopted full two-bracket model reduces mean RMSE from 8.656 to 1.624 N·m for training and 7.588 to 1.490 N·m for validation. RMSE means root-mean-square error. Figure traces illustrate the pattern; they do not replace per-test statistics. The nominal panel labels 46 cm and 40 cm tendon correspond to actual lengths 45.7 cm and 40.6 cm, as explained by the current caption.','corrected','Held-out root-mean-square error: 7.588 to 1.490 N·m');
text(s,'Torque (N·m)',45,170,220,35,24,C.text,true);

s=base('Corrections reduce held-out error in both pinned-knee configurations','Joint torque: validation','Revision C · Ch. 4, Tables FlxPinGof and ExtPinGof',
 '1 minute. Native table uses the exact reported validation means, N·m. Baseline means the characterized BPA model without transmission corrections. This table does not compare the original legacy force model. Flexor n=2 validation tests; extensor n=3. All six used extensor tests improve in RMSE, but the 43.2 cm tendon training test retains fraction of variance unexplained above one. Avoid stating every fit is uniformly good.');
table(s,[['Root-mean-square error (N·m)','No corrections','Corrected'],['Flexor, 2 validation tests','7.588','1.490'],['Extensor, 3 validation tests','2.554','1.098']],220,260);

await fig('The earlier single-actuator flexor falls short of its human target','Joint torque: reason to redesign','Revision C · Ch. 4, Fig. FlxGrp, panel B',
 '1 minute. This is the earlier 41.5 cm single 20 mm actuator configuration, not the later 56.2 cm pair. The maximum-isometric human reference is model-derived, not a person tested on this rig. It falls short at roughly 60–100 degrees of flexion. Lower-pressure points expanded accessible motion because full-pressure loading threatened the printed femur. This motivates route redesign rather than presenting the earlier capacity as the final result.','oldflex','Earlier geometry: one 20 mm actuator');

await pair('The optimized cam route became a two-actuator test assembly','Route redesign: geometry and hardware','Revision C · Ch. 4, Cam route, panel A; assembled hardware photograph',
 '1 minute. Route optimization changes attachment geometry under the adopted correction terms and geometric constraints. Credit Morrow’s collaboration in optimization. The cam extensor design of record uses two 20 mm BPAs. The assembly photograph establishes construction, while subsequent torque measurements establish its tested capacity. The image is a crop of existing panel A at −120 degrees; no geometry is redrawn.','camroute','hardware',['Route at −120°; coordinates in m','Assembled cam-extensor hardware']);

await fig('The cam model predicts at least 30.83% margin over the required envelope','Route redesign: model prediction','Revision C · Ch. 4, redesigned-route results; cam result of record, 14 October 2026',
 '1 minute. The required envelope is the selected vastus-medialis maximum-isometric target plus 5 percent. Minimum predicted margin: +30.83 percent over the required envelope, or +39.77 percent over the plain target. These are model predictions, not measured margins. Pressure needed to meet the required envelope is predicted at 269–507 kPa. This existing MATLAB plot uses the optimized 620 kPa curve. Distinguish it from the variable-pressure holds next.','camtorque','Prediction at 620 kPa; required envelope includes 5% overshoot');

await fig('Nine cam-extensor holds test torque across −119° to +5°','Route redesign: physical evidence','Revision C · Ch. 4, Fig. extensor_physical_final; MATLAB Ext20mm_compare_20261022.png',
 '1.5 minutes. Blue measured points, yellow pressure-specific model predictions, pink optimized 620 kPa curve, dotted human target. Marker area indicates pressure. The test set has nine points spanning −119 to +5 degrees. Measurements meet or exceed the human-reference magnitude over most, not all, tested angles. Pressure varies between holds; this does not measure a 620 kPa envelope at every angle. Use the final MATLAB reduction in Revision C, not the earlier paper paragraph that still says the reduction is pending, and never substitute Excel sanity-check cells.','extphysical','Pressure varied between holds');

await fig('All 15 flexor holds meet or exceed the human muscle target','Route redesign: physical evidence','Knee-torque study · Fig. Flx20mm; final MATLAB reduction',
 '1.5 minutes. Final 23 October MATLAB reduction with measured reaction-point geometry. Two 20 mm actuators, resting length 56.2 cm. Fifteen holds from −119 to +5 degrees. Magnitude margins over the model-derived biceps femoris short-head reference range from +0.42 N·m at −105 degrees to +17.0 N·m at −14 degrees. These are discrete isometric points, not a continuous experimentally verified envelope. Recorded pressures span 226.5–638 kPa, including points above the nominal 620 kPa model curve. Model peaks and measured points must not be compared as if pressures were matched. Plot regenerated from the original numeric MATLAB source with unchanged calculation, enlarged labels, and external legend.','meets_final','15 holds at 226.5–638 kPa; margins +0.42 to +17.0 N·m');

await fig('The corrected 620 kPa prediction still exceeds many measured flexor torques','Route redesign: remaining uncertainty','Knee-torque paper · Fig. Flx20mm, panel B; unchanged MATLAB calculation',
 '1 minute. Rigid prediction peak magnitude 190.3 N·m and corrected prediction 91.9 N·m are at 620 kPa. Measurements use different pressures, so this plot establishes neither a matched-pressure residual nor universal accuracy. The corrections substantially change the design prediction, but remaining mismatch motivates controlled matched-pressure tests and better transmission modeling. Distinguish the demonstrated capacity from validated accuracy of the full predicted curve.','rigid_final','Fixed-pressure curves; variable-pressure measurements');

await fig('SADb links sensory-feedback literature to inspectable source records','Sensory Afferent Database: gap and contribution','Revision C · Ch. 5, The Resource; SADb application snapshot',
 '1 minute. Gap: controller-design evidence is dispersed across species, preparations, tasks, and papers. I curate linked paper, pathway, model, review, and recorder records and build the read-only browser application. Latest distributed snapshot: 1,393 total records, including 943 papers and 450 other linked records. This is a curated research collection, not a systematic review. Offline embedded data support searching, pivoting, and citation exploration; the separate online-search tab requires connectivity. The existing screenshot illustrates the application, not a benchmark of retrieval accuracy.','clusters');
await fig('SADb exposes classification coverage without supplying controller gains','Sensory Afferent Database: use and limits','Revision C · Ch. 5, Fig. sadb_pivot and Use in Neural-Controller Design',
 '1 minute. A pivot cell counts tagged papers, not neurons or independent physiological observations. A paper can appear in more than one cell. Selecting a cell retrieves the underlying records. Afferent classes Ia, Ib, and II denote spindle-primary, tendon-organ, and spindle-secondary channels. Their properties depend on preparation and task. Source conditions, sign, and phase guide interpretation; translating them into robot encoders and gains remains an explicit design and verification step. AI-assisted drafting does not count as scientific curation without source checking. The corpus is uneven and not exhaustively classified.','pivot','Tagged papers retain their experimental context');

await pair('The initial AnimatLab controller used one rhythm generator for suspended stepping','Neural simulation: motivating limitation','Revision C · Ch. 6, The AnimatLab Bipedal Walker; Li (2017), Deng (2019), Bolen (2021)',
 '1 minute. Morrow and I combined Li’s mechanics with Deng’s two-level controller. The initial controller had ONE rhythm generator controlling both sides and a separate pattern-formation layer. I added a passive toe joint. Suspended alternating activity and stepping-like movement preceded subsequent feedback and bilateral-RG extensions. Hand tuning did not achieve weight-bearing walking. Do not attribute this earlier result to the later two-RG network. The body view and left-joint traces are the historical AnimatLab reference images. The angle units are not explicit in the source capture, so no numeric amplitude is asserted or converted.','historicalbody','historicaltrace',['Historical mechanical model','Original left-joint traces']);

s=base('I later interconnected separate left and right rhythm generators','Neural simulation: distinct contribution','Revision C · Ch. 6, The AnimatLab Bipedal Walker',
 '45 seconds. RG means rhythm generator. This editable conceptual schematic summarizes the documented bilateral architecture, not the full connectome. I constructed right-side subsystems from the left organization and connected the two RGs through commissural pathways. No walking result follows from assembly alone. The initial suspended stepping preceded this extension. Published Deng architecture underlies the controller; the bilateral extension is my contribution.');
function node(label,x,y,w=420){const a=s.shapes.add({geometry:'rect',position:{left:x,top:y,width:w,height:90},fill:'#FFFFFF',line:{fill:C.teal,width:2}});a.text=label;a.text.style={typeface:'Calibri',fontSize:29,color:C.text,bold:true,alignment:'center'};return a;}
const l=node('Left rhythm generator',90,210),r=node('Right rhythm generator',770,210);
const lp=node('Left pattern formation',90,425),rp=node('Right pattern formation',770,425);
s.shapes.connect(l,r,{line:{fill:C.teal,width:3},head:{type:'arrow'},tail:{type:'arrow'}});
s.shapes.connect(l,lp,{line:{fill:C.teal,width:3},tail:{type:'arrow'}});
s.shapes.connect(r,rp,{line:{fill:C.teal,width:3},tail:{type:'arrow'}});
text(s,'Commissural pathways',470,325,350,45,26,C.teal,true);
text(s,'Conceptual architecture; walking remains unverified',80,590,1120,50,29,C.teal,true);

await pair('My interfaces connect existing mechanics and neural formulations','Neural simulation: tool contribution','Revision C · Ch. 6, Model Conversion and Graphical Circuit Construction',
 '1 minute. OpenSim and MyoConverter supply geometry and initial conversion; MuJoCo supplies physics; SNS-Toolbox supplies conductance-based neural equations. My additions repair and audit custom muscle paths and joint signs, couple the quasi-static BPA-force model, and implement editable Simulink blocks and circuit interfaces. Neural and mechanical models exchange signals every 2 ms. There is no new neuron formulation or demonstrated integrated valve dynamics. Mechanical imports remain distinct from standalone reflex demonstrations. No Hunt/Young library comparison or speed superiority is claimed.','robot','library',['Repaired and audited model coupling','Editable implementations of existing equations']);

await fig('The graphical library implements defined Ia and Ib reflex pathways','Neural simulation: representative circuit','Revision C · Ch. 6, Fig. simulink_reflex_topology',
 '1 minute. Reduced-order antagonist knee model. Ia spindle excitation drives its homonymous motor pool and reciprocal inhibition of the antagonist. Ib tendon-organ input inhibits the homonymous pool in this demonstration. These are chosen circuit rules, not a universal claim about Ib function. Mechanical parameters and generic BPA-block parameters remain provisional. In the separate beverage-arm demonstration, reflex enabled sag was 2.00 degrees versus 7.85 degrees disabled. That reduced model does not establish reflex control of imported full-body mechanisms.','reflex','Spindle (Ia) and tendon-organ (Ib) inputs in a provisional model');

await fig('Matched rhythm-core implementations reproduce the same driven bursts','Neural simulation: matched verification','Revision C · Ch. 6, Fig. w2l_crossengine; matched NumPy and Simulink protocol',
 '1.5 minutes. The contact-driven W2L rhythm-core test receives matched inputs. Ten expected bursts appear, burst period 1.000 s, full-trace correlation at least 0.99996, with a residual 2 ms lag of one solver step. Removing phasic heel-contact inputs makes both realizations latch. This agreement is scoped to the circuit and protocol. It is neither autonomous oscillation proof nor walking verification. The separate 92-actuator tonic-drive test has a discrepancy and follows next.','traces','1.000 s periods; trace correlation ≥ 0.99996');

s=base('The bridge exchanges signals while full-network dynamics still disagree','Neural simulation: verification boundary','Revision C · Ch. 6, Closed-Loop Model Coupling and capability table',
 '1 minute. The bridge was checked against independent Python integration with single-muscle activation, different control rates, and sensor readback. Readback is bit-exact over the scored window. The exported neural network drove all 92 actuators for 5 seconds. However, tonic drive caused the Simulink realization to settle while the NumPy reference oscillated. Transfer of tuned controller parameters remains unresolved. Do not generalize the rhythm-core agreement to the full network. Hardware deployment and unsupported walking remain future work.');
table(s,[['Verification level','Observed result'],['Sensor and command exchange','Bit-exact scored readback'],['92-actuator closed-loop run','Signals exchanged for 5 seconds'],['Full network under tonic drive','Simulink settles; NumPy oscillates']],195,340);
text(s,'Interface correctness and controller equivalence need separate tests',55,590,1170,55,27,C.teal,true);

s=base('The evidence supports static torque capacity and specific simulation interfaces','Conclusions','Revision C · Chs. 3–7; final flexor MATLAB reduction, 23 October 2026',
 '1 minute. Return to the four original gaps. Force: supported by characterized lengths and isometric data. Torque: identified corrections improve held-out predictions and redesigned routes have physical tests. Sensory evidence: searchable, linked, partially curated literature resource. Neural tools: matched circuit and interface checks, with explicit full-network discrepancy. Neither static isometric capacity nor supported stepping demonstrates a walking robot. No harness-free walking exists in the program records.');
table(s,[['Contribution','Evidence','Boundary'],['Actuator force model','Measured force and held-out fits','Static, tested actuator range'],['Joint model and routes','Held-out error and physical holds','Variable pressure; residual mismatch'],['SADb','Inspectable linked records','Uneven source coverage'],['Neural interfaces','Matched tests and signal readback','Full-network discrepancy; no free walking']],170,410);

await fig('Matched-pressure tests and controller reconciliation are the next checks','Next experiments','Revision C · Ch. 7; Ch. 6, remaining limits',
 '45 seconds. The next experiments address the uncertainty exposed by the results: matched-pressure torque holds over the range, and diagnosing numerical/network differences before transferring tuned control. The dashed hardware path in this existing pipeline is proposed. Do not present a Jetson deployment, completed valve dynamics, or whole-body walking as achieved. Physical pulley tests and simultaneous biarticular verification are future work.','pipeline','Close the measured gaps before transferring control to hardware');

s=base('Questions','Acknowledgments','Ben P. Bolen · Mechanical Engineering · Portland State University',
 'Thank the advisor M. Anthony Hunt, the committee, the AARL, and the dissertation collaborators. Recognize M. Elzein, L. Pang, L. Burgess, C. Morrow, J. Scharzenberger, C. Steele, and J. Chung. This slide ends the main presentation. Hidden backup follows.');
text(s,'Thank you',65,240,1130,90,64,C.title,true);
text(s,'Advisor: M. Anthony Hunt\nCommittee, AARL, and research collaborators',65,370,1120,120,30);

await fig('The force model combines normalized force with resting-length scaling','Backup: actuator equation','Revision C · Ch. 3, force equation; Bolen et al. (2026)',
 'Define epsilon as contraction, P as pressure, l0 as resting length. Diameter-specific coefficients enter normalized force and maximum-force scaling. Maximum contraction is measured for each actuator. This model is quasi-static; neither valve dynamics nor hysteresis are included.','forceeq','Separate coefficients apply to each diameter',true);
await fig('Earlier biomimetic-extensor measurements retain residual prediction error','Backup: transfer to biomimetic knee','Revision C · Ch. 4, Fig. Ext10mm_52cm; see adopted shared-pair table for final statistics',
 'This is the earlier 10 mm, 51.8 cm configuration, not the new two-20 mm cam route. Do not mix earlier standalone fits with final shared-pair statistics. Revision C prose retains conflicting earlier 0.8920 RMSE and later adopted-pair 1.17 RMSE descriptions. The latter is explicitly the shared-pair comparison. Resolve the analysis scope before quoting either as a universal final fit.','oldext','Earlier 10 mm configuration',true);
await fig('Pulley transmission predictions still require physical verification','Backup: route extension','Revision C · Ch. 4, Mono- and Bi-Pulley Tendon Routing',
 'These surfaces are optimizer predictions, not physical measurements. Parallel actuator number and tackle line parts alter force/travel tradeoffs. Biarticular evaluation fixes the adjacent joint, so it does not establish simultaneous performance across both joints. Different route objective values are not directly comparable. Ongoing transmission studies are not finished hardware results.','pulley','Model predictions; adjacent joint fixed in biarticular planes',true);
await fig('Supported actuator-driven motion provides a coupling test','Backup: simulation scope','Revision C · Ch. 6 and Appendix, Simulation Verification',
 'Historical representative BPA-actuated run, retained-support/air-stepping scope. Neural output drives actuator forces and mechanical signals return to the controller. No harness-free walking. Do not use frozen or NaN runs as walking, and do not infer current tuned performance from a legacy score. Later wean studies reduce tonic drive, not mechanical support.','supported','Mechanical support remains',true);
await fig('Spiking and non-spiking realizations can produce different timing','Backup: neural coding','Revision C · Appendix, Simulation Verification; knee spiking-twin figure',
 'Representative reduced knee demonstration. Continuous motor activation interface is retained while the neural coding changes. This comparison is specific to selected parameters. It does not establish superiority, biological fidelity, full-network equivalence, or successful parameter transfer. Generic mechanical parameters remain provisional.','spiking','A matched interface does not ensure matched dynamics',true);
await fig('The human comparison uses one selected muscle reference','Backup: target definition','Knee-torque paper · Fig. Flx20mm, panel A; OpenSim Gait2392',
 'Maximum-isometric biceps femoris short-head target, negative flexion torque. This is one modeled monoarticular muscle, not total voluntary human knee torque, a patient outcome, or an experimentally measured person. Extensor optimization uses its separately defined selected vastus reference.','human_final','Biceps femoris short head; maximum isometric model torque',true);

for (let i=0;i<specs.length;i++){const refs=specs[i].assets.map(k=>manifest[k]?.source??path.join(root,'assets',k+'.png'));p.slides.items[i].speakerNotes.text += '\n\nFigure source paths:\n'+refs.join('\n');}
await fs.writeFile(path.join(root,'slide_plan.json'),JSON.stringify(specs,null,2));
await fs.mkdir(path.join(dir,'render'),{recursive:true});
const draft=path.join(dir,'candidate.pptx');
await (await PresentationFile.exportPptx(p)).save(draft);
console.log('Exported',specs.length,'slides');
const SKILL='C:/Users/Ben Bolen/.codex/plugins/cache/openai-primary-runtime/presentations/26.1020.11127/skills/presentations';
const {finalizePresentation}=await import(pathToFileURL(path.join(SKILL,'container_tools/artifact_tool_utils.mjs')).href);
const result=await finalizePresentation({workspaceDir:root,candidatePath:draft,finalPath:path.join(root,'output/Bolen_Defense_Assertion_Evidence_20261023_v6.pptx'),pythonExecutable:'C:/Users/Ben Bolen/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe',integrityValidatorPath:path.join(SKILL,'container_tools/inspect_presentation_package_integrity.py'),layoutValidatorPath:path.join(SKILL,'container_tools/inspect_presentation_layout_geometry.py'),layoutArgs:['--expected-slide-size-emu','12192000,6858000'],requiredNativeTableOwnerSlides:[],requiredNativeChartOwnerSlides:[],fontPolicy:{basis:'design',families:['Calibri']},materializeLiteralChartWorkbooks:true,verifyArtifactToolImport:true,receiptPath:path.join(dir,'validation_v6.json')});
console.log(JSON.stringify(result));








