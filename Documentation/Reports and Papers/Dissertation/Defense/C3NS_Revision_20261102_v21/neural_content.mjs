import path from 'node:path';
import fs from 'node:fs';

export async function apply(h) {
  const {d,rec,begin,newSlide,im,tx,line,table,changed,GREEN,INK,a}=h;
  const N=p=>a(path.join('neural',p)), T=p=>N(path.join('tools',p));
  const campaign='D:/Users/Ben Bolen/NeuroCampaign_20261102_9h';
  const source=campaign+'/RESULTS_AND_FILES.md';
  const scope='Selected human Gait2392 / graded biw2lvar-Hunt-PT with a four-cell RF/HAM PF module per leg and the separately enabled autogenic Renshaw repair. Ground9 RC S4 is guided stepping in place with five right/four left qualifying unloaded, whole-foot-mesh-clear swings during 5–13 s. Net vertical assistance is 48.848356% BW and horizontal restraint RMS is 115.560892 N. These measurements use the actual post-physics sample alignment and pelvis-subtree weight 768.756726 N. It is a separate post-submission result, not W2L108 or historical Syn6. ';
  const foot='Dissertation Ch. 6, Table 6.1';
  const citeNew='Gait2392 / graded Hunt–PT · selected ground9 RC · 2 Nov 2026';
  const backup=[];
  function box(n,text,x,y,w,hgt,size=24) {
    d.slides.getItem(n-1).shapes.add({geometry:'rect',position:{left:x,top:y,width:w,height:hgt},fill:'#F4F7F5',line:{fill:GREEN,width:1.5}});
    tx(n,text,x+12,y+10,w-24,hgt-16,size,true,INK,'center');
  }
  function note(n,extra) {d.slides.getItem(n-1).speakerNotes.text+= '\n'+extra;changed.push(n);}

  begin(57,'Detailed neuromechanical models require a faster, repeatable testing workflow',
    'The previous suspended-versus-dropped AnimatLab comparison motivates a workflow gap. AnimatLab allowed graphical construction and interactive inspection, but the transfer to more detailed mechanics required robust programmatic coupling, repeatable parameter and stimulus experiments, and faster execution. Its unsupported development status motivated moving the workflow, rather than implying that AnimatLab lacked graphical inspection. The figure is a conceptual controller/body signal loop using an authentic source neural schematic and the actual human reference model. The required measurements include neural activity, motion, loading, and all external assistance. The following tools address this gap in the order of dissertation Table 6.1. Sources: submitted Chapter 6, sections Knowledge Gaps and the New Tools; Hunt INRC2 lecture; SESSION1_NEURAL_CONTROLLER_HANDOFF.md.',
    foot+' · Knowledge gap after the suspended / ground comparison');
  im(57,T('hunt_neural_schematic.png'),40,140,365,285);
  im(57,T('gait2392_model.png'),665,123,210,300);
  tx(57,'Inspectable pathways',35,106,390,32,24,true,GREEN,'center');
  tx(57,'New mechanics',620,106,300,32,24,true,GREEN,'center');
  line(57,432,210,625,210,true);
  tx(57,'Motor commands',432,169,193,35,21,false,INK,'center');
  line(57,625,330,432,330,true);
  tx(57,'Sensory feedback',427,341,203,35,21,false,INK,'center');
  tx(57,'Controller fitting needed repeatable automated trials',70,435,840,42,25,true,GREEN,'center');

  begin(58,'The model library links simplified walkers with detailed human mechanics',
    'The first two rows of Table 6.1 supply complementary test bodies. Li–Deng/W2L is a simplified sagittal controller-transfer example. Gait2392 supplies detailed human musculoskeletal geometry; a separate robot counterpart uses BPA routes and the measured force law. Six floating-base coordinates and seventeen internal coordinates belong to the human model. Do not interpret a small sagittal controller as demonstrated control of every remaining degree of freedom. These source images illustrate model implementations, not one identical controller across bodies. Sources: submitted Chapter 6; MODEL_DIMENSION_AND_CONTROL.md.',
    foot+' · 1–2. Walker and Gait2392-based models');
  im(58,T('w2l_model.png'),35,140,420,270);
  im(58,T('gait2392_model.png'),600,115,240,310);
  tx(58,'Simplified W2L',40,431,410,36,25,true,GREEN,'center');
  tx(58,'Human musculoskeletal model',490,431,440,36,25,true,GREEN,'center');

  begin(59,'Neural blocks connect sensory input to actuator commands',
    'SNS-Simscape implements graphical neuron, synapse, sensory, activation, biological-muscle, and BPA-force components. This authentic library figure comes from the submitted component PDF. These blocks support circuit inspection and mechanical coupling; block and reflex demonstrations validate declared component functions. They do not certify an entire walking controller. The loaded-arm response and matched implementation results are retained in questions. Sources: dissertation Table 6.1 and SNS-Simscape library description.',
    foot+' · 3. SNS-Simscape block library and editor');
  im(59,T('sns_library.png'),60,110,840,345);

  begin(60,'CAD assemblies bring robot geometry and actuator mechanics into the workflow',
    'The SOLIDWORKS-derived robot is a separate integration test. CAD geometry, joint definitions, inertias, and actuator routes are assembled for mechanical simulation; SNS-Simscape supplies neural and actuator components. The displayed source assembly predates the laptop 3-DoF-foot upgrade. Existing robot results retain this earlier model identity. The submitted supported CAD demonstration uses ideal reference tensions/hybrid actuation. A later pressure-pair example uses only two pressure-modeled BPAs and 52 reference-tension routes; neither is a completed 54-BPA neural walking controller. Sources: submitted Chapter 6; current campaign RESULTS_AND_FILES.md CAD media section.',
    foot+' · 4. SOLIDWORKS-derived robot assembly');
  im(60,T('cad_model.png'),40,110,605,330,{left:.5,top:.08,right:0,bottom:0});
  im(60,T('sns_library.png'),650,215,260,170,{left:0,top:.66,right:.5,bottom:.08});
  tx(60,'Geometry, joints, and muscle routes',35,440,610,35,24,true,GREEN,'center');
  tx(60,'BPA force blocks',650,402,260,40,23,true,GREEN,'center');

  begin(61,'The bridge exchanges muscle commands and sensory measurements',
    'The MuJoCo–Simulink interface connects block-programmed neural control to MuJoCo mechanics. Activation commands enter the mechanical model; sensor values return to the controller. The authentic readback plot is the submitted short interface check, with five channels, 80 scored steps, and zero maximum deviation in that test. The question table records 81 samples on the 5 ms grid; scored transitions and saved samples are distinct counts. This is interface verification, not a gait or timing benchmark. Source: submitted Chapter 6, MuJoCo–Simulink integration.',
    foot+' · 5. MuJoCo–Simulink interface');
  box(61,'SNS-Simscape\ncontroller',40,130,300,85,25);
  box(61,'MuJoCo\nmechanical model',40,340,300,85,25);
  line(61,132,232,132,322,true);
  line(61,256,322,256,232,true);
  tx(61,'Commands',37,264,95,44,20);
  tx(61,'Sensory\nmeasurements',264,253,124,64,20);
  im(61,T('bridge_readback.png'),405,135,520,290);
  tx(61,'Matched sensor traces verify the interface',405,435,520,35,23,true,GREEN,'center');

  begin(62,'Staged experiments separate rhythm, motor output, and local feedback',
    'Hunt’s four-part instructional method is functional calibration, not the physical curriculum. A: select a free-running CPG with approximately twice the target stepping period. B: replay expected afferents and fit timing at the target cadence. C: freeze timing and fit CPG/PF-to-MN output. D: freeze feedforward settings and fit local Ia/Ib feedback to refine the waveform. These operations are implemented through saved configurations, scripted inputs, parameter fitting, and recorded neural and mechanical responses. In the selected November campaign, intrinsic/timing fits were reused and checked, then four shared local feedback targets and three paired motor-output settings were refined. Physical air, standing, and ground stages remain separate tests. Citation: Hunt, A. J. (2019, February 15). Real-time control of dynamic walking in a dog-like robot. Presented at the 2nd Intel Neuromorphic Research Community Meeting. Diagram: source controller schematic retained from the submitted deck.',
    'A. J. Hunt, INRC2 Meeting, Intel Labs (Feb 15, 2019) · Table 6.1 experiment drivers');
  im(62,T('hunt_neural_schematic.png'),30,125,395,315);
  const steps=[['1','Intrinsic rhythm'],['2','Timing with expected feedback'],['3','Motor-output pattern'],['4','Local sensory refinement']];
  steps.forEach(([num,text],i)=>{tx(62,num,475,135+i*74,40,42,30,true,GREEN);tx(62,text,532,137+i*74,385,55,27,true,INK);});
  tx(62,'Each stage freezes the solution needed by the next',45,446,870,31,23,true,GREEN,'center');

  begin(63,'Literature circuits give controlled stimuli a defined neural target',
    'The source schematic is attributed to Deng et al. (2019). It establishes the hierarchy used to select a rhythm-generator intervention. The adjacent trace is the separate reduced Deng-derived circuit experiment, with a +50 nA RG pulse lasting 3.876 s. This is a neural mechanism test; it is not the public rat120 body controller or the selected human Hunt/PT result. Detailed motor deletions, perturbation layers, and spiking variants remain in questions. Literature models and variant implementations are tools for testing specified mechanisms, not interchangeable physiological ground truth. Sources: dissertation Chapter 6 and Appendix D; Deng et al. (2019); NON_GAIT_AUTHOR_CHOICE_RECORD.md.',
    foot+' · 7. Literature circuits and variants · Deng et al. (2019)');
  im(63,T('deng2019_source_circuit.png'),35,140,360,270,{left:0,top:0,right:.66,bottom:0});
  im(63,T('deng_rg_stimulus.png'),420,184,505,228,{left:0,top:0,right:.50,bottom:0});
  tx(63,'Specified intervention: RG +50 nA for 3.876 s',395,116,540,55,23,true,GREEN,'center');
  tx(63,'Time from stimulus onset (s)',475,423,430,32,22,false,INK,'center');
  tx(63,'Source architecture',35,432,350,34,24,true,GREEN,'center');

  begin(64,'Compiled execution reduced complete-process time by 6.82×',
    'This paired benchmark compares source_native_s1 with source_baseline_s1, not the final ground9 body trial. Whole-process times are 199.623573 and 29.284280 s, including fresh-process startup and capture/diagnostic costs. Accumulated neural calls are 179.079112 and 9.779516 s. Whole-process acceleration is 6.817×; neural calls accelerate 18.312×. Forty-five numeric fields across six saved files agree exactly, including recorded neural, mechanical, and force arrays. Current fully instrumented Gait trials are not universally faster than real time. The old suspended-W2L benchmark retains its own configuration in questions. Source: '+campaign+'/figures/BACKEND_AND_VALIDATOR_COMPARISON.md.',
    'Matched Gait source trial · full process and neural-call scopes reported separately');
  table(64,[['Runtime scope','Native (s)','Accelerated (s)','Speed-up'],['Complete process','199.624','29.284','6.82×'],['Neural updates','179.079','9.780','18.31×']],55,165,850,185,[250,190,230,180],25);
  tx(64,'Identical recorded trajectories in the paired test',65,395,830,55,29,true,GREEN,'center');

  begin(65,'RF and hamstring pathways combine inputs from two joints',
    scope+'The authentic executed right-leg subnetwork has four cells: an RF integrator, a HAM integrator, and their two inhibitory relays. Hip-flexor plus knee-extensor PF inputs converge on RF. Hip-extensor plus knee-flexor PF inputs converge on HAM. Each integrator excites its own relay, which inhibits the opposing integrator. RF drives rectus femoris; HAM drives semimembranosus, semitendinosus, and biceps femoris long head. Sartorius and gracilis retain separate modules. The displayed diagram preserves all twelve incident edges of these cells; full-precision source fields are in the adjacent ledger. This is an engineering adaptation rather than a literal reproduction of the rat biarticular circuit. Source: '+campaign+'/figures/rfham_selected_subnetwork/CAPTION.md. New EB475WS4 delivery is used only if its provenance is independently confirmed.',
    citeNew+' · Actual selected RF/HAM subnetwork');
  im(65,N('rfham_selected/RFHAM_FOURCELL_SELECTED.svg'),45,103,870,367);

  begin(66,'The selected circuit links rhythm to patterned motor output',
    scope+'The displayed right-leg RG, four-cell RF/HAM PF, and knee-MN traces are three selections from the same 7.663–13.000 s recording. The landscape derivative preserves all original data paths, actual graded-model rest-relative voltages, color/line styles, and y scales. Only axes placement and labels changed. The RF/HAM integration and relay cells are not necessarily equal in amplitude, and these traces are not EMG. Both legs and all hip/knee/ankle MN selections remain in the original source plate and supporting question material. Source: '+campaign+'/figures/ground9_rc_s4/plots/PRESENTATION_NEURAL.pdf.',
    citeNew+' · Recorded activity, right leg');
  tx(66,'RG: rhythm generator   PF: pattern formation   MN: motor neuron',45,112,870,33,21);
  im(66,N('slide_panels/NEURAL_RIGHT_HIERARCHY.svg'),25,165,910,292);

  begin(67,'The human model takes repeated steps with measured support',
    scope+'The GIF replays saved states at 1× over 7.663–13.000 s, approximately three cycles. It performs no new integration. The foot-load panel is an exact vector crop of the matching saved mechanics PDF, with force normalized by the same body-weight denominator. Assessment counts apply to 5–13 s, rather than only the shorter replay window. A minimum support threshold and forward progression were not established. Source: '+source+'; '+campaign+'/figures/ground9_rc_s4/plots/SUPPORT_DECOMPOSITION_5_TO_13.json.',
    citeNew+' · Guided stepping in place');
  im(67,N('ground9_rc_s4/PRESENTATION_REPLAY.gif'),25,150,510,295,{left:0,top:.1334,right:0,bottom:0});
  im(67,N('slide_panels/FOOT_LOADS_AND_ASSISTANCE.png'),555,135,380,295);
  tx(67,'5 right / 4 left clear swings',30,107,490,36,26,true,GREEN,'center');
  tx(67,'48.8% BW vertical assistance; 115.6 N horizontal guidance',40,440,880,35,24,true,GREEN,'center');

  begin(68,'Quarter-strength contact feedback restored clearance in this test',
    scope+'This is the S5 comparison, distinct from the S4 movie. The same RC-enabled controller and body are compared at contact-to-PF scale 0.25 and 1.0, with matching apparatus settings. The whole-foot mesh/unloading assessment finds five right/four left clear swings at quarter strength and two/two at original strength over 5–13 s. Actual right ankle PF channels and bilateral vertical foot forces are shown with their original amplitudes, time stamps, colors, and line styles. PF/DF in the plot means plantarflexion/dorsiflexion; R/L means right/left. Toe-to-dorsiflexion-PF inhibition is stronger in the original condition. The result motivates phase-dependent feedback testing; no such new policy is installed here. Source: '+campaign+'/figures/stage5_contact_strength_review/CAPTION.md.',
    'Selected Hunt–PT S5 · matched contact-to-PF strength comparison');
  im(68,N('slide_panels/CONTACT_PF_AND_FOOT_LOADING.svg'),110,102,740,333);
  tx(68,'5 right / 4 left clear swings',65,438,410,36,25,true,GREEN,'center');
  tx(68,'2 right / 2 left clear swings',495,438,410,36,25,true,GREEN,'center');

  // Keep the mechanical overview/future slides. Update the neural future work
  // that the newly completed campaign makes stale.
  begin(71,'Further tests must reduce support while preserving clear, alternating steps',
    scope+'The next lower-lift quarter-strength S5 condition failed the clearance test with only two qualifying swings per foot. A minimum support level has not been established. The proposed next experiments reduce prescribed lift and horizontal/rotational guidance separately while keeping explicit whole-foot clearance, unloading, knee safety, and repeatability criteria. This future work does not relabel guided stepping as independent balance. Source: current campaign RESULTS_AND_FILES.md and COMMITTEE_QA.md.',
    'Future work · Support reduction with explicit movement criteria');
  im(71,N('ground9_rc_s4/PRESENTATION_STILL.png'),35,130,520,315,{left:.5,top:.1334,right:0,bottom:0});
  tx(71,'Reduce lift',605,152,300,55,30,true,GREEN);
  tx(71,'Reduce guidance',605,248,320,55,30,true,GREEN);
  tx(71,'Retain clearance\nand unloading',605,344,325,95,28,true,GREEN);
  tx(71,'Test vertical support and lateral balance separately',45,447,870,33,24,true,GREEN,'center');

  begin(72,'Contact feedback can be tested at specific gait phases',
    'Proposed experiment informed by the selected Hunt/PT S5 comparison. First identify stance and swing from the declared foot-load/unloading rules. Compare unchanged always-enabled contact-PF feedback with phase-gated variants, keeping the remaining controller, physical model, and apparatus fixed. Record RG phase, dorsiflexor PF current, MN activity, whole-foot clearance, and loading. No phase-gated controller is implemented or dynamically validated by this slide. This replaces the older W2L-only future-test proposal so that the next experiment follows the new measured failure. Source: current campaign COMMITTEE_QA.md and stage5_contact_strength_review/CAPTION.md.',
    'Future experiment · Test the observed contact-feedback mechanism');
  tx(72,'Current result',60,122,360,38,27,true,GREEN);
  tx(72,'Full contact feedback\nsuppressed swing clearance',60,184,360,100,28);
  tx(72,'Proposed comparison',490,122,410,38,27,true,GREEN);
  table(72,[['Feedback timing','Measure'],['Always enabled','Reference response'],['Stance only','Support and clearance'],['Swing only','Dorsiflexor suppression']],490,182,410,215,[180,230],20);
  tx(72,'Hold mechanics, support, and the remaining controller fixed',55,444,850,38,24,true,GREEN,'center');

  begin(74,'This work advances actuator prediction, joint design, and neural testing',
    'The force model accounts for initial length, pressure, and contraction. Experimentally identified design factors improve joint-torque prediction and guide knee designs that meet the selected human targets at measured angles. The neural tools enable inspectable controller construction and repeated, measured body experiments. The third panel shows the selected post-submission guided Gait2392 result: five right/four left qualifying clear unloaded swings, 48.8% BW net vertical assistance, and 115.6 N horizontal guidance. It is not unsupported walking. The contact-strength comparison demonstrates a specific controller-body experiment enabled by the workflow. Mechanical calibration, functional circuit checks, and supported stepping supply complementary evidence. Hardware/dynamic validation and further support reduction remain future work. Sources: dissertation conclusion and '+source+'.',
    'Dissertation contributions · New neural demonstration from the November 2 campaign');
  im(74,T('conclusion_force.png'),30,145,280,235);
  im(74,T('conclusion_torque.png'),340,155,280,220);
  im(74,N('ground9_rc_s4/PRESENTATION_STILL.png'),650,132,280,254,{left:.5,top:.1334,right:0,bottom:0});
  tx(74,'Supported Gait2392',650,101,280,30,21,true,GREEN,'center');
  tx(74,'Predict\nBPA force',30,410,280,64,26,true,GREEN,'center');
  tx(74,'Meet selected\nknee-torque targets',340,410,280,64,26,true,GREEN,'center');
  tx(74,'Test neural control\nwith measured support',650,408,280,70,24,true,GREEN,'center');
  tx(74,'Thank you',390,480,180,30,25,true,GREEN,'center');

  let n=newSlide('neural_stage_evidence','The same selected controller passed the first four physical stages',
    scope+'The nominal RC-enabled saved controller passed fresh air stages 1–2, supported standing 3, and supported ground stepping 4. The quarter-strength S5 extension also passes, but its contact-PF scale differs from S4. Full-strength S5 and the next lower-lift rung fail clearance. Physical stage numbers do not identify Hunt A–D calibration operations. The result is a feasible tested candidate, not proof of global optimality or normal walking. Source: '+source+'; ALL_CONTROLLER_DISPOSITION_20261102.md.',
    'Selected Hunt–PT controller · Saved-condition validation');
  table(n,[['Physical condition','Result'],['1. Air, rhythm test','Pass'],['2. Air, sensory feedback','Pass'],['3. Supported standing','Pass'],['4. Supported ground steps','Pass, 5 right / 4 left'],['5. Contact-PF scale 0.25','Pass, 5 right / 4 left']],65,130,830,300,[520,310],22);
  tx(n,'Stage 5 changes contact-feedback strength',90,444,780,34,23,true,GREEN,'center');backup.push(n);

  n=newSlide('neural_support','Measured vertical assistance includes lift and restraint',
    scope+'Mean superior support is 40.051043% BW explicit lift plus 8.797312% BW signed vertical restraint, totaling 48.848356% BW. The selected denominator is the 78.3646 kg pelvis subtree, weight 768.756726 N. Horizontal restraint RMS is 115.560892 N (15.03% BW), separate from vertical assistance. The original source anatomical masses total 75.1646 kg; 32 numerical follower frames have 0.1 kg each. No mechanics were changed for this deck. Source: SUPPORT_DECOMPOSITION_5_TO_13.json and reproducibility/MASS_PROVENANCE.json.',
    citeNew+' · Actual aligned 5–13 s force record');
  table(n,[['Vertical contribution','% body weight'],['Prescribed lift','40.051'],['Signed restraint','8.797'],['Net assistance','48.848']],40,140,450,230,[270,180],22);
  im(n,N('slide_panels/FOOT_LOADS_AND_ASSISTANCE.png'),515,125,410,285);
  tx(n,'Horizontal restraint RMS: 115.6 N',60,433,840,40,27,true,GREEN,'center');backup.push(n);

  n=newSlide('neural_implementation_match','NumPy and Simulink reproduce the same driven W2L rhythm',
    'Preserved submitted result, separate from the new human ground9 study. With matched parameters and input signals, the selected driven W2L neural core has a 1.000 s period, ten expected bursts, and trace correlation at least 0.99996. This demonstrates agreement for that declared core/protocol, rather than arbitrary full-controller equivalence. Source: submitted Chapter 6.',
    'Submitted W2L neural-core comparison');
  im(n,T('w2l_numpy_simulink.png'),35,118,890,345);backup.push(n);

  n=newSlide('neural_reflex_response','Sensory feedback reduces forearm sag in the loaded-arm test',
    'Preserved submitted SNS-Simscape component demonstration. Under the specified cup-loading protocol, feedback reduces final forearm sag from 7.85 degrees to 2.00 degrees. This verifies a coupled sensory-to-motor response in the implemented blocks, rather than locomotor balance. Source: dissertation Chapter 6 loaded-arm reflex example.',
    'SNS-Simscape component verification · 7.85° to 2.00° final sag');
  im(n,T('loaded_arm_reflex.png'),35,150,350,245,{left:.50,top:.04,right:0,bottom:.75});
  im(n,T('loaded_arm_reflex.png'),445,125,480,320,{left:.25,top:.325,right:.27,bottom:.35});backup.push(n);

  n=newSlide('neural_motor_maps','The circuit keeps motor conversion and force feedback explicit',
    'New EB475WS4 drawing, local RF crop of the named historical w2lvar standing export. This is not the current ground9 Hunt/PT controller. The source ledger SHA256 is 49b946c627222e030aa7bc741d7de20faeb804646d791b51548e6d3d704d4e2f and was checked against the local source. MN voltage passes through an adapter/map to muscle activation. Muscle force passes through a separate adapter to Ib and the local inhibitory circuit; Renshaw feedback is also shown. AD is a map, not another neuron. Connections continuing above/right leave this local viewport; the full paired RF/semimembranosus drawing is preserved beside the source SVG. Gray dashed descending PF-to-MN routes have zero effective conductance in this saved standing condition. Original SVG shapes and paths are unchanged. Source: functional_subnetwork_review/drawings/full_motor_loops/README.md and VERIFICATION.json.',
    'EB475WS4 schematic · historical w2lvar standing export · RF local crop');
  im(n,N('slide_panels/RF_MOTOR_FORCE_CROP.svg'),55,140,850,285);
  tx(n,'AD: adapter/map     RF: rectus femoris     MN: motor neuron',55,439,850,35,23,true,GREEN,'center');backup.push(n);

  n=newSlide('neural_sensory_maps','Separate sensory maps preserve Ia and II input identities',
    'New EB475WS4 drawing, local RF sensory crop of the historical w2lvar standing export, not the current ground9 controller. The spindle symbol identifies encoded muscle length and velocity, without adding new mechanistic receptor dynamics. Separate adapter/maps inject the saved Ia and II currents; Ia, II excitatory, and II inhibitory branches remain distinct. Brown paths denote physical signals/current injection and black paths denote neural connections. Connections continuing below/right lead to motor and antagonist pathways outside this viewport. Source body and complete edge/physical-map tables are preserved in assets/neural/eb475ws4_motor_loops. Source ledger hash 49b946c627222e030aa7bc741d7de20faeb804646d791b51548e6d3d704d4e2f.',
    'EB475WS4 schematic · historical w2lvar standing export · RF sensory crop');
  im(n,N('slide_panels/RF_LENGTH_SENSORY_CROP.svg'),55,125,850,315);
  tx(n,'Brown: physical signals    Black: neural connections',65,445,830,34,23,true,GREEN,'center');backup.push(n);

  return {main:Array.from({length:12},(_,i)=>57+i),backup};
}
