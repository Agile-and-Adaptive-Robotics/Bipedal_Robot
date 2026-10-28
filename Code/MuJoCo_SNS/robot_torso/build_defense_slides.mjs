import fs from 'node:fs/promises';
import path from 'node:path';
import crypto from 'node:crypto';
import {pathToFileURL} from 'node:url';
const runtime='C:/Users/Ben Bolen/.cache/codex-runtimes/codex-primary-runtime/dependencies';
process.env.RUNTIME_NODE_MODULES=runtime+'/node/node_modules';
const skill='C:/Users/Ben Bolen/.codex/plugins/cache/openai-primary-runtime/presentations/26.1017.22136/skills/presentations';
const {Presentation,PresentationFile}=await import(pathToFileURL(runtime+'/node/node_modules/@oai/artifact-tool/dist/artifact_tool.mjs').href);
const {finalizePresentation}=await import(pathToFileURL(skill+'/container_tools/artifact_tool_utils.mjs').href);
const repo='D:/GitHub/Bipedal_Robot';
const workspace=repo+'/Documentation/Reports and Papers/Dissertation/Defense/slideshow/AH001_20261028';
const reference=repo+'/Documentation/Reports and Papers/Dissertation/Defense/AE_Revision_20261023/output/Bolen_Defense_Assertion_Evidence_20261023_v6.pptx';
const robot=repo+'/Code/MuJoCo_SNS/robot_torso';
const on=JSON.parse(await fs.readFile(robot+'/runs/abd_defense_exact_abduction_assessment.json','utf8'));
const off=JSON.parse(await fs.readFile(robot+'/runs/abd_defense_control_off_abduction_assessment.json','utf8'));
await fs.mkdir(workspace+'/.build',{recursive:true});await fs.mkdir(workspace+'/presentation',{recursive:true});
const presentation=Presentation.create({slideSize:{width:1280,height:720}});
presentation.theme.defaultFont='Calibri';
const brown='#B47C50',teal='#008579',ink='#333333';
function text(slide,content,x,y,w,h,size=24,color=ink,bold=false){
 const item=slide.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
 item.text=content;item.text.style={typeface:'Calibri',fontSize:size,color,bold};return item;
}
function shell(title,section,number){
 const slide=presentation.slides.add();slide.background.fill='#FFFFFF';
 text(slide,section,32,20,1216,28,18,teal,true);text(slide,title,32,58,1216,96,42,brown,true);
 text(slide,String(number),1228,687,28,24,18,'#666666');return slide;
}
const notes='Recorded physics from '+robot+'/runs/abd_defense_exact.npz. Model '+robot+'/ah001_abduction_diagnostic_v2.xml. 75 kg; 54 reference-tension actuators, 262-neuron hybrid SNS/postural controller, 2 ms step. 80% actual vertical support after 15 s, with lateral and angular pelvis assistance. Provisional relative-knee passive stops and opposite-leg CAD contacts are active. Abductor control combines separation feedback and a 6% allocation floor on currently abducting routes; global activation floor is 5%. This is not independent balance, calibrated BPA force validation, or accepted walking. Foot toe-out remains. Detailed comparison is in the source assessment JSONs. Fresh reference and accelerated runs are bit-identical.';
const slide1=shell('Abductor-enhanced control keeps the supported legs separated','CAD-DERIVED MODEL: SUPPORTED CONTROL',1);
slide1.images.add({blob:new Uint8Array(await fs.readFile(workspace+'/ah001_abductor_control.gif')),contentType:'image/gif',alt:'Recorded anterior and oblique views of the supported CAD robot with separated legs',fit:'cover',crop:{left:0,top:.0625,right:0,bottom:.0625},position:{left:208,top:164,width:864,height:504}});
text(slide1,'80% vertical support · hybrid control · reference tensions · knee/contact assumptions',32,684,1180,28,18,'#555555');
slide1.speakerNotes.text=notes+' The animation shows 20–32 s of the 40 s run. Right is blue, left gold, and selected hip-abducting routes are green.';
const slide2=shell('The matched control comparison quantifies leg separation','CAD-DERIVED MODEL: MEASURED COMPARISON',2);
slide2.images.add({blob:new Uint8Array(await fs.readFile(workspace+'/ah001_abductor_control_still.png')),contentType:'image/png',alt:'Same recorded robot state at 24 seconds, viewed from anterior and oblique cameras',fit:'contain',position:{left:32,top:170,width:744,height:434}});
const values=[['15–40 s interval','Off','On'],['Ankle order reversed',(100*off.reversed_ankle_order_fraction).toFixed(0)+'%',(100*on.reversed_ankle_order_fraction).toFixed(0)+'%'],['Minimum spacing',Math.round(1000*off.minimum_ankle_spacing_m)+' mm',Math.round(1000*on.minimum_ankle_spacing_m)+' mm'],['Right TFL activation',Math.round(100*off.mean_abductor_activation.tfl_r)+'%',Math.round(100*on.mean_abductor_activation.tfl_r)+'%'],['Vertical support','80%','80%'],['Abductor floor','5%','6%']];
const table=slide2.tables.add({rows:values.length,columns:3,left:800,top:180,width:448,height:350,columnWidths:[240,104,104],values});
table.borders.assign({style:'solid',fill:'#AAAAAA',width:1});
for(let r=0;r<values.length;r++)for(let c=0;c<3;c++){
 const cell=table.getCell(r,c);cell.fill='#FFFFFF';cell.text.style={typeface:'Calibri',fontSize:19,color:r===0?teal:ink,bold:r===0};
}
text(slide2,'All other conditions matched. Lateral and angular pelvis support remain active.',800,550,448,80,21,ink);
text(slide2,'This supported diagnostic does not establish independent balance or BPA-controlled walking.',32,646,1216,35,21,teal,true);
slide2.speakerNotes.text=notes+' The still is from 24.00 s. The ablation removes the extra separation gain and reduces the abductor floor to the 5% global floor. Other controls, geometry, constraints, support protocol, and initial pose are held fixed. Ankle order alone is not a full collision test; the mesh-contact and loop-closure audit is separate.';
const candidate=workspace+'/.build/candidate.pptx';await(await PresentationFile.exportPptx(presentation)).save(candidate);
const finalPath=workspace+'/presentation/AH001_supported_control_20261028.pptx';
const result=await finalizePresentation({workspaceDir:workspace,candidatePath:candidate,finalPath,
 explicitTotalSlideCount:2,requiredNativeTableOwnerSlides:[2],
 pythonExecutable:runtime+'/python/python.exe',
 integrityValidatorPath:skill+'/container_tools/inspect_presentation_package_integrity.py',
 layoutValidatorPath:skill+'/container_tools/inspect_presentation_layout_geometry.py',
 layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-heading-fit','--require-native-table-slide','2'],
 fontPolicy:{basis:'reference',families:['Calibri'],referencePath:reference,referenceSha256:crypto.createHash('sha256').update(await fs.readFile(reference)).digest('hex')},
 verifyArtifactToolImport:true,receiptPath:workspace+'/.build/presentation_validation.json'});
console.log(JSON.stringify(result));
