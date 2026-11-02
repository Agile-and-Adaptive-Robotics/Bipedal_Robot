import fs from 'node:fs';
import path from 'node:path';

// Seven question-slide repairs owned by the independent presentation supervisor.
export async function fixSupervisor(deck, workdir) {
  const snapshot = await deck.inspect({select:'$.slides[*]',include:['id','index','elements.id','elements.kind','elements.text','elements.position']});
  const records = new Map(snapshot.records.map(r=>[r.index+1,r]));
  const GREEN='#174D24', INK='#202A25';
  const asset = name => path.join(workdir,'assets',name);
  const slide = n => deck.slides.getItem(n-1);
  function clear(n) {
    const ids=records.get(n).elements.filter(e=>e.position && e.position.top>=88 && e.position.top<478).map(e=>e.id);
    if(ids.length)deck.delete(ids);
  }
  function tx(n,t,x,y,w,h,size=25,bold=false,color=INK,align='left') {
    const s=slide(n).shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
    s.text=t;s.text.insets=0;s.text.style={fontFamily:'Arial',fontSize:size,bold,color,alignment:align};return s;
  }
  function im(n,file,x,y,w,h,crop) {
    slide(n).images.add({blob:fs.readFileSync(file),contentType:'image/png',alt:path.basename(file),fit:'contain',position:{left:x,top:y,width:w,height:h},...(crop?{crop}:{})});
  }
  function title(n,t) {
    const e=records.get(n).elements.find(e=>e.kind==='textbox' && e.position?.top<80 && e.position?.width>650);
    if(e) { const s=deck.resolve(e.id);s.text=t;s.text.style={fontFamily:'Arial',fontSize:31,bold:true,color:GREEN}; }
  }
  function cite(n,t) {
    const e=records.get(n).elements.find(e=>e.kind==='textbox' && e.position?.top>=478 && e.position?.top<510 && e.position?.width>700);
    if(e)deck.resolve(e.id).text=t;
  }
  function append(n,t){slide(n).speakerNotes.text=String(slide(n).speakerNotes.text||'')+'\n\n'+t;}

  // Both axes, legend, and the named supported condition remain explicit.
  clear(76);
  im(76,asset('s076_img0.png'),30,125,900,290);
  tx(76,'Gait cycle (%)',345,434,300,38,25,false,INK,'center');
  cite(76,'Dissertation Appendix D · Supported human-model comparison with OpenSim IK');
  append(76,'The horizontal coordinate is normalized gait cycle in percent. Blue shading is cycle-to-cycle variation in the model; the dashed trace is the OpenSim inverse-kinematics reference. Both plots show hip angles, not forward progression.');

  clear(77);
  tx(77,'Solid blue: model     Dashed black: OpenSim IK',165,106,700,34,23,false,INK,'center');
  tx(77,'Right',150,157,300,35,25,true,GREEN,'center');
  tx(77,'Left',585,157,300,35,25,true,GREEN,'center');
  im(77,asset('s077_img0.png'),30,190,900,225);
  tx(77,'Duty: 0.71 left · 0.65 right · reference 0.61',45,448,875,36,24,true,GREEN);
  cite(77,'Dissertation Appendix D · Supported human-model comparison with OpenSim IK');
  append(77,'Panel order follows the complete source figure: right, then left. Stance fraction is the fraction of aligned cycles in contact at each phase, rather than a single-cycle binary contact signal.');

  // Use the same source, preserving plot colors. Compare one close and one poor match.
  clear(78);
  title(78,'Matched-input force agreement varies across the converted muscles');
  im(78,asset('s078_img0.png'),70,112,820,45,{left:.20,top:.046,right:.20,bottom:.918});
  tx(78,'Vastus lateralis: R² = 0.94',45,177,430,36,24,true,GREEN,'center');
  tx(78,'Semimembranosus: R² = 0.57',485,177,440,36,24,true,GREEN,'center');
  im(78,asset('s078_img0.png'),35,215,890,190,{left:0,top:.405,right:0,bottom:.35});
  // The two source panels share the bottom-row x scale. Native labels avoid
  // repeating that row's baseline and zero-y ticks in a separate image strip.
  for(const start of [122.8,554.2]) for(let k=0;k<6;k++)
    tx(78,String(k*20),start+k*67.3-17,407,34,25,18,false,INK,'center');
  tx(78,'% of analyzed walk',145,435,315,28,20,true,INK,'center');
  tx(78,'% of analyzed walk',575,435,315,28,20,true,INK,'center');
  cite(78,'Dissertation Appendix D · Matched poses and OpenSim static-optimization activations');
  append(78,'Shown are two panels from the six-muscle comparison. Vastus lateralis has R² 0.94 and reported lag one frame; semimembranosus has R² 0.57 and lag zero frames. The source uses identical poses and OpenSim static-optimization activations for the comparison. The normalized traces do not establish equality of absolute force. The complete six-panel source remains in the reproducible assets.');

  clear(79);
  title(79,'The graphical editor exposes neural connection parameters');
  im(79,asset('s084_img0.png'),90,110,780,354);
  cite(79,'Dissertation Chapter 6 · SNS-Simscape graphical editor');
  append(79,'This is the implemented editor, showing the palette, example circuit, and selected connection properties. Neural parameters remain inspectable in the property panel before a Simulink model is built.');

  clear(82);
  title(82,'The BPA blocks evaluate the measured force model during simulation');
  im(82,asset('s082_img0.png'),70,115,820,64);
  im(82,asset('s082_img1.png'),355,195,570,275);
  im(82,asset('sns_library_grid_current.png'),35,225,275,180,{left:.255,top:.665,right:.525,bottom:.105});
  cite(82,'Dissertation Chapter 6 · SNS-Simscape BPA blocks; Bolen et al. (2026)');
  append(82,'The displayed BPA10mm block is from the actual library export. Diameter selects the fitted coefficients; current length, pressure, and initial length determine force using the same normalized static model as the actuator chapter. This static force block does not establish valve or pneumatic transient dynamics.');

  clear(83);
  title(83,'Afferent blocks convert muscle state into sensory input');
  im(83,asset('sns_library_grid_current.png'),70,160,355,240,{left:0,top:.38,right:.78,bottom:.38});
  im(83,asset('s083_img0.png'),535,160,355,240);
  tx(83,'Muscle length and velocity',45,425,420,40,25,true,GREEN,'center');
  tx(83,'Muscle force',510,425,420,40,25,true,GREEN,'center');
  cite(83,'Dissertation Chapter 6 · Implemented Ia spindle and Ib tendon-organ blocks');
  append(83,'The left image is the actual IaMuscleSpindle block and the right is IbGolgiTendon. These are sensory encoders: they convert the supplied mechanical quantities to neural inputs. Correct encoding does not itself establish the physiological targets of every downstream synapse.');

  clear(87);
  title(87,'Coupling a CAD model requires actuator and sensory interfaces');
  im(87,asset('s088_img0.png'),35,115,320,345,{left:0,top:.075,right:.5,bottom:0});
  tx(87,'Length and velocity',380,130,525,40,25,false,INK,'center');
  tx(87,'→',383,190,75,55,44,true,GREEN);
  tx(87,'Neural controller',500,200,380,45,29,true,GREEN,'center');
  tx(87,'↓',695,268,65,55,44,true,GREEN);
  tx(87,'Activation / pressure',500,333,395,40,25,false,INK,'center');
  tx(87,'Actuator force model',490,388,410,45,29,true,GREEN,'center');
  tx(87,'←',380,382,80,55,44,true,GREEN);
  cite(87,'Dissertation Chapter 6 · CAD model and neural–mechanical interfaces');
  append(87,'This is an interface diagram, not a claim that every CAD route is driven by a neural BPA. The submitted CAD demonstration used hybrid ideal tensions and support; the particular actuator implementation must be identified for each trial. Sensory values return from mechanics, while neural commands pass through the selected actuator model.');
  return [76,77,78,79,82,83,87];
}
