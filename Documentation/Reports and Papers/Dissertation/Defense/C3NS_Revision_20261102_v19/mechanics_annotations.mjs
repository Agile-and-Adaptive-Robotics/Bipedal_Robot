import fs from 'node:fs';
import path from 'node:path';

// Source definitions: 20-muscles.tex 20,78; 30-torque.tex 93–100,661–662.
// Owns13,32–40 except untouched slide30. Distances on photos exclude fittings.
export async function fixMechanicsAnnotations(deck, workdir=process.cwd()) {
 const GREEN='#174D24', INK='#202A25', GRAY='#75827B', YELLOW='#FFE399', MAGENTA='#DD00DD';
 const snap=await deck.inspect({select:'$.slides[*]',include:['id','index','elements.id','elements.kind','elements.position']});
 const records=new Map(snap.records.map(r=>[r.index+1,r]));
 function wipe(n){
  const ids=records.get(n).elements.filter(e=>e.position&&e.position.top<478).map(e=>e.id);
  if(ids.length)deck.delete(ids);
 }
 function tx(n,text,x,y,w,h,size=25,bold=false,color=INK){
  const s=deck.slides.getItem(n-1).shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
  s.text=text;s.text.insets=0;s.text.style={fontFamily:'Arial',fontSize:size,bold,color,alignment:'left'};return s;
 }
 function line(n,x1,y1,x2,y2,color=INK,width=2,arrow=false){
  const s=deck.slides.getItem(n-1).shapes.add({geometry:'connector',kind:'straight',position:{left:Math.min(x1,x2),top:Math.min(y1,y2),width:Math.abs(x2-x1),height:Math.abs(y2-y1),horizontalFlip:x2<x1,verticalFlip:y2<y1},line:{fill:color,width},...(arrow?{tail:{type:'triangle',width:'med',length:'med'}}:{})});s.bringToFront();return s;
 }
 function dim(n,x1,x2,y,color=INK){line(n,x1,y,x2,y,color);line(n,x1,y-6,x1,y+6,color);line(n,x2,y-6,x2,y+6,color);}
 function im(n,file,x,y,w,h,crop){
  const blob=fs.readFileSync(path.join(workdir,'assets',file));
  if(crop){const iw=blob.readUInt32BE(16),ih=blob.readUInt32BE(20),cw=iw*(1-crop.left-crop.right),ch=ih*(1-crop.top-crop.bottom),s=Math.min(w/cw,h/ch),nw=cw*s,nh=ch*s;x+=(w-nw)/2;y+=(h-nh)/2;w=nw;h=nh;}
  return deck.slides.getItem(n-1).images.add({blob,contentType:'image/png',fit:'contain',alt:file,position:{left:x,top:y,width:w,height:h},...(crop?{crop}:{})});
 }
 function title(n,t){tx(n,t,18,19,850,81,31,true,GREEN);}

 wipe(13);title(13,'BPA length is measured between the inner edges of the clamps');
 im(13,'s013_img0.png',160,170,640,255);
 tx(13,'Contracted\nat 620 kPa',25,208,155,70,24,true,GREEN);
 tx(13,'Uninflated',25,343,145,35,24,true,GREEN);
 // Original 668-pixel image: inner boundaries x245/537 (top), x125/540 (bottom).
 const topA=160+640*245/668,topB=160+640*537/668,botA=160+640*125/668,botB=160+640*540/668;
 dim(13,topA,topB,153);tx(13,'l₆₂₀',520,116,90,37,29,true,INK);
 line(13,topA,160,topA,195,GRAY,1);line(13,topB,160,topB,195,GRAY,1);
 dim(13,botA,botB,445);tx(13,'l₀',475,448,60,33,28,true,INK);
 line(13,botA,405,botA,438,GRAY,1);line(13,botB,405,botB,438,GRAY,1);
 tx(13,'lₘ = current\nBPA length',810,304,135,80,23,false,INK);
 deck.slides.getItem(12).speakerNotes.text='Resting length l0 is measured between the inner edges of the hose clamps, excluding both end fittings. l620 is measured between the same boundaries after free contraction at 620 kPa. lm is the current inner-clamp distance in a general operating condition. The photographed product has integrated end fittings; the dimension ends are placed where the rubber section emerges from them. Source: submitted dissertation Chapter 3, definitions and test procedure; Bolen et al., Actuators 15(5):230.';

 wipe(35);title(35,'A constant offset corrects the predicted BPA length');
 im(35,'s013_img0.png',115,120,650,135,{left:0,top:.5,right:0,bottom:0});
 // Measured active portion aligns with the pictured rubber; ideal reference is longer.
 const xL=115+650*125/668,xM=115+650*540/668,xG=xM+48;
 line(35,xL,244,xL,400,GRAY,1);line(35,xM,244,xM,400,GRAY,1);line(35,xG,276,xG,357,GRAY,1);
 dim(35,xL,xG,296,GRAY);tx(35,'Geometric BPA length',235,257,430,34,25,false,GRAY);
 dim(35,xL,xM,389,GREEN);tx(35,'Corrected lₘ',300,400,300,35,27,true,GREEN);
 // Positive chi0 points INWARD from geometric endpoint to corrected endpoint.
 line(35,xG,343,xM,343,GREEN,3,true);tx(35,'χ₀ > 0',xM-20,305,130,35,27,true,GREEN);
 tx(35,'lₘ = l',255,449,100,36,28,true,INK);tx(35,'LMT',318,465,47,24,17,true,INK);
 tx(35,'− 2 l',365,449,85,36,28,true,INK);tx(35,'fitting',419,465,57,24,17,true,INK);tx(35,'− χ₀',477,449,105,36,28,true,INK);
 deck.slides.getItem(34).speakerNotes.text='The constant term corrects the predicted current BPA length, not the measured resting length l0. For positive chi0, subtracting the offset moves the right endpoint inward and reduces lm. The diagram uses the positive flexor example and is not to scale. The adopted flexor chi0 is +3.94 mm; the extensor chi0 is −6.38 mm, which increases the corresponding corrected length. This is an effective fitted offset, not an additional end-cap dimension. Source: dissertation Chapter 4, equation Lm and AdoptedXi table.';

 wipe(37);title(37,'Offsets, compliance, and bending correct the predicted BPA length');
 tx(37,'χ₀: length offset',30,120,275,35,25,true,GREEN);
 im(37,'s037_img0.png',30,186,270,55);
 const a=80,b=248,c=282;
 dim(37,a,c,291,GRAY);tx(37,'Geometric length',40,257,260,30,21,false,GRAY);
 dim(37,a,b,366,GREEN);tx(37,'Corrected lₘ',70,378,230,35,24,true,GREEN);
 line(37,c,327,b,327,GREEN,2.5,true);tx(37,'χ₀ > 0',210,292,105,32,23,true,GREEN);
 line(37,a,244,a,372,GRAY,1);line(37,b,244,b,372,GRAY,1);line(37,c,283,c,334,GRAY,1);
 tx(37,'Positive flexor example',30,438,280,30,21,false,GRAY);

 tx(37,'χ₁, χ₂: system stiffness',335,120,300,40,25,true,GREEN);
 // Rotate the source photo so its femur is vertical, with proximal end above knee.
 const photo=im(37,'knee_cutout_approved.png',330,239,300,142);photo.rotation=270;
 tx(37,'χ₁, χ₂',355,171,135,42,27,true,GREEN);
 line(37,418,204,480,215,GREEN,3,true);

 tx(37,'χ₃: bending loss',665,120,275,35,25,true,GREEN);
 im(37,'s037_img2.png',660,153,250,310);
 // Length loss follows the wrapped R*alpha arc, rather than a radial arrow.
 const cx=660+250*.635,cy=153+310*.572,r=250*.252;
 let previous=null;
 for(let i=0;i<=14;i++){const phi=.12+i*(1.05-.12)/14;const p=[cx+r*Math.cos(phi),cy+r*Math.sin(phi)];if(previous)line(37,previous[0],previous[1],p[0],p[1],MAGENTA,4);previous=p;}
 tx(37,'Δl',882,351,65,37,27,true,MAGENTA);
 tx(37,'Δl = χ₃ R α (1 − ε*)²',652,455,300,31,23,true,INK);
 deck.slides.getItem(36).speakerNotes.text='Chi0 is a signed constant correction to geometric BPA length; the pictured positive flexor example shortens the corrected lm. Chi1 and chi2 are effective system stiffnesses for attachment deflection, lumping brackets, fixtures, artificial tendons, and the loading apparatus. The callout identifies the attachment system; its leader is not a bracket coordinate axis or a displacement vector. Exact principal-axis assignments depend on configuration and are given in Appendix C. Chi3 scales length loss due to bending along the wrapped arc R alpha. The adopted chi0 signs differ between flexor and extensor. Photograph: original test apparatus, with user-approved background removal. Source: dissertation Chapter 4 and Appendix C.';

 wipe(32);title(32,'The same attachment locations produce different knee torques with BPAs');
 tx(32,'Human muscle',75,126,365,40,28,true,GREEN);tx(32,'BPA',550,126,365,40,28,true,GREEN);
 im(32,'s032_img0.png',30,172,440,293,{left:.045,top:.333,right:.51,bottom:.333});
 im(32,'s032_img0.png',495,172,440,293,{left:.535,top:.333,right:.005,bottom:.333});
 deck.slides.getItem(31).speakerNotes.text='The knee torque surfaces compare biceps femoris long head in the human model with the corresponding BPA route. Corresponding attachment locations do not reproduce the human torque surface because the force-length properties differ. The color axes have different numerical ranges; compare their values and surface shapes, not the same color across panels. The original figure also includes hip torques, which are omitted here for focus. Source: Bolen and Hunt, 2019, Figure 8; submitted dissertation Chapter 2.';

 wipe(33);title(33,'The test stand measures knee torque at a fixed joint angle');
 im(33,'s033_img0.png',28,120,590,353,{left:.035,top:.591,right:.42,bottom:0});
 const rig=im(33,'s037_img1.png',555,209,350,165.7);rig.rotation=270;
 deck.slides.getItem(32).speakerNotes.text='The fixed femur is mounted to the test frame. The load cell connects the tibia to the swing arm; the winch string fixes the arm at each angle. The actuator is inflated and deflated with the winch locked, making the measurement as isometric as possible. The knee axis is coaxial with gravity, so gravity does not contribute knee torque. Keep swing arm, winch string, load-cell string, leg, BPA, and pneumatic lines visible. Source: submitted dissertation Chapter 4, Test Setup and Procedure.';

 const CORAL='#FA8775',PURPLE='#CC33B5',BLUE='#4040FF',GOLD='#FFE140';
 function legend(n,x,y,items){items.forEach(([label,color,kind],i)=>{const yy=y+i*43;
  if(kind==='point'){deck.slides.getItem(n-1).shapes.add({geometry:'ellipse',position:{left:x+10,top:yy+7,width:12,height:12},fill:color,line:{fill:'none',width:0}});}else if(label==='Corrected'){line(n,x,yy+13,x+40,yy+13,color,3);}else{for(const dx of [0,15,30])line(n,x+dx,yy+13,x+dx+9,yy+13,color,3);}
  tx(n,label,x+50,yy,205,34,23,true,INK);
 });}
 for(const n of [34,38,40]){
  wipe(n);title(n,n===34?'Improved BPA force prediction leaves a large knee-torque error':n===38?'The design factors improve predictions in the four pinned-flexor tests':'The same stiffness pair improves pinned-extensor torque predictions');
  im(n,`s${String(n).padStart(3,'0')}_img0.png`,15,110,755,357);
  const list=n===34?[['Measured',BLUE,'point'],['Hybrid',GOLD,'point'],['Hunt model','#666666','line'],['BPA model',CORAL,'line']]:n===38?[['Measured',BLUE,'point'],['BPA model',CORAL,'line'],['Corrected',PURPLE,'line']]:[['Measured',BLUE,'point'],['Hybrid',GOLD,'point'],['BPA model',CORAL,'line'],['Corrected',PURPLE,'line']];
  legend(n,690,125,list);
  const ph=im(n,`s${String(n).padStart(3,'0')}_img1.png`,812,318,122,150);if(n===40)ph.rotation=180;
  deck.slides.getItem(n-1).speakerNotes.text=(n===34?'The improved BPA force model is insufficient to explain measured joint torque. ':n===38?'Compare optimized predictions with the same measured flexor data shown before correction. ':'Extensor predictions use the same effective stiffness pair as the flexor, with extensor-specific chi0 and chi3. ')+'Hybrid means the actuator force model uses measured muscle length and measured moment arm, rather than their geometric predictions. The original four-panel flexor layout and six-panel extensor layout are retained. Effective stiffness includes the complete setup and is not a literal bracket material measurement. Source: dissertation Chapter 4 and Appendix C.';
 }

 wipe(36);title(36,'Attachment deflection changes the BPA length under load');
 const setup=im(36,'knee_cutout_approved.png',85,202,365,173);setup.rotation=270;
 // Rotated image has its femoral attachment above the knee attachment.
 tx(36,'Femur attachment',555,150,345,40,27,true,GREEN);
 line(36,550,179,255,177,GREEN,2,true);
 tx(36,'Tibia attachment',555,315,345,40,27,true,GREEN);
 line(36,550,343,258,335,GREEN,2,true);
 tx(36,'χ₁ and χ₂ describe effective\nstiffness of the complete setup',515,406,420,66,25,true,GREEN);
 deck.slides.getItem(35).speakerNotes.text='Actuator force deflects the attachments, fixtures, and tendons. Chi1 and chi2 describe the effective system stiffness needed to predict the resulting BPA length. The two labeled attachment locations orient the audience in the setup; the photograph itself is not a measurement of deflection. The winch string, swing arm, load-cell string, leg, BPA, and air lines have been preserved. Source: dissertation Chapter 4.';

 wipe(39);title(39,'The corrections reduce mean torque error on held-out configurations');
 tx(39,'Pinned flexor',75,130,380,44,30,true,GREEN);tx(39,'Pinned extensor',540,130,385,44,30,true,GREEN);
 tx(39,'Mean validation RMSE (N·m)',245,193,580,40,27,false,INK);
 tx(39,'7.588  →  1.490',75,267,390,65,39,true,GREEN);tx(39,'2.554  →  1.098',540,267,395,65,39,true,GREEN);
 tx(39,'80% reduction',90,376,360,55,33,true,GREEN);tx(39,'57% reduction',555,376,370,55,33,true,GREEN);
 tx(39,'Before → after correction',300,444,500,35,24,false,INK);
 deck.slides.getItem(38).speakerNotes.text='These are mean RMSE values on the adopted held-out sets. Flexor mean RMSE falls from 7.588 to 1.490 N m; extensor from 2.554 to 1.098 N m. Reported whole-percent reductions are 80% and 57%. They are not individual-test reductions. The adopted flexor split trains on two configurations and validates on two; extensor uses three and three, and shares the same effective stiffness pair. One flexor configuration and three extensor configurations were not used in identification. Source: dissertation Chapter 4 Tables FlxPinGof and ExtPinGof, lines538–541 and580–583; summary line651.';
}
