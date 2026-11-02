import fs from 'node:fs/promises';
export async function apply(h){
 const {d,rec,begin,im,tx,a,old,changed,GREEN,INK}=h;
 const snap=await d.inspect({select:'$.slides[*]',include:['id','index','elements.id','elements.kind','elements.position','elements.text']});
 rec.clear();for(const r of snap.records)rec.set(r.index+1,r);
 await fs.writeFile('build/pass3_elements.json',JSON.stringify(snap.records));
 const find=(n,rx)=>rec.get(n).elements.find(e=>e.kind==='textbox'&&rx.test(e.text??''));
 function replace(n,rx,text,pos){const e=find(n,rx);if(!e)throw Error('Missing '+n+' '+rx);const ob=d.resolve(e.id);ob.text=text;if(pos)ob.position=pos;changed.push(n);return ob;}
 function remove(n,p){const ids=rec.get(n).elements.filter(p).map(e=>e.id);if(ids.length)d.delete(ids);changed.push(n);}
 function white(n,x,y,w,h){d.slides.getItem(n-1).shapes.add({geometry:'rect',position:{left:x,top:y,width:w,height:h},fill:'#FFFFFF',line:{fill:'none',width:0}});}
 replace(3,/^BPAs provide/,'Braided pneumatic actuators provide lightweight, compliant actuation');
 replace(9,/Controlled interventions/,'Robot experiments generate\nbiological predictions');
 replace(9,/Improve autonomous/,'Goal: Better autonomous movement on unstructured terrain');
 replace(24,/Normalization separates/,'Normalization combines measurements across initial lengths');
 d.slides.getItem(23).speakerNotes.text='Force is normalized by the maximum reference force, and contraction by the measured free contraction. The two diameter classes retain their own fitted coefficients. The published study used 537 measured data points, with 80% for fitting and 20% for testing. The selected-length comparison is not a separately held-out-actuator study. Source: Bolen et al. (2026), Actuators15(5):230.';
 replace(125,/^Force displaces/,'Force displaces the attachments',{left:475,top:133,width:440,height:70});
 replace(125,/^Displacement changes/,'Displacement changes\nBPA length and force',{left:475,top:235,width:440,height:75});
 tx(125,'δ: displacement along the BPA path\nfrom the series compliance',475,329,440,57,22,false,INK);
 tx(125,'Normalized BPA force is dimensionless',475,391,440,28,22,false,INK);
 replace(129,/Fit the parameters/,'Fit torque residuals');
 d.slides.getItem(128).speakerNotes.text+='\nThe identification objectives are baseline-normalized per-test RMSE, fraction of variance unexplained, and maximum absolute residual, averaged across the fitting tests. gamultiobj returns Pareto candidates; validation behavior also informed their acceptance. This is a selected, cross-configuration fit, not blind external validation. The stiffness pair represents effective system compliance and is not uniquely identifiable as two bracket material properties.';changed.push(129);
 replace(43,/Predict torque at every joint angle/,'Predict torque across the joint range\nin the required direction',{left:40,top:325,width:510,height:70}).text.style={fontFamily:'Arial',fontSize:25,color:GREEN,bold:true,alignment:'center'};
 // These small masks cover only the clipped source-axis fragment outside motion.
 white(43,660,112,37,28);white(127,228,111,41,30);changed.push(127);
 tx(46,'BPA path',380,206,210,36,24,true,'#0000FF','center');
 tx(46,'Black points:\ncam contacts',380,268,210,64,23,true,INK,'center');
 tx(46,'Circles: clearance\nand bend constraints',365,353,240,67,21,false,INK,'center');changed.push(46);
 // Keep the apparatus image below the original legend and outside the plot axes.
 const pin=rec.get(40).elements.find(e=>e.kind==='image'&&e.position?.left>810&&e.position?.top>150);
 if(!pin)throw Error('Pinned-extensor inset missing');d.resolve(pin.id).position={left:826,top:320,width:105,height:150};
 remove(40,e=>e.kind==='textbox'&&/χ₁, χ₂ fixed/.test(e.text??''));
 // The first overview uses correction factors; it is not the later route redesign.
 white(42,715,317,207,34);tx(42,'Rigid-body model',716,319,206,17,13.5,true,INK);tx(42,'Corrected model',716,335,206,17,13.5,true,INK);changed.push(42);
 // Restore the upper source heading rather than leaving a clipped fragment.
 const rg=rec.get(53).elements.find(e=>e.kind==='image'&&e.position?.left<310);
 if(!rg)throw Error('RG crop missing');d.delete(rg.id);
 im(53,old('s053_img0.png'),35,130,270,285,{left:.065,top:0,right:.80,bottom:.55});changed.push(53);
 replace(78,/Each curve is divided/,'Each curve is divided by its own peak',{left:150,top:460,width:700,height:19}).text.style={fontFamily:'Arial',fontSize:16,color:INK,alignment:'center'};
 remove(101,e=>e.kind==='textbox'&&e.position?.top>435&&e.position?.top<478);
 tx(101,'Phase shift: 0.000 · Blue: stimulus · Dashed: baseline',70,445,830,33,23,true,GREEN,'center');
 begin(70,'This work advances actuator prediction, joint design, and neural testing','The force model accounts for initial length, pressure, and contraction. Experimentally identified design factors improve joint-torque predictions and guide knee designs that meet the selected human targets at the measured angles. The neural tools support explicit, inspectable model interventions; the recovered-W2L feedback comparison illustrates one such test. Mechanical calibration, circuit implementation checks, supported stepping, and feedback interventions are distinct evidence. The next work is dynamic physical testing with measured support, pressure, force, joint motion, and contact.','Dissertation Conclusion');
 im(70,old('s025_img0.png'),30,145,280,235);
 im(70,old('s045_img0.png'),340,155,280,220);
 tx(70,'Recovered W2L',650,102,280,24,19,false,INK,'center');
 tx(70,'Feedback intact',650,130,280,26,21,true,GREEN,'center');
 im(70,old('s056_img0.png'),650,161,280,94,{left:0,top:.755,right:.46,bottom:0});
 tx(70,'Feedback removed',650,265,280,28,21,true,GREEN,'center');
 im(70,old('s056_img0.png'),650,295,280,100,{left:.535,top:.755,right:0,bottom:0});
 tx(70,'Predict\nBPA force',30,411,280,64,26,true,GREEN,'center');
 tx(70,'Meet selected\nknee-torque targets',340,411,280,64,26,true,GREEN,'center');
 tx(70,'Test feedback\nmechanisms',650,411,280,64,26,true,GREEN,'center');
 tx(70,'Thank you',390,480,180,30,25,true,GREEN,'center');
}
