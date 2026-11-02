import fs from 'node:fs';
export async function apply(h) {
 const {d,rec,begin,tx,line,im,changed,replaceText,GREEN,INK,a}=h;
 // Native text boxes are sized against PowerPoint's rendered text bounds.
 for(const [n,rx,extra] of [[11,/^Force across|^Joint torque|^Tests of coupled/,8],[13,/Free contraction/,15],[18,/Measure F/,8]]) {
  for(const e of rec.get(n).elements.filter(e=>e.kind==='textbox'&&rx.test(e.text??'')&&e.position.top>100)) {d.resolve(e.id).position={...e.position,height:e.position.height+extra};changed.push(n);}
 }
 replaceText(61,/^Commands$/,'Commands',{left:27,top:261,width:102,height:44});
 replaceText(61,/^Sensory\nmeasurements$/,'Feedback',{left:264,top:261,width:130,height:44});
 const timing=rec.get(62).elements.find(e=>/^Timing with/.test(e.text??''));d.resolve(timing.id).position={...timing.position,height:66};changed.push(62);
 replaceText(42,/^The same stiffness pair/,'The same stiffness pair improves pinned-extensor torque predictions',{left:18,top:19,width:860,height:81});
 replaceText(93,/^Finite motion with/,'Motion with external assistance');
 const thanks=rec.get(74).elements.find(e=>e.text==='Thank you');if(thanks)d.delete([thanks.id]);changed.push(74);
 const rev=h.added.reverse_pulley;replaceText(rev,/^Design factors/,'Design factors can guide a more biomimetic knee mechanism');
 // The complete 12-edge selected RF/HAM graph is laid out for a wide screen.
 const notes=d.slides.getItem(64).speakerNotes.text;
 begin(65,'Rectus femoris and hamstring pathways combine hip and knee inputs',notes+' The presentation schematic uses native editable shapes and preserves the four cells and all twelve incident connections of the selected source ledger. Excitatory arrowheads and inhibitory terminal circles are defined on the slide.', 'Selected ground9 RC · Actual RF/HAM subnetwork · Same twelve connections as the source ledger');
 const edge=(points,inhib=false)=>{for(let i=1;i<points.length;i++)line(65,...points[i-1],...points[i],i===points.length-1&&!inhib,'#222222',2);if(inhib){const [x,y]=points.at(-1);d.slides.getItem(64).shapes.add({geometry:'ellipse',position:{left:x-4,top:y-4,width:8,height:8},fill:'#222222',line:{fill:'#222222',width:1}});}};
 const node=(text,x,y,w,ht,ellipse=false)=>{d.slides.getItem(64).shapes.add({geometry:ellipse?'ellipse':'rect',position:{left:x,top:y,width:w,height:ht},fill:ellipse?'#F4F7F5':'#FFFFFF',line:{fill:'#37443B',width:1.8}});tx(65,text,x+6,y+9,w-12,ht-12,ellipse?20:20,true,INK,'center');};
 // Inputs, reciprocal inhibition, and motor outputs; no connections inferred.
 edge([[215,177],[293,195]]);edge([[215,257],[293,225]]);
 edge([[410,210],[525,210]]);edge([[587.5,240],[587.5,315]],true);
 edge([[525,345],[410,345]]);edge([[347.5,315],[347.5,240]],true);
 edge([[347.5,180],[347.5,153],[737,153]]);
 edge([[325,420],[560,370]]);edge([[535,420],[615,370]]);
 edge([[650,345],[706,345],[706,278],[737,278]]);
 edge([[650,345],[737,345]]);
 edge([[650,345],[706,345],[706,412],[737,412]]);
 node('Hip flexor PF',35,151,180,52);node('Knee extensor PF',35,231,180,52);
 node('RF\nintegrator',285,180,125,60,true);node('RF\nrelay',525,180,125,60,true);
 node('HAM\nintegrator',525,315,125,60,true);node('HAM\nrelay',285,315,125,60,true);
 node('Hip extensor PF',235,420,180,50);node('Knee flexor PF',445,420,180,50);
 node('Rectus femoris MN',737,127,190,52);
 node('Semimembranosus',737,252,190,52);node('Semitendinosus',737,319,190,52);node('Biceps femoris\nlong head',737,386,190,62);
 tx(65,'PF: pattern formation     RF: rectus femoris     HAM: hamstrings',35,105,880,32,20);
 tx(65,'Motor neurons',737,213,190,31,20,true,GREEN,'center');
 tx(65,'Arrows: excitation\nCircles: inhibition',30,343,235,65,20);
 // Replace derivatives only after source-preserving regeneration has completed.
 const ham=a('neural/slide_panels/NEURAL_RIGHT_HAM_HIERARCHY.svg');
 if(fs.existsSync(ham)) {const imgs=rec.get(66).elements.filter(e=>e.kind==='image'&&e.position.top>140);d.delete(imgs.map(e=>e.id));im(66,ham,25,165,910,292);d.slides.getItem(65).speakerNotes.text=d.slides.getItem(65).speakerNotes.text.replace('knee-MN traces','biceps-femoris-long-head MN traces').replace('all hip/knee/ankle MN selections','available hip/knee/ankle MN selections');changed.push(66);}
 const stim=a('neural/slide_panels/DENG_RG_STIMULUS_COMPLETE.svg');
 if(fs.existsSync(stim)) {d.delete(rec.get(63).elements.filter(e=>(e.kind==='image'&&e.position.left>400)||/^Time from stimulus/.test(e.text??'')).map(e=>e.id));im(63,stim,405,175,525,268);changed.push(63);}
 for(const n of [67,h.added.neural_support]) {const e=rec.get(n).elements.find(e=>e.kind==='image'&&e.position.left>500);if(e){d.delete([e.id]);im(n,a('neural/slide_panels/FOOT_LOADS_AND_ASSISTANCE_WIDE.svg'),e.position.left,e.position.top,e.position.width,e.position.height);changed.push(n);}}
 const ct=rec.get(68).elements.find(e=>e.kind==='image');d.delete([ct.id]);im(68,a('neural/slide_panels/CONTACT_PF_AND_FOOT_LOADING_WIDE.svg'),45,102,870,333);changed.push(68);
 const shortfall=a('mechanics/single_bpa_shortfall_clean.png');
 if(fs.existsSync(shortfall)) {d.delete(rec.get(45).elements.filter(e=>e.kind==='image'&&e.position.top>100).map(e=>e.id));im(45,shortfall,90,105,780,369);changed.push(45);}
 const hy=rec.get(33).elements.find(e=>e.kind==='image'&&e.position.top>100);d.delete([hy.id]);im(33,a('mechanics/pinned48_preoptimization_torque.png'),75,105,810,342);tx(33,'Hybrid: measured BPA length, pressure, and moment arm',75,449,810,28,21,false,INK,'center');d.slides.getItem(32).speakerNotes.text+=' Hybrid calculation: BPA force is calculated from the measured inflated length and pressure, then multiplied by the measured ICR-to-muscle distance (moment arm). Source: Plot_KneeFlxPin10mm_48cm.m, lines 83–87. It is distinct from a pressure-only substitution.';changed.push(33);
}
