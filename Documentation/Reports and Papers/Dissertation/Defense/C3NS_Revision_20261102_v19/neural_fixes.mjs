import fs from 'node:fs';
import path from 'node:path';

// Review fixes only: preserve the submitted measurements and plot colors.
export async function fixNeural(deck, workdir) {
  const snapshot = await deck.inspect({select:'$.slides[*]',include:['id','index','elements.id','elements.kind','elements.text','elements.position']});
  const records = new Map(snapshot.records.map(r=>[r.index+1,r]));
  const green='#174D24', ink='#202A25';
  const slide=n=>deck.slides.getItem(n-1);
  function clearBody(n) {
    const ids=records.get(n).elements.filter(e=>e.position && e.position.top>=88 && e.position.top<478).map(e=>e.id);
    if(ids.length)deck.delete(ids);
  }
  function text(n,value,x,y,w,h,size=24,bold=false,color=ink) {
    const s=slide(n).shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
    s.text=value;s.text.insets=0;s.text.style={fontFamily:'Arial',fontSize:size,bold,color};return s;
  }
  function title(n,value) {
    const e=records.get(n).elements.find(e=>e.kind==='textbox'&&e.position?.top<80&&e.position?.width>650);
    if(e){const s=deck.resolve(e.id);s.text=value;s.text.style={fontFamily:'Arial',fontSize:31,bold:true,color:green};}
  }
  function cropImage(n,file,bounds,crop) {
    const blob=fs.readFileSync(file), width=blob.readUInt32BE(16),height=blob.readUInt32BE(20);
    const cw=width*(1-crop.left-crop.right),ch=height*(1-crop.top-crop.bottom),scale=Math.min(bounds.width/cw,bounds.height/ch);
    const w=cw*scale,h=ch*scale;
    return slide(n).images.add({blob,contentType:'image/png',alt:'Dissertation model mean, standard deviation, and OpenSim IK comparison',fit:'contain',crop,position:{left:bounds.left+(bounds.width-w)/2,top:bounds.top+(bounds.height-h)/2,width:w,height:h}});
  }

  clearBody(108);
  const table=slide(108).tables.add({rows:5,columns:3,left:45,top:130,width:870,height:300,columnWidths:[325,255,290],values:[
    ['Protocol','RG shift (cycles)','Motor shift (cycles)'],
    ['RG, two periods','0.211','0.211'],
    ['PF, two periods','0.000','0.000'],
    ['MN, two periods','0.000','0.000'],
    ['PF, one second','0.000','0.208'],
  ]});
  table.styleOptions={headerRow:false,bandedRows:false};
  table.borders.assign({style:'solid',fill:'#738078',width:1});
  for(let r=0;r<5;r++)for(let c=0;c<3;c++){
    const cell=table.getCell(r,c);cell.fill='#FFFFFF';
    cell.text.style={fontFamily:'Arial',fontSize:r===0?23:26,bold:r===0,color:r===0?green:ink};
  }
  slide(108).speakerNotes.text='The two-period stimuli delete two expected motor bursts. RG stimulation also changes the rhythm phase by 0.211 cycles. PF and MN stimulation in that protocol preserve the original schedule. A separate one-second PF pulse shifts motor timing by 0.208 cycles while the RG clock is unchanged. Thus stimulus duration and the measured layer both matter. Source: submitted dissertation Chapter 6 and Appendix D, reduced Deng-derived circuit; stimulus_results_all.csv.';

  clearBody(111);
  title(111,'Joint motion and stance differ from the OpenSim walking reference');
  const file=path.join(workdir,'assets','s111_img0.png');
  text(111,'Blue: model mean ± SD     Dashed black: OpenSim IK',100,112,810,38,23,false);
  cropImage(111,file,{left:25,top:152,width:442,height:298},{left:0,top:80/2400,right:0,bottom:.479});
  cropImage(111,file,{left:490,top:152,width:442,height:298},{left:0,top:.522,right:0,bottom:0});
  text(111,'Gait cycle (%)',150,447,300,30,22,false);
  slide(111).speakerNotes.text='The same eight original panels are rearranged into two groups for legibility: hip and knee on the left, ankle and stance fraction on the right. The blue curve and band show the model mean and standard deviation; black dashed curves show the OpenSim inverse-kinematics reference. These differences illustrate the limits of the supported result. A finite trace comparison is not proof of biological equivalence or global walking stability. Source: submitted dissertation Chapter 6 and Appendix D; original IK comparison figure. No plot values or colors changed.';

  for(const e of records.get(118).elements){
    if(e.kind==='textbox'&&e.text?.includes('Sensory Anatomy Database'))deck.resolve(e.id).text=e.text.replace('Sensory Anatomy Database','Sensory Afferent Database');
  }
  slide(95).speakerNotes.text='The historical syn6 recipe has nonzero heel/toe drive and zero pooled central Ib gain. General direct central Ia/II projections are absent, while local sensory/reflex routes remain. Heel excites both ipsilateral RG interneurons, whose outputs inhibit opposite half-centers. This matters in the historical contact-driven recipe, not only the later silent standing export. The intended combined phase response needs a functional test. Absence of direct central Ia/II is not a universal defect: Di Russo uses phase oscillators with heel-event resets. Syn6 combines NMF channels with half-centers and is an engineered hybrid, not a literal Di Russo replication. Source: GAIT_DISSERTATION_RECIPE_REVIEW.md and DEFENSE_DESIGN_DECISIONS.md, November 2, 2026.';
}
