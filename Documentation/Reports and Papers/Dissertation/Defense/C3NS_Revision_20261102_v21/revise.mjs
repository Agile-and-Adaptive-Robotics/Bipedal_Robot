import fs from 'node:fs/promises';
import fss from 'node:fs';
import path from 'node:path';
import {FileBlob,PresentationFile} from 'file:///C:/Users/Ben%20Bolen/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/@oai/artifact-tool/dist/artifact_tool.mjs';
const W=path.resolve('.'),R=path.resolve('../../../../..'),V20=path.resolve('../C3NS_Revision_20261102_v20');
const d=await PresentationFile.importPptx(await FileBlob.load(path.join(V20,'output/Bolen_PhD_Defense_RoseGreen_20261102_v20.pptx')));
const rec=new Map(),changed=[],GREEN='#174D24',INK='#202A25';let count=127;
const added={},after=new Map(),excluded=new Set([53]),backupExtra=[];
async function refresh(){const s=await d.inspect({select:'$.slides[*]',include:['id','index','elements.id','elements.kind','elements.position','elements.text']});rec.clear();for(const r of s.records)rec.set(r.index+1,r);return s;}
await refresh();await fs.writeFile('build/baseline_elements.json',JSON.stringify([...rec.values()]));
const a=n=>path.join(W,'assets',n),old=n=>path.join(V20,'assets',n);
function tx(n,text,x,y,w,h,size=25,bold=false,color=INK,alignment='left'){
 const s=d.slides.getItem(n-1).shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});s.text=text;s.text.insets=0;s.text.style={fontFamily:'Arial',fontSize:size,bold,color,alignment};return s;
}
function im(n,file,x,y,w,h,crop){
 const blob=fss.readFileSync(file),ext=path.extname(file).toLowerCase();
 if(ext==='.svg') { const v=blob.toString('utf8').match(/viewBox=["']\s*([-\d.eE]+)[ ,]+([-\d.eE]+)[ ,]+([\d.eE]+)[ ,]+([\d.eE]+)\s*["']/);if(!v)throw Error('SVG missing viewBox: '+file);const iw=Number(v[3]),ih=Number(v[4]),sc=Math.min(w/iw,h/ih),nw=iw*sc,nh=ih*sc;x+=(w-nw)/2;y+=(h-nh)/2;w=nw;h=nh; }
 if(ext==='.png'||ext==='.gif'){const iw=(ext==='.png'?blob.readUInt32BE(16):blob.readUInt16LE(6))*(1-(crop?.left??0)-(crop?.right??0)),ih=(ext==='.png'?blob.readUInt32BE(20):blob.readUInt16LE(8))*(1-(crop?.top??0)-(crop?.bottom??0)),sc=Math.min(w/iw,h/ih);const nw=iw*sc,nh=ih*sc;x+=(w-nw)/2;y+=(h-nh)/2;w=nw;h=nh;}
 const ct={'.jpg':'image/jpeg','.jpeg':'image/jpeg','.gif':'image/gif','.svg':'image/svg+xml','.png':'image/png'}[ext];if(!ct)throw Error('Unsupported image '+file);
 return d.slides.getItem(n-1).images.add({blob,contentType:ct,fit:'contain',alt:path.basename(file),position:{left:x,top:y,width:w,height:h},...(crop?{crop}:{})});
}
function line(n,x1,y1,x2,y2,arrow=true,color=GREEN,width=2.5){return d.slides.getItem(n-1).shapes.add({geometry:'connector',kind:'straight',position:{left:Math.min(x1,x2),top:Math.min(y1,y2),width:Math.abs(x2-x1),height:Math.abs(y2-y1),horizontalFlip:x2<x1,verticalFlip:y2<y1},line:{fill:color,width},...(arrow?{tail:{type:'triangle',width:'med',length:'med'}}:{})});}
function begin(n,title,notes,cite){const ids=(rec.get(n)?.elements??[]).filter(e=>e.position&&e.position.top<478).map(e=>e.id);if(ids.length)d.delete(ids);tx(n,title,18,19,850,82,title.length>85?30:32,true,GREEN);d.slides.getItem(n-1).speakerNotes.text=notes;const source=rec.get(n)?.elements.find(e=>e.position?.top>=478&&e.position?.top<510&&e.position?.width>650&&e.kind==='textbox');if(cite&&source)d.resolve(source.id).text=cite;changed.push(n);}
function newSlide(key,title,notes,cite){const n=++count;d.slides.add({width:960,height:540});rec.set(n,{elements:[]});added[key]=n;begin(n,title,notes,cite);tx(n,cite??'',34,481,850,28,10,false,'#52615A');tx(n,String(n),918,518,30,20,13,false,'#FFFFFF','center');return n;}
function replaceText(n,rx,value,pos){const e=rec.get(n).elements.find(e=>e.kind==='textbox'&&rx.test(e.text??''));if(!e)throw Error('Missing text on '+n+':'+rx);const s=d.resolve(e.id);s.text=value;if(pos)s.position=pos;changed.push(n);return s;}
function table(n,values,x,y,w,h,widths,size=22){const t=d.slides.getItem(n-1).tables.add({rows:values.length,columns:values[0].length,left:x,top:y,width:w,height:h,columnWidths:widths,values});t.cells.block({row:0,column:0,rowCount:values.length,columnCount:values[0].length}).assign({textStyle:{fontFamily:'Arial',fontSize:size,color:INK,alignment:'center'},margins:{left:5,right:5,top:4,bottom:4}});t.borders.assign({fill:'#D4DFD6',width:.7});for(let r=0;r<values.length;r++)for(let c=0;c<values[0].length;c++){t.getCell(r,c).fill=r===0?GREEN:(r%2===0?'#F0F6F1':'#FFFFFF');t.getCell(r,c).text.style={fontFamily:'Arial',fontSize:size,color:r===0?'#FFFFFF':INK,alignment:c===0?'left':'center',bold:r===0};}return t;}
function addAfter(n,...ids){after.set(n,[...(after.get(n)??[]),...ids]);}
const h={d,rec,refresh,begin,newSlide,im,tx,line,table,replaceText,changed,GREEN,INK,W,R,V20,a,old,addAfter,excluded,backupExtra,added};
await (await import('./force_content.mjs')).apply(h);
await refresh();await (await import('./mechanics_content.mjs')).apply(h);
await refresh();const neural=await (await import('./neural_content.mjs')).apply(h);
await refresh();await (await import('./backup_content.mjs')).apply(h);
await refresh();await (await import('./visual_fixes.mjs')).apply(h);
const main=[];
for(let n=1;n<=74;n++){
 // Keep the requested force-comparison tables next to the final force result;
 // the existing torque divider follows them rather than interrupting that evidence.
 if(n===27&&after.has(n))main.push(...after.get(n));
 if(n===57&&neural?.main){main.push(...neural.main);}
 if(neural?.main&&n>=57&&n<=68)continue;
 if(!excluded.has(n))main.push(n);
 if(n!==27&&after.has(n))main.push(...after.get(n));
}
const questions=[...Array.from({length:46},(_,i)=>75+i),...(neural?.backup??[]),...backupExtra,...Array.from({length:7},(_,i)=>121+i)];
const order=[...main,...questions];if(new Set(order).size!==order.length)throw Error('Duplicate slide in order');
await fs.writeFile('build/slide_order.json',JSON.stringify(order));
await fs.writeFile('build/changed_slides.json',JSON.stringify([...new Set(changed)]));
await fs.writeFile('build/structure.json',JSON.stringify({main,questions,added,source_count:count,total:order.length},null,2));
const blob=await PresentationFile.exportPptx(d);await blob.save('build/artifact_draft.pptx');
console.log(JSON.stringify({main:main.length,questions:questions.length,total:order.length,changed:[...new Set(changed)]}));
