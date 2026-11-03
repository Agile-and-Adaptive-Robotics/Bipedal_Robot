import fs from 'node:fs/promises';
import path from 'node:path';
import {FileBlob,PresentationFile} from 'file:///C:/Users/Ben%20Bolen/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/@oai/artifact-tool/dist/artifact_tool.mjs';
const W=process.cwd(),B=path.join(W,'build/plot_rotation_180');
const input=path.join(W,'output/Bolen_PhD_Defense_RoseGreen_20261102_v21.pptx');
const d=await PresentationFile.importPptx(await FileBlob.load(input));
const scan=await d.inspect({select:'$.slides[*]',include:['id','index','elements.id','elements.kind','elements.position']});
const blob=await fs.readFile(path.join(W,'assets/mechanics/historical_bflh_common_scales.png'));
const changes=[];
for(const n of [8,33]){
 const row=scan.records.find(r=>r.index===n-1);
 const matches=row.elements.filter(e=>e.kind==='image' && e.position.width>700 && e.position.height>250);
 if(matches.length!==1)throw Error(`Expected one comparison image on slide ${n}`);
 const image=d.resolve(matches[0].id);
 const saved={frame:image.frame,crop:image.crop,geometry:image.geometry,borderRadius:image.borderRadius,rotation:image.rotation,flipHorizontal:image.flipHorizontal,flipVertical:image.flipVertical,lockAspectRatio:image.lockAspectRatio};
 image.replace({blob,contentType:'image/png',alt:image.alt??'Human and robot torque surfaces',fit:image.fit??'contain'});
 Object.assign(image,saved);
 changes.push({slide:n,image_id:matches[0].id,position:matches[0].position});
}
await (await PresentationFile.exportPptx(d)).save(path.join(B,'artifact_images.pptx'));
await fs.writeFile(path.join(B,'artifact_edit.json'),JSON.stringify(changes,null,2));
console.log(JSON.stringify(changes));
