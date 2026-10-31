const fs=require('node:fs');
const path=require('node:path');
const crypto=require('node:crypto');
const sharp=require('C:/Users/Ben Bolen/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/sharp');
const current=__dirname;
const previous=path.join(path.dirname(current),'C3NS_Revision_20261030_v14');
const hash=b=>crypto.createHash('sha256').update(b).digest('hex');
(async()=>{
  const slides=[];
  for(let n=1;n<=83;n++){
    const file=`slide_${String(n).padStart(2,'0')}.png`;
    const [a,b]=await Promise.all([sharp(path.join(previous,'build/render',file)).ensureAlpha().raw().toBuffer(),sharp(path.join(current,'build/render',file)).ensureAlpha().raw().toBuffer()]);
    const same=hash(a)===hash(b);
    slides.push({slide:n,pixels_identical_to_individually_reviewed_v14:same,visual_review:same?'Exact pixel match to previously inspected slide':'Inspect revised slide at full size'});
  }
  const changed=slides.filter(s=>!s.pixels_identical_to_individually_reviewed_v14).map(s=>s.slide);
  const report={rendered:83,unchanged_pixel_matches:83-changed.length,changed_slides:changed,slides};
  fs.writeFileSync(path.join(current,'build/visual_review.json'),JSON.stringify(report,null,2));
  console.log(JSON.stringify({rendered:83,unchanged_pixel_matches:83-changed.length,changed_slides:changed}));
  if(changed.some(n=>![1,19].includes(n)))process.exitCode=1;
})().catch(e=>{console.error(e);process.exitCode=1});
