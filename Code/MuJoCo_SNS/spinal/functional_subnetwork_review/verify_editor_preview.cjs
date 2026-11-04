const assert=require('node:assert/strict'),fs=require('node:fs'),crypto=require('node:crypto');
const api=require('./schematic_editor.js'),root='D:/Users/Ben Bolen/NeuroGraphics_20261102/editor_preview/';
const hash=x=>crypto.createHash('sha256').update(JSON.stringify(x)).digest('hex');
const results=[];
for(let i=1;i<=4;i++){
 const doc=JSON.parse(fs.readFileSync(root+`example_${i}.json`,'utf8'));api.validate(doc);
 const sourceHash=hash(doc.sourceGraph),copy=api.roundTrip(doc);
 assert.deepEqual(copy.sourceGraph,doc.sourceGraph);assert.deepEqual(copy.views,doc.views);
 copy.views[0].glyphs[0].x+=43;copy.views[0].glyphs[0].label='display-only edit';
 copy.views[0].wires[0].bends.push({x:37,y:91});
 assert.equal(hash(api.roundTrip(copy).sourceGraph),sourceHash);
 let count=0;
 for(const w of doc.views[0].wires){
  const points=api.verts(doc.views[0],w);assert(points.every(p=>Number.isFinite(p.x)&&Number.isFinite(p.y)));
  const t=api.terminal(points,w.sign);assert.equal(t.shape,w.sign==='inh'?'circle':w.sign==='physical'?'polyline':'polygon');
  if(t.shape==='polygon'){const b=t.points.slice(1);assert(Math.abs((b[0].x+b[1].x)/2-points.at(-1).x)<1e-7);assert(Math.abs((b[0].y+b[1].y)/2-points.at(-1).y)<1e-7);}
  const path=api.routedPath(points,[]);assert(path.startsWith('M '));count+=w.records.length;
 }
 results.push({example:i,source_identity:doc.identity.title,view_glyphs:doc.views[0].glyphs.length,represented_records:count,source_graph_exact_round_trip:true,source_unchanged_after_layout_edit:true});
}
assert(api.routedPath([{x:0,y:50},{x:100,y:50}],[[{x:50,y:0},{x:50,y:100}]]).includes('A 7 7'));
assert(api.routedPath([{x:0,y:0},{x:100,y:0},{x:100,y:100}],[]).includes('Q '));
assert.equal(fs.readFileSync('Code/MuJoCo_SNS/spinal/functional_subnetwork_review/schematic_editor.js','utf8').includes('localStorage'),false);
const report={examples:results,terminal_base_on_destination:true,semicircular_crossings:true,rounded_bends:true,no_browser_storage_api:true,circuits_executed:0};
fs.writeFileSync(root+'VERIFICATION.json',JSON.stringify(report,null,2));console.log(JSON.stringify(report));
