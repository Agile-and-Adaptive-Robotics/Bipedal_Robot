const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict'),path=require('node:path'),crypto=require('node:crypto');
const root=path.resolve(__dirname,'..'),html=fs.readFileSync(path.resolve(root,'../spinal/connectome_block_editor.html'),'utf8'),original=html.slice(html.lastIndexOf('<script>')+8,html.lastIndexOf('</script>'));
function Fake(){this.children=[];this.dataset={};this.style={};this.classList={add(){},toggle(){}};}Fake.prototype.setAttribute=function(){};Fake.prototype.addEventListener=function(){};Fake.prototype.append=function(...x){this.children.push(...x);};Fake.prototype.appendChild=function(x){this.children.push(x);};Fake.prototype.remove=function(){};Fake.prototype.querySelector=function(){return null;};Fake.prototype.querySelectorAll=function(){return[];};Fake.prototype.replaceChildren=function(){};
const sandbox={console,crypto:{randomUUID:crypto.randomUUID},setTimeout,clearTimeout,window:{},document:{createElement:()=>new Fake(),createElementNS:()=>new Fake(),getElementById:()=>new Fake(),querySelectorAll:()=>[]},NS:'svg',TYPES:{MN:{shape:'circle'}},nodes:[],edges:[],models:[{name:'test',data:{}}],activeTab:0,ctxStack:[],selSet:new Set(),sel:null,hiddenG:new Set(),grpNames:{},uid:1,lastMouseWorld:{x:2000,y:2000},undoStack:[],redoStack:[],status(){},prompt:()=> 'Packed',confirm:()=>true,addNode(){},renderNode(){},enterSub(){},exitSub(){},redrawEdges(){},renderTree(){},clearSel(){sandbox.sel=null;sandbox.selSet.clear();},applySelClasses(){},updatePanel(){},persist(){},startDrag(){},renderTabs(){},applyView(){},view:{},svg:{getBoundingClientRect:()=>({width:1200,height:800})},isHiddenNode:()=>false};
sandbox.selectEdge=id=>{sandbox.sel={kind:'edge',id};};sandbox.selectNode=id=>{sandbox.sel={kind:'node',id};};
vm.createContext(sandbox);vm.runInContext(original.slice(original.indexOf('function modelSpec()'),original.indexOf('function snapshot()')),sandbox);vm.runInContext(fs.readFileSync(path.join(root,'desktop_editor.js'),'utf8'),sandbox);
vm.runInContext(original.slice(original.indexOf('function uniqueLabel('),original.indexOf('function copySel(')),sandbox);
vm.runInContext(original.slice(original.indexOf('function nodeById('),original.indexOf('// ---------------- node rendering')),sandbox);
sandbox.snapshot=()=>sandbox.undoStack.push(JSON.stringify(vm.runInContext('serialize()',sandbox)));
vm.runInContext('function serialize(){return specOf(nodes,edges,grpNames,hiddenG)}',sandbox);
const fixture=JSON.parse(fs.readFileSync(path.join(root,'examples/paired_mn.sns.json'),'utf8'));
fixture.inputs=[{name:'explicit afferent port',node_id:fixture.nodes[0].id}];
function load(g){sandbox.g=g;vm.runInContext('constLoaded=materialize(g);nodes=constLoaded.nodes;edges=constLoaded.edges;',sandbox);}
function documentNow(){return JSON.parse(vm.runInContext('JSON.stringify(serialize())',sandbox));}
function semantics(g){const out=[];function walk(g){const by=new Map(g.nodes.map(n=>[n.id,n]));function resolve(e,s){let n=by.get(e[s+'Id']);for(const id of e[s+'Path']||[])n=n.sub.nodes.find(n=>n.id===id);assert(n);return n.id;}for(const e of g.edges)out.push([e.id,resolve(e,'from'),resolve(e,'to'),e.synapse,e.source_record,e.map]);for(const n of g.nodes)if(n.sub)walk(n.sub);}walk(g);return out.sort((a,b)=>a[0].localeCompare(b[0]));}
load(fixture);const before=documentNow(),meaning=semantics(before);assert.equal(before.nodes.length,32);assert.equal(before.edges.length,50);
const neuronFields=before.nodes.map(n=>[n.id,n.model,n.source_record]);
sandbox.selSet=new Set([sandbox.nodes[0].id,sandbox.nodes[3].id]);vm.runInContext('packSub()',sandbox);let packed=documentNow();assert.deepEqual(semantics(packed),meaning);assert.deepEqual(packed.source,before.source);
assert.equal(packed.inputs[0].node_path[0],fixture.nodes[0].id);assert.equal(packed.inputs[0].node_id,sandbox.nodes.find(n=>n.sub).id);
const host=sandbox.nodes.find(n=>n.sub);sandbox.selSet=new Set([host.id,sandbox.nodes[0].id]);vm.runInContext('packSub()',sandbox);const twice=documentNow();assert.deepEqual(semantics(twice),meaning);
load(twice);sandbox.sel={kind:'node',id:sandbox.nodes.find(n=>n.sub).id};vm.runInContext('unpackSub()',sandbox);sandbox.sel={kind:'node',id:sandbox.nodes.find(n=>n.sub).id};vm.runInContext('unpackSub()',sandbox);assert.deepEqual(semantics(documentNow()),meaning);
const after=documentNow();assert.deepEqual(after.nodes.map(n=>[n.id,n.model,n.source_record]).sort(),neuronFields.sort());
assert.equal(after.inputs[0].node_id,fixture.nodes[0].id);assert.deepEqual(after.inputs[0].node_path,[]);
load(packed);sandbox.selSet=new Set(sandbox.nodes.map(n=>n.id));vm.runInContext('copySel();pasteClip()',sandbox);const copied=documentNow();function ids(g){return g.nodes.flatMap(n=>[n.id,...(n.sub?ids(n.sub):[])]);}assert.equal(ids(copied).length,new Set(ids(copied)).size);assert.equal(semantics(copied).length,meaning.length*2);
const sourceParameter=fixture.nodes.find(n=>n.model.kind==='graded').model.capacitance_nF;assert.equal(documentNow().nodes.find(n=>n.model?.kind==='graded').model.capacitance_nF,sourceParameter);
assert.equal(copied.inputs.length,2);for(const p of copied.inputs){let n=copied.nodes.find(n=>n.id===p.node_id);for(const id of p.node_path||[])n=n.sub.nodes.find(n=>n.id===id);assert(n);}
// Saving from a still-focused parameter field must consume its latest input.
vm.runInContext("testEditedEdge=edges.find(e=>e.synapse);testInput=desktopField(document.createElement('div'),'g',testEditedEdge.synapse.gmax_uS,v=>testEditedEdge.synapse.gmax_uS=v);testInput.value='0.375';testInput.oninput();",sandbox);
assert(documentNow().edges.some(e=>e.synapse?.gmax_uS===0.375));
console.log('PASS: all 50 connections survive nested pack/save/reload/unpack; complete node/source fields preserved; subtree copy remaps boundary IDs without aliasing.');
