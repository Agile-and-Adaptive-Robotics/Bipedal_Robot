const{sandbox:s,load,documentNow,vm}=require('./test_editor.cjs');
const assert=require('node:assert/strict');
const cell=(id,type)=>({id,type,label:id,x:0,y:0,model:{kind:'graded'}});
const layer={id:'RG',type:'SUB',label:'RG',x:700,y:240,scope:{layer:'RG'},sub:{nodes:['Left','Right'].map((side,i)=>({id:side,type:'SUB',label:side+' · 1 cells',scope:{layer:'RG',side},x:240+i*380,y:180,sub:{nodes:[cell('rg'+side,'MN')],edges:[]}})),edges:[{id:'cross',fromId:'Left',fromPath:['rgLeft'],toId:'Right',toPath:['rgRight'],synapse:{gmax_uS:.75,reversal_mV:-5,threshold_mV:0,saturation_mV:5},pts:[{x:400,y:180}]}]}};
const doc={nodes:[layer,cell('Ia','SN-Ia')],edges:[{id:'feedback',fromId:'Ia',toId:'RG',toPath:['Left','rgLeft'],synapse:{gmax_uS:.125,reversal_mV:8,threshold_mV:0,saturation_mV:5},pts:[]}],inputs:[{name:'input',node_id:'RG',node_path:['Left','rgLeft']}]};
load(doc);s.sel=null;s.selSet=new Set(['RG']);vm.runInContext('unpackSub()',s);const result=documentNow();
assert.deepEqual(result.nodes.filter(n=>n.scope?.layer==='RG').map(n=>n.scope.side),['Left','Right']);
assert(result.nodes.find(n=>n.id==='Left').label.startsWith('Left RG'));
const edge=result.edges.find(e=>e.id==='feedback');assert.equal(edge.toId,'Left');assert.deepEqual(edge.toPath,['rgLeft']);assert.deepEqual(edge.synapse,doc.edges[0].synapse);
assert.deepEqual(result.edges.find(e=>e.id==='cross').pts,[{x:670,y:240}]);
assert.equal(result.inputs[0].node_id,'Left');assert.deepEqual(result.inputs[0].node_path,['rgLeft']);
// Expanding the left RG one further step must keep Ia on the same RG neuron.
s.sel=null;s.selSet=new Set(['Left']);vm.runInContext('unpackSub()',s);const expanded=documentNow();const feedback=expanded.edges.find(e=>e.id==='feedback');assert.equal(feedback.toId,'rgLeft');assert.deepEqual(feedback.toPath,[]);assert.equal(expanded.inputs[0].node_id,'rgLeft');
console.log('PASS: unpacking an old RG wrapper exposes labeled left/right RG blocks, preserves afferent and explicit input bindings, translates internal wire bends, and supports selection-box unpacking.');
