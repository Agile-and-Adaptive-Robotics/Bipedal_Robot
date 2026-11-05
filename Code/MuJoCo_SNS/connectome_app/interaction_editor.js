// Selection and routing are presentation edits; canonical endpoints stay fixed.
let selectedWirePoints=new Set(), wirePointDrag=null, perimeterDrag=null, nodePointDrag=null;
const wirePointKey=(e,k)=>e.id+':'+k;
function selectedPointRecords(){return edges.flatMap(e=>(e.pts||[]).flatMap((p,k)=>selectedWirePoints.has(wirePointKey(e,k))?[{e,k,x:p.x,y:p.y}]:[]));}
function selectRectangle(r,add=false){
 if(!add){selSet.clear();selectedWirePoints.clear();}
 for(const n of nodes){const g=desktopGeom(n);if(!isHiddenNode(n)&&n.x+g.w/2>=r.x0&&n.x-g.w/2<=r.x1&&n.y+g.h/2>=r.y0&&n.y-g.h/2<=r.y1)selSet.add(n.id);}
 for(const e of edges)if(!isHiddenNode(nodeById(e.from))&&!isHiddenNode(nodeById(e.to)))for(const[p,k]of(e.pts||[]).map((p,k)=>[p,k]))if(p.x>=r.x0&&p.x<=r.x1&&p.y>=r.y0&&p.y<=r.y1)selectedWirePoints.add(wirePointKey(e,k));
 sel=selSet.size===1?{kind:'node',id:[...selSet][0]}:null;applySelClasses();redrawEdges();updatePanel();status(`${selSet.size} neurons, ${selectedWirePoints.size} wire points selected`);
}
function perimeterPort(n,p){
 const g=desktopGeom(n),x=p.x-n.x,y=p.y-n.y;
 if(g.shape==='circle'||g.shape==='ellipse')return{angle:Math.atan2(y/(g.h/2),x/(g.w/2))};
 const candidates=[{side:'left',d:Math.abs(x+g.w/2),offset:(y+g.h/2)/g.h},{side:'right',d:Math.abs(x-g.w/2),offset:(y+g.h/2)/g.h},{side:'top',d:Math.abs(y+g.h/2),offset:(x+g.w/2)/g.w},{side:'bottom',d:Math.abs(y-g.h/2),offset:(x+g.w/2)/g.w}];
 const hit=candidates.sort((a,b)=>a.d-b.d)[0];return{side:hit.side,offset:Math.max(0,Math.min(1,hit.offset))};
}
const perimeterAnchor=anchor;
anchor=function(n,tx,ty,port){if(port?.angle!==undefined){const g=desktopGeom(n);return{x:n.x+g.w/2*Math.cos(port.angle),y:n.y+g.h/2*Math.sin(port.angle)};}return perimeterAnchor(n,tx,ty,port);};
function startSelectedPointDrag(ev,e,k){
 if(ev.button!==undefined&&ev.button!==0)return;
 ev.stopPropagation();const key=wirePointKey(e,k);
 if(!selectedWirePoints.has(key)){if(!ev.ctrlKey){selectedWirePoints.clear();selSet.clear();}selectedWirePoints.add(key);}
 wirePointDrag={origin:svgPoint(ev),points:selectedPointRecords(),neurons:[...selSet].map(id=>nodeById(id)).filter(Boolean).map(n=>({n,x:n.x,y:n.y})),moved:false};
 window.addEventListener('mousemove',onMove);window.addEventListener('mouseup',endDrag);redrawEdges();applySelClasses();
}
const routingRenderEdge=renderEdge;
renderEdge=function(e,hm){routingRenderEdge(e,hm);if(!e.el)return;e.el.querySelectorAll('.bend').forEach(h=>h.remove());
 for(const[p,k]of(e.pts||[]).map((p,k)=>[p,k]))if(selectedWirePoints.has(wirePointKey(e,k))||(sel?.kind==='edge'&&sel.id===e.id)){
  const h=dEl('circle',{cx:p.x,cy:p.y,r:6,fill:selectedWirePoints.has(wirePointKey(e,k))?'#CD34B5':'white',stroke:'#A01A91','stroke-width':2,class:'bend'});
  h.addEventListener('mousedown',ev=>startSelectedPointDrag(ev,e,k));h.addEventListener('click',ev=>ev.stopPropagation());e.el.append(h);
 }
 if(sel?.kind==='edge'&&sel.id===e.id){const G=edgeVerts(e);for(const side of['from','to']){const p=side==='from'?G.vs[0]:G.vs.at(-1),h=dEl('rect',{x:p.x-5,y:p.y-5,width:10,height:10,fill:'#FFF',stroke:'#246CB2','stroke-width':2,class:'bend'});
  h.addEventListener('mousedown',ev=>{ev.stopPropagation();perimeterDrag={e,side,moved:false};window.addEventListener('mousemove',onMove);window.addEventListener('mouseup',endDrag);});h.addEventListener('click',ev=>ev.stopPropagation());e.el.append(h);
 }}
};
const routingStartDrag=startDrag;
startDrag=function(ev,id){const n=nodeById(id);nodePointDrag=!ev.altKey&&!ev.shiftKey&&selSet.has(id)?{id,x:n.x,y:n.y,points:selectedPointRecords()}:null;routingStartDrag(ev,id);};
const routingMove=onMove;
onMove=function(ev){
 if(perimeterDrag){const d=perimeterDrag;if(!d.moved)snapshot();d.moved=true;d.e[d.side+'Port']=perimeterPort(nodeById(d.e[d.side]),svgPoint(ev));redrawEdges();return;}
 if(wirePointDrag){const d=wirePointDrag,p=svgPoint(ev),dx=p.x-d.origin.x,dy=p.y-d.origin.y;if(!d.moved&&Math.hypot(dx,dy)<=2)return;if(!d.moved)snapshot();d.moved=true;for(const a of d.points)a.e.pts[a.k]={x:a.x+dx,y:a.y+dy};for(const a of d.neurons){a.n.x=a.x+dx;a.n.y=a.y+dy;a.n.el?.remove();renderNode(a.n,true);}redrawEdges();return;}
 routingMove(ev);
 if(nodePointDrag&&drag?.moved){const d=nodePointDrag,n=nodeById(d.id);for(const a of d.points)a.e.pts[a.k]={x:a.x+n.x-d.x,y:a.y+n.y-d.y};redrawEdges();}
};
const routingEnd=endDrag;
endDrag=function(ev){
 if(perimeterDrag||wirePointDrag){const moved=perimeterDrag?.moved||wirePointDrag?.moved;perimeterDrag=null;wirePointDrag=null;window.removeEventListener('mousemove',onMove);window.removeEventListener('mouseup',endDrag);if(moved){suppressClick=true;persist();}updatePanel();return;}
 if(band?.moved){const r={x0:Math.min(band.x0,band.x1),x1:Math.max(band.x0,band.x1),y0:Math.min(band.y0,band.y1),y1:Math.max(band.y0,band.y1)};bandRect?.remove();bandRect=null;band=null;suppressClick=true;window.removeEventListener('mousemove',onMove);window.removeEventListener('mouseup',endDrag);selectRectangle(r,ev.ctrlKey);return;}
 const moved=drag?.moved;routingEnd(ev);if(moved){suppressClick=true;persist();}nodePointDrag=null;
};
const routingClear=clearSel;
clearSel=function(){selectedWirePoints.clear();routingClear();};
const routingSelectNode=selectNode,routingSelectEdge=selectEdge;
selectNode=function(id){if(!selSet.has(id))selectedWirePoints.clear();routingSelectNode(id);};
selectEdge=function(id){selectedWirePoints.clear();routingSelectEdge(id);};
const routingDelete=deleteSel;
deleteSel=function(){if(selectedWirePoints.size){snapshot();for(const e of edges)e.pts=(e.pts||[]).filter((p,k)=>!selectedWirePoints.has(wirePointKey(e,k)));selectedWirePoints.clear();redrawEdges();persist();if(!selSet.size)return;}routingDelete();};
const routingReload=reloadCurrent;
reloadCurrent=function(spec,keep){selectedWirePoints.clear();routingReload(spec,keep);};
