
"use strict";
// ---------------- type catalog ----------------
const TYPES = {
  "SN-Ia":      {shape:"circle", color:"#FFD700", w:2.5, grp:"Afferents"},
  "SN-II":      {shape:"circle", color:"#FFD700", w:2.5, grp:"Afferents"},
  "SN-Ib":      {shape:"circle", color:"#FFD700", w:2.5, grp:"Afferents"},
  "SN-heel":    {shape:"circle", color:"#FFD700", w:2.5, grp:"Afferents"},
  "SN-toe":     {shape:"circle", color:"#FFD700", w:2.5, grp:"Afferents"},
  "PORT-load":  {shape:"rect",   color:"#FFB14E", w:2,   grp:"Afferents"},
  "IN-V0D":     {shape:"circle", color:"#FA8775", w:2,   grp:"V/C INs"},
  "IN-V0V":     {shape:"circle", color:"#FA8775", w:2,   grp:"V/C INs"},
  "IN-V1":      {shape:"circle", color:"#FA8775", w:2,   grp:"V/C INs"},
  "IN-V2a":     {shape:"circle", color:"#FA8775", w:2,   grp:"V/C INs"},
  "IN-V2b":     {shape:"circle", color:"#FA8775", w:2,   grp:"V/C INs"},
  "IN-V3":      {shape:"circle", color:"#FA8775", w:2,   grp:"V/C INs"},
  "IN-C":       {shape:"circle", color:"#FA8775", w:2,   grp:"V/C INs"},
  "HC-RG-E":    {shape:"circle", color:"#9D02D7", w:2.5, grp:"Rhythm"},
  "HC-RG-F":    {shape:"circle", color:"#9D02D7", w:2.5, grp:"Rhythm"},
  "IN-InE":     {shape:"circle", color:"#FA8775", w:2,   grp:"Rhythm"},
  "IN-InF":     {shape:"circle", color:"#FA8775", w:2,   grp:"Rhythm"},
  "HC-PF-E":    {shape:"circle", color:"#9D02D7", w:2,   grp:"Pattern"},
  "HC-PF-F":    {shape:"circle", color:"#9D02D7", w:2,   grp:"Pattern"},
  "IN-PF":      {shape:"circle", color:"#FA8775", w:2,   grp:"Pattern"},
  "IN-IaIN":    {shape:"circle", color:"#FA8775", w:2,   grp:"Reflex INs"},
  "IN-IbIN":    {shape:"circle", color:"#FA8775", w:2,   grp:"Reflex INs"},
  "IN-IIe":     {shape:"circle", color:"#FA8775", w:2,   grp:"Reflex INs"},
  "IN-IIi":     {shape:"circle", color:"#FA8775", w:2,   grp:"Reflex INs"},
  "IN-Ib+":     {shape:"circle", color:"#FA8775", w:2,   grp:"Reflex INs"},
  "IN-LBIN":    {shape:"circle", color:"#EA5F94", w:2,   grp:"Reflex INs"},
  "IN-KINH":    {shape:"circle", color:"#EA5F94", w:2,   grp:"Reflex INs"},
  "RC":         {shape:"circle", color:"#FA8775", w:2,   grp:"Reflex INs"},
  "MN":         {shape:"circle", color:"#0000FF", w:3,   grp:"Output"},
  "MUSCLE":     {shape:"ellipse",color:"#FFB14E", w:2,   grp:"Output"},
  "SUB":        {shape:"circle", color:"#555555", w:3,   grp:"Subsystems"},
};

const svg = document.getElementById("svg");
const NS = "http://www.w3.org/2000/svg";
const status = m => document.getElementById("status").textContent = m;

// ---------------- world group: pan + zoom (v3.0) ----------------
// Everything drawn lives in #world; the view is translate+scale on it.
// Node coordinates stay in world units for localStorage/JSON/undo, so
// zooming never mutates a model.
const world = document.createElementNS(NS, "g");
world.setAttribute("id", "world");
svg.appendChild(world);
let view = { x: 0, y: 0, k: 1 };
function applyView() {
  world.setAttribute("transform",
    "translate(" + view.x + "," + view.y + ") scale(" + view.k + ")");
  document.getElementById("zoomLabel").textContent =
    Math.round(view.k * 100) + "%";
}
function svgPoint(ev) {   // client coords -> world coords
  const r = svg.getBoundingClientRect();
  return { x: (ev.clientX - r.left - view.x) / view.k,
           y: (ev.clientY - r.top - view.y) / view.k };
}
function zoomAt(cx, cy, f) {
  const k2 = Math.min(4, Math.max(0.05, view.k * f));
  f = k2 / view.k;
  const r = svg.getBoundingClientRect();
  const px = cx - r.left, py = cy - r.top;   // screen px in canvas box
  view.x = px - (px - view.x) * f;
  view.y = py - (py - view.y) * f;
  view.k = k2;
  applyView();
}
function fitView() {
  const vis = nodes.filter(n => !isHiddenNode(n));
  if (!vis.length) { view = { x: 0, y: 0, k: 1 }; applyView(); return; }
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
  vis.forEach(n => { x0 = Math.min(x0, n.x - 60); y0 = Math.min(y0, n.y - 50);
                     x1 = Math.max(x1, n.x + 60); y1 = Math.max(y1, n.y + 50); });
  const r = svg.getBoundingClientRect();
  const k = Math.min(4, Math.max(0.05,
    Math.min(r.width / (x1 - x0), r.height / (y1 - y0))));
  view = { x: (r.width - (x1 - x0) * k) / 2 - x0 * k,
           y: (r.height - (y1 - y0) * k) / 2 - y0 * k, k };
  applyView();
}
svg.addEventListener("wheel", ev => {
  ev.preventDefault();
  zoomAt(ev.clientX, ev.clientY, Math.exp(-ev.deltaY * 0.0012));
}, { passive: false });
document.getElementById("bZoomIn").onclick = () => {
  const r = svg.getBoundingClientRect();
  zoomAt(r.left + r.width / 2, r.top + r.height / 2, 1.25);
};
document.getElementById("bZoomOut").onclick = () => {
  const r = svg.getBoundingClientRect();
  zoomAt(r.left + r.width / 2, r.top + r.height / 2, 1 / 1.25);
};
document.getElementById("bZoomFit").onclick = fitView;
document.getElementById("bZoom100").onclick = () =>
  { view = { x: 0, y: 0, k: 1 }; applyView(); };
// middle-button drag pans; left-drag on empty canvas still rubber-bands
let pan = null, lastMouseWorld = { x: 300, y: 200 };
svg.addEventListener("mousedown", ev => {
  if (ev.button !== 1) return;
  ev.preventDefault();
  pan = { x: ev.clientX, y: ev.clientY, vx: view.x, vy: view.y };
  const up = () => { pan = null;
    window.removeEventListener("mousemove", mv);
    window.removeEventListener("mouseup", up); };
  const mv = e2 => { if (!pan) return;
    view.x = pan.vx + (e2.clientX - pan.x);
    view.y = pan.vy + (e2.clientY - pan.y);
    applyView(); };
  window.addEventListener("mousemove", mv);
  window.addEventListener("mouseup", up);
});
svg.addEventListener("mousemove", ev => { lastMouseWorld = svgPoint(ev); });
applyView();

// ---------- palette (was dropped in the v2.1 rewrite — restored) ----
const pal = document.getElementById("palette");
const palGroups = {};
Object.entries(TYPES).forEach(([t, d]) => {
  (palGroups[d.grp] = palGroups[d.grp] || []).push([t, d]);
});
Object.entries(palGroups).forEach(([g, items]) => {
  const h = document.createElement("h3");
  h.textContent = g;
  pal.appendChild(h);
  items.forEach(([t, d]) => {
    const el = document.createElement("div");
    el.className = "pitem";
    el.style.borderColor = d.color;
    el.textContent = t;
    el.onclick = () => addNode(t, lastMouseWorld.x, lastMouseWorld.y);
    pal.appendChild(el);
  });
});

// ---------------- tabs / models ----------------
let models = [{ name: "Ben", data: { nodes: [], edges: [] } }];
let activeTab = 0;
try {
  const saved = JSON.parse(localStorage.getItem("cme_models_v1"));
  if (saved && saved.models && saved.models.length) {
    models = saved.models;
    activeTab = Math.min(saved.active || 0, models.length - 1);
  }
} catch (e) {}

function cur() { return models[activeTab].data; }
function persist() {
  // v3.0 fix: sync the CURRENT tab from the canvas before saving.
  // (v2.2 only synced on tab-switch, so the active tab's localStorage
  // copy went stale until you switched away -- closing the browser
  // right after editing lost that tab's latest changes.)
  // v3.3: persist the FLATTENED model (nested subsystem write-back).
  try {
    models[activeTab].data = modelSpec();
    localStorage.setItem("cme_models_v1", JSON.stringify(
      { models, active: activeTab }));
  } catch (e) {}
}
function renderTabs() {
  const t = document.getElementById("tabs");
  t.innerHTML = "";
  models.forEach((m, i) => {
    const el = document.createElement("span");
    el.className = "tab" + (i === activeTab ? " active" : "");
    el.textContent = m.name;
    el.onclick = () => switchTab(i);
    el.ondblclick = () => {
      const nm = prompt("Model name:", m.name);
      if (nm) { m.name = nm; renderTabs(); persist(); }
    };
    const x = document.createElement("span");
    x.className = "x"; x.textContent = "×";
    x.onclick = ev => {
      ev.stopPropagation();
      if (models.length === 1) return;
      while (ctxStack.length) exitSub();
      models.splice(i, 1);
      activeTab = Math.min(activeTab, models.length - 1);
      undoStack.length = 0; redoStack.length = 0;
      loadIntoCanvas(cur()); renderTabs(); persist();
    };
    el.appendChild(x);
    t.appendChild(el);
  });
  const b = document.createElement("button");
  b.id = "newtab"; b.textContent = "+";
  b.onclick = () => {
    const nm = prompt("New model name:", "model " + (models.length + 1));
    if (!nm) return;
    while (ctxStack.length) exitSub();
    models[activeTab].data = serialize();
    models.push({ name: nm, data: { nodes: [], edges: [] } });
    activeTab = models.length - 1;
    undoStack.length = 0; redoStack.length = 0;
    loadIntoCanvas(cur()); renderTabs(); persist();
  };
  t.appendChild(b);
}
function switchTab(i) {
  if (i === activeTab) return;
  while (ctxStack.length) exitSub();   // never switch from inside a sub
  models[activeTab].data = serialize();
  activeTab = i;
  // v3.0 fix: undo history is per-tab in spirit — never let one tab's
  // snapshots restore into another (the old global stack did exactly
  // that once tabs multiplied).
  undoStack.length = 0; redoStack.length = 0;
  loadIntoCanvas(cur());
  renderTabs(); persist();
}

// ---------------- state ----------------
let nodes = [], edges = [];
let sel = null, linkFrom = null, drag = null, linkDrag = null;
// v3.2 layers: nodes may carry a `grp` id (stamped by the walker
// templates); nodes WITHOUT one fall back to type buckets, so every
// tab gets a tree. hiddenG = layer ids hidden on canvas (persisted as
// data.hidden — a VIEW state that survives tabs/reload); grpNames =
// pretty labels from the spec's `groups` dict; selGrp/openGrps are
// pure tree-UI state.
let hiddenG = new Set();
let grpNames = {};
let selGrp = null;
const openGrps = new Set();
const gidOf = n => n.grp || ("type:" + n.type);
const isHiddenNode = n => hiddenG.has(gidOf(n));
// v3.3 subsystems: a node may carry a nested spec (node.sub =
// {nodes, edges, ...}); double-click ENTERS it -- the current context
// is pushed on ctxStack and the sub's nodes/edges become the canvas,
// fully editable. F4 (or double-click empty canvas) writes the sub
// back into the host node and pops the context. modelSpec() flattens
// the stack so persist/export always save the true nested model.
let ctxStack = [];
// v2.2 group selection: selSet = node ids highlighted + moved
// together; band = rubber-band gesture in progress; suppressClick
// eats the click event that follows a completed drag gesture so it
// cannot clear/replace the selection the drag just produced.
let selSet = new Set(), band = null, bandRect = null;
let suppressClick = false;
let uid = 1;
let undoStack = [], redoStack = [];

// svgPoint() now lives with the world/view code at the top of the
// script (v3.0): it converts client coords to WORLD coords through
// the current zoom/pan. All gestures below use it.

function specOf(ns, es, gn, hid) {
  // serialize an arbitrary (nodes, edges) context -- also used for
  // nested subsystem write-back on exitSub/modelSpec. Dangling edges
  // (endpoint not in ns) are skipped, not fatal.
  const byId = new Map(ns.map(n => [n.id, n]));
  return JSON.parse(JSON.stringify({ nodes: ns.map(n => ({
    id: n.id, type: n.type, label: n.label, x: Math.round(n.x),
    y: Math.round(n.y), grp: n.grp || undefined,
    sub: n.sub || undefined })), edges: es.filter(e =>
      byId.has(e.from) && byId.has(e.to)).map(e => ({
    from: byId.get(e.from).label, to: byId.get(e.to).label,
    sign: e.sign, gain: e.gain, tag: e.tag,
    pts: (e.pts || []).map(p => ({
      x: Math.round(p.x * 10) / 10, y: Math.round(p.y * 10) / 10 })) })),
    groups: Object.keys(gn).length ? gn : undefined,
    hidden: hid.size ? [...hid] : undefined }));
}
function serialize() {
  return specOf(nodes, edges, grpNames, hiddenG);
}
// flatten the context stack: write each frame's live nodes/edges back
// into its host node's sub, innermost first, then serialize the
// outermost context -- the TRUE model for persist/export.
function modelSpec() {
  for (let i = ctxStack.length - 1; i >= 0; i--) {
    const f = ctxStack[i];
    f.hostNode.sub = specOf(f === ctxStack[ctxStack.length - 1]
      ? nodes : f.nodes,
      f === ctxStack[ctxStack.length - 1] ? edges : f.edges,
      f.grpNames, f.hiddenG);
  }
  return serialize();
}
function snapshot() {
  undoStack.push(JSON.stringify(serialize()));
  if (undoStack.length > 100) undoStack.shift();
  redoStack.length = 0;
  persist();
}
function undo() {
  if (!undoStack.length) { status("nothing to undo"); return; }
  redoStack.push(JSON.stringify(serialize()));
  loadIntoCanvas(JSON.parse(undoStack.pop()), true);
  status("undo");
}
function redo() {
  if (!redoStack.length) { status("nothing to redo"); return; }
  undoStack.push(JSON.stringify(serialize()));
  loadIntoCanvas(JSON.parse(redoStack.pop()), true);
  status("redo");
}

function nodeById(id) { return nodes.find(n => n.id === id); }
function edgeById(id) { return edges.find(e => e.id === id); }

// ---------------- node rendering ----------------
function nodeRadius() { return 18; }
function anchor(nd, tx, ty) {
  // point on the node body boundary toward (tx,ty)
  const dx = tx - nd.x, dy = ty - nd.y;
  const L = Math.hypot(dx, dy) || 1;
  const r = nodeRadius() + 2;
  return { x: nd.x + dx / L * r, y: nd.y + dy / L * r };
}
function renderNode(nd, keepSel) {
  if (isHiddenNode(nd)) {          // hidden layer: no DOM at all
    if (nd.el) { nd.el.remove(); nd.el = null; }
    return;
  }
  const d = TYPES[nd.type];
  const g = document.createElementNS(NS, "g");
  g.setAttribute("class", "node");
  g.dataset.id = nd.id;
  let shape;
  if (d.shape === "rect") {
    shape = document.createElementNS(NS, "rect");
    shape.setAttribute("x", nd.x - 26); shape.setAttribute("y", nd.y - 13);
    shape.setAttribute("width", 52); shape.setAttribute("height", 26);
  } else if (d.shape === "ellipse") {
    shape = document.createElementNS(NS, "ellipse");
    shape.setAttribute("cx", nd.x); shape.setAttribute("cy", nd.y);
    shape.setAttribute("rx", 30); shape.setAttribute("ry", 13);
  } else {
    shape = document.createElementNS(NS, "circle");
    shape.setAttribute("cx", nd.x); shape.setAttribute("cy", nd.y);
    shape.setAttribute("r", 18);
  }
  shape.setAttribute("fill", "#fff");
  shape.setAttribute("stroke", d.color);
  shape.setAttribute("stroke-width", d.w);
  g.appendChild(shape);
  if (nd.sub) {
    // subsystem cue: dashed outer ring (AnimatLab-style module)
    const ring = document.createElementNS(NS, "circle");
    ring.setAttribute("cx", nd.x); ring.setAttribute("cy", nd.y);
    ring.setAttribute("r", 24);
    ring.setAttribute("fill", "none");
    ring.setAttribute("stroke", d.color);
    ring.setAttribute("stroke-width", "1.2");
    ring.setAttribute("stroke-dasharray", "3,3");
    ring.setAttribute("class", "subring");
    g.appendChild(ring);
  }
  const txt = document.createElementNS(NS, "text");
  txt.setAttribute("x", nd.x); txt.setAttribute("y", nd.y - 24);
  txt.setAttribute("text-anchor", "middle");
  txt.textContent = nd.label;
  g.appendChild(txt);
  g.addEventListener("mousedown", ev => startDrag(ev, nd.id));
  g.addEventListener("dblclick", ev => {
    ev.stopPropagation();
    if (nd.sub) enterSub(nd);    // ENTER the subsystem (view its parts)
    else status("no subsystem inside " + nd.label +
                " — select it and Ctrl+G to pack one");
  });
  g.addEventListener("click", ev => {
    ev.stopPropagation();
    if (suppressClick) { suppressClick = false; return; }
    if (drag && drag.moved) return;
    onClickNode(nd.id, ev.shiftKey || ev.ctrlKey, ev.altKey);
  });
  nd.el = g;
  world.appendChild(g);
  if (keepSel && selSet.has(nd.id))
    g.classList.add("sel");
}

function addNode(type, x, y, label, silent, grp, sub) {
  const nd = { id: "n" + (uid++), type, x, y,
               label: label || (type + "_" + uid), grp };
  if (sub) nd.sub = sub;
  nodes.push(nd);
  renderNode(nd);      // renderNode honors hiddenG via nd.grp
  if (!silent) { snapshot(); renderTree(); }
  return nd;
}

// ---------------- edges ----------------
// An edge is a polyline: node A anchor -> e.pts[] bend points -> node
// B anchor. Bend points (ALT+click a segment) live in world coords,
// drag by their pink handles, and serialize as "pts" so undo, tabs,
// copy/paste and JSON export keep the routing. Where two edges cross,
// the later-created one hops in a semicircle (circuit-diagram style).
const HOP_R = 7;        // hop radius (world units)
const HIT_W = 12;       // invisible selection line width
const HIT_INSET = 7;    // hit lines stop short of the node bodies
const r2 = v => Math.round(v * 100) / 100;

function edgeVerts(e) {
  const a = nodeById(e.from), b = nodeById(e.to);
  if (!a || !b) return null;
  const pts = e.pts || [];
  const first = pts.length ? pts[0] : { x: b.x, y: b.y };
  const last = pts.length ? pts[pts.length - 1] : { x: a.x, y: a.y };
  const vs = [anchor(a, first.x, first.y),
              ...pts.map(p => ({ x: p.x, y: p.y })),
              anchor(b, last.x, last.y)];
  const segs = [];
  for (let i = 0; i < vs.length - 1; i++)
    segs.push({ x1: vs[i].x, y1: vs[i].y,
                x2: vs[i + 1].x, y2: vs[i + 1].y });
  return { vs, segs };
}

// proper interior crossing of two segments (null when parallel,
// collinear, or touching only at an endpoint)
function segXint(s, o) {
  const d1x = s.x2 - s.x1, d1y = s.y2 - s.y1;
  const d2x = o.x2 - o.x1, d2y = o.y2 - o.y1;
  const den = d1x * d2y - d1y * d2x;
  if (Math.abs(den) < 1e-9) return null;
  const t = ((o.x1 - s.x1) * d2y - (o.y1 - s.y1) * d2x) / den;
  const u = ((o.x1 - s.x1) * d1y - (o.y1 - s.y1) * d1x) / den;
  if (t <= 1e-4 || t >= 1 - 1e-4 || u <= 1e-4 || u >= 1 - 1e-4)
    return null;
  return { x: s.x1 + t * d1x, y: s.y1 + t * d1y };
}

// point at fraction f of the total polyline length (sign marker/label)
function pointAtFrac(segs, f) {
  let L = 0;
  segs.forEach(s => L += Math.hypot(s.x2 - s.x1, s.y2 - s.y1));
  let target = L * f;
  for (let i = 0; i < segs.length; i++) {
    const s = segs[i];
    const l = Math.hypot(s.x2 - s.x1, s.y2 - s.y1);
    if (target <= l || i === segs.length - 1) {
      const k = l ? target / l : 0;
      return { x: s.x1 + (s.x2 - s.x1) * k,
               y: s.y1 + (s.y2 - s.y1) * k };
    }
    target -= l;
  }
}

// which edges hop where: pairwise AABB cull, then exact segment
// tests. The edge with the LARGER uid (created later) is the hopper,
// so the choice is stable no matter what order things render in.
function computeHops() {
  const geos = edges.map(e => {
    if (isHiddenNode(nodeById(e.from)) || isHiddenNode(nodeById(e.to)))
      return null;                   // hidden wires make no crossings
    const g = edgeVerts(e);
    if (!g) return null;
    let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
    g.segs.forEach(s => {
      x0 = Math.min(x0, s.x1, s.x2); y0 = Math.min(y0, s.y1, s.y2);
      x1 = Math.max(x1, s.x1, s.x2); y1 = Math.max(y1, s.y1, s.y2);
    });
    return { g, x0, y0, x1, y1 };
  });
  const hops = new Map();   // edge.id -> { segIndex: [ { t } ] }
  const add = (id, seg, t) => {
    if (!hops.has(id)) hops.set(id, {});
    const m = hops.get(id);
    (m[seg] = m[seg] || []).push({ t });
  };
  for (let i = 0; i < edges.length; i++) {
    const gi = geos[i];
    if (!gi) continue;
    for (let j = i + 1; j < edges.length; j++) {
      const gj = geos[j];
      if (!gj || gj.x0 > gi.x1 || gj.x1 < gi.x0 ||
          gj.y0 > gi.y1 || gj.y1 < gi.y0) continue;
      const w = +(edges[i].id.slice(1)) > +(edges[j].id.slice(1))
        ? { e: edges[i], g: gi.g, o: gj.g }
        : { e: edges[j], g: gj.g, o: gi.g };
      for (let a = 0; a < w.g.segs.length; a++)
        for (let b = 0; b < w.o.segs.length; b++) {
          const p = segXint(w.g.segs[a], w.o.segs[b]);
          if (!p) continue;
          const s = w.g.segs[a];
          const L = Math.hypot(s.x2 - s.x1, s.y2 - s.y1) || 1;
          add(w.e.id, a, Math.hypot(p.x - s.x1, p.y - s.y1) / L);
        }
    }
  }
  return hops;
}

function renderEdge(e, hopFor) {
  const G = edgeVerts(e);
  if (!G) return;
  if (isHiddenNode(nodeById(e.from)) || isHiddenNode(nodeById(e.to)))
    return;                          // an edge into a hidden layer goes too
  const segs = G.segs;
  const g = document.createElementNS(NS, "g");
  g.setAttribute("class", "edge");
  g.dataset.id = e.id;
  // visible wire: straight runs + one semicircular arc per crossing
  const hm = hopFor || {};
  let d = "M" + r2(G.vs[0].x) + " " + r2(G.vs[0].y);
  segs.forEach((s, i) => {
    const L = Math.hypot(s.x2 - s.x1, s.y2 - s.y1) || 1e-9;
    const ux = (s.x2 - s.x1) / L, uy = (s.y2 - s.y1) / L;
    let done = 0;
    ((hm[i] || []).slice().sort((p, q) => p.t - q.t)).forEach(h => {
      const c = h.t * L;
      if (c - HOP_R <= done + 1 || c + HOP_R > L - 1) return;
      const c0 = c - HOP_R, c1 = c + HOP_R;
      d += " L" + r2(s.x1 + c0 * ux) + " " + r2(s.y1 + c0 * uy) +
           " A" + HOP_R + " " + HOP_R + " 0 0 1 " +
           r2(s.x1 + c1 * ux) + " " + r2(s.y1 + c1 * uy);
      done = c1;
    });
    d += " L" + r2(s.x2) + " " + r2(s.y2);
  });
  const path = document.createElementNS(NS, "path");
  path.setAttribute("d", d);
  path.setAttribute("fill", "none");
  path.setAttribute("stroke", e.sign === "inh" ? "#B2182B" : "#0072B2");
  path.setAttribute("stroke-width", 2);
  path.setAttribute("pointer-events", "none");
  g.appendChild(path);
  const mp = pointAtFrac(segs, 0.55);
  if (e.sign === "inh") {
    const dot = document.createElementNS(NS, "circle");
    dot.setAttribute("cx", mp.x); dot.setAttribute("cy", mp.y);
    dot.setAttribute("r", 6); dot.setAttribute("fill", "#1a1a1a");
    dot.setAttribute("pointer-events", "none");
    g.appendChild(dot);
  } else {
    const tri = document.createElementNS(NS, "polygon");
    tri.setAttribute("points", `${mp.x},${mp.y - 8} ${mp.x - 7},${mp.y + 6} ` +
      `${mp.x + 7},${mp.y + 6}`);
    tri.setAttribute("fill", "#fff");
    tri.setAttribute("stroke", "#1a1a1a");
    tri.setAttribute("stroke-width", "1.4");
    tri.setAttribute("pointer-events", "none");
    g.appendChild(tri);
  }
  const last = segs[segs.length - 1];
  const Ll = Math.hypot(last.x2 - last.x1, last.y2 - last.y1) || 1;
  const ux = (last.x2 - last.x1) / Ll, uy = (last.y2 - last.y1) / Ll;
  const head = document.createElementNS(NS, "polygon");
  head.setAttribute("points",
    `${r2(last.x2)},${r2(last.y2)} ` +
    `${r2(last.x2 - ux * 9 - uy * 4)},${r2(last.y2 - uy * 9 + ux * 4)} ` +
    `${r2(last.x2 - ux * 9 + uy * 4)},${r2(last.y2 - uy * 9 - ux * 4)}`);
  head.setAttribute("fill", "#1a1a1a");
  head.setAttribute("pointer-events", "none");
  g.appendChild(head);
  const lab = document.createElementNS(NS, "text");
  lab.setAttribute("x", mp.x); lab.setAttribute("y", mp.y - 10);
  lab.setAttribute("text-anchor", "middle");
  lab.setAttribute("font-size", "10");
  lab.setAttribute("pointer-events", "none");
  lab.textContent = (e.tag ? e.tag + " " : "") + e.gain;
  g.appendChild(lab);
  // invisible wide hit lines: one per segment, pulled in from the
  // ends so node rims stay grabbable. The seg index is how ALT+click
  // names which segment gets the bend point.
  segs.forEach((s, i) => {
    const L = Math.hypot(s.x2 - s.x1, s.y2 - s.y1);
    if (L < HIT_INSET * 2 + 2) return;
    const ux = (s.x2 - s.x1) / L, uy = (s.y2 - s.y1) / L;
    const hl = document.createElementNS(NS, "line");
    hl.setAttribute("x1", r2(s.x1 + HIT_INSET * ux));
    hl.setAttribute("y1", r2(s.y1 + HIT_INSET * uy));
    hl.setAttribute("x2", r2(s.x2 - HIT_INSET * ux));
    hl.setAttribute("y2", r2(s.y2 - HIT_INSET * uy));
    hl.setAttribute("stroke", "transparent");
    hl.setAttribute("stroke-width", HIT_W);
    hl.setAttribute("pointer-events", "stroke");
    hl.dataset.seg = i;
    g.appendChild(hl);
  });
  // bend-point handles (draggable)
  (e.pts || []).forEach((p, k) => {
    const h = document.createElementNS(NS, "circle");
    h.setAttribute("cx", p.x); h.setAttribute("cy", p.y);
    h.setAttribute("r", 4.5);
    h.setAttribute("fill", "#fff");
    h.setAttribute("stroke", "#CD34B5");
    h.setAttribute("stroke-width", "2");
    h.setAttribute("class", "bend");
    h.addEventListener("mousedown", ev => {
      ev.stopPropagation();
      startBendDrag(ev, e, k);
    });
    g.appendChild(h);
  });
  g.addEventListener("click", ev => {
    ev.stopPropagation();
    if (suppressClick) { suppressClick = false; return; }
    if (ev.altKey && ev.target.tagName === "line") {
      addBend(e, +(ev.target.dataset.seg));
      return;
    }
    selectEdge(e.id);
  });
  e.el = g;
  world.appendChild(g);
}
function redrawEdges() {
  const hops = computeHops();
  edges.forEach(e => { if (e.el) e.el.remove(); e.el = null; });
  edges.forEach(e => renderEdge(e, hops.get(e.id)));
}

// ALT+click a segment = ONE bend point at that segment's midpoint.
// A held ALT cannot machine-gun points: creation is click-driven
// (one per physical click) with a cooldown on top, and each further
// bend needs a fresh ALT+click on one of the new segments.
// far past, so the very first ALT+click after load is never eaten
let lastBendAt = -1e9;
function addBend(e, seg) {
  const now = performance.now();
  if (now - lastBendAt < 250) return;
  lastBendAt = now;
  const G = edgeVerts(e);
  if (!G || !(seg >= 0) || seg >= G.segs.length) return;
  const s = G.segs[seg];
  snapshot();
  if (!e.pts) e.pts = [];
  e.pts.splice(seg, 0, { x: (s.x1 + s.x2) / 2, y: (s.y1 + s.y2) / 2 });
  redrawEdges();
  selectEdge(e.id);
  status("bend point added — drag the pink handle to route the wire; " +
         "ALT+click any segment for another");
}

let bendDrag = null;
function startBendDrag(ev, e, k) {
  bendDrag = { e, k, moved: false };
  window.addEventListener("mousemove", onMove);
  window.addEventListener("mouseup", endDrag);
}

// ---------------- selection / panel ----------------
// One place decides the visible highlight: every node in selSet
// plus (if the primary selection is an edge) that edge.
function applySelClasses() {
  nodes.forEach(n => n.el &&
    n.el.classList.toggle("sel", selSet.has(n.id)));
  edges.forEach(e => e.el &&
    e.el.classList.toggle("sel", !!(sel && sel.kind === "edge" &&
                                    sel.id === e.id)));
}
function clearSel() {
  sel = null;
  selSet.clear();
  applySelClasses();
}
function selectNode(id) {
  sel = { kind: "node", id };
  selSet = new Set([id]);
  applySelClasses();
  showNodePanel(nodeById(id));
}
function selectEdge(id) {
  sel = { kind: "edge", id };
  selSet.clear();
  applySelClasses();
  showEdgePanel(edgeById(id));
}
// SHIFT+click: toggle one node in/out of the group selection.
function toggleNodeSel(id) {
  if (selSet.has(id)) {
    selSet.delete(id);
    if (sel && sel.kind === "node" && sel.id === id) {
      const next = selSet.size ? [...selSet][selSet.size - 1] : null;
      sel = next ? { kind: "node", id: next } : null;
    }
  } else {
    selSet.add(id);
    if (!(sel && sel.kind === "node" && selSet.has(sel.id))) {
      sel = { kind: "node", id };
    }
  }
  applySelClasses();
  updatePanel();
  status(selSet.size ? selSet.size + " node(s) selected — drag one" +
         " to move the group" : "selection cleared");
}
// Panel follows the primary selection; a group of >1 shows a hint.
function updatePanel() {
  if (sel && sel.kind === "node" && nodeById(sel.id)) {
    showNodePanel(nodeById(sel.id));
  } else if (sel && sel.kind === "edge" && edgeById(sel.id)) {
    showEdgePanel(edgeById(sel.id));
  } else if (selSet.size > 1) {
    document.getElementById("pTitle").textContent =
      selSet.size + " nodes selected";
    document.getElementById("pBody").innerHTML =
      `<div class="sub">Drag any selected node to move the whole ` +
      `group (relative positions kept). SHIFT+click a node to add ` +
      `or remove it. Click empty canvas or press Esc to clear.</div>`;
  } else {
    document.getElementById("pTitle").textContent = "Nothing selected";
    document.getElementById("pBody").innerHTML =
      `<span class="sub">Click a node or synapse.</span>`;
  }
}
function showNodePanel(nd) {
  document.getElementById("pTitle").textContent = nd.type;
  document.getElementById("pBody").innerHTML =
    `<label>label</label><input id="pLabel" value="${nd.label}">` +
    `<div class="sub">SHIFT/ALT-drag to another node = synapse; ` +
    `SHIFT+click = group-select toggle; ALT-click A then B = ` +
    `two-click link.</div>`;
  document.getElementById("pLabel").addEventListener("change", ev => {
    snapshot();
    nd.label = ev.target.value;
    nd.el.remove(); renderNode(nd, true); redrawEdges();
  });
}
function showEdgePanel(e) {
  document.getElementById("pTitle").textContent = "Synapse";
  document.getElementById("pBody").innerHTML =
    `<label>sign</label><select id="pSign">
       <option value="exc"${e.sign === "exc" ? " selected" : ""}>
       excitatory</option>
       <option value="inh"${e.sign === "inh" ? " selected" : ""}>
       inhibitory</option></select>
     <label>gain (this synapse only)</label>
     <input id="pGain" type="number" step="0.05" value="${e.gain}">
     <label>tuning tag (group synapses that share one gain)</label>
     <input id="pTag" value="${e.tag}" placeholder="e.g. ib_inh">
     <div class="sub" style="margin-top:6px">ALT+click the wire = bend
     point (drag the pink handle to reroute).</div>
     ${e.pts && e.pts.length
       ? `<button id="pStraight">Remove ${e.pts.length} bend ` +
         `point(s)</button>` : ""}`;
  const upd = () => {
    snapshot();
    e.sign = document.getElementById("pSign").value;
    e.gain = parseFloat(document.getElementById("pGain").value) || 0;
    e.tag = document.getElementById("pTag").value;
    redrawEdges();
  };
  ["pSign", "pGain", "pTag"].forEach(id =>
    document.getElementById(id).addEventListener("change", upd));
  const ps = document.getElementById("pStraight");
  if (ps) ps.addEventListener("click", () => {
    snapshot();
    e.pts = [];
    redrawEdges();
    showEdgePanel(e);
    status("edge straightened");
  });
}

// ---------------- layers tree (v3.2) ----------------
// Left-hand Layers panel: one row per layer (explicit `grp`, else the
// node's TYPE so every model gets a tree). Click selects, double-click
// drills the view to that layer, F2 (or the eye dot) hides/shows it.
function groupsOf() {
  const m = new Map();
  nodes.forEach(n => {
    const gid = gidOf(n);
    if (!m.has(gid))
      m.set(gid, { label: grpNames[gid] ||
                   (n.grp ? n.grp : n.type), ids: [] });
    m.get(gid).ids.push(n.id);
  });
  return m;
}
function renderTree() {
  const t = document.getElementById("tree");
  t.innerHTML = "";
  groupsOf().forEach((ginfo, gid) => {
    const wrap = document.createElement("div");
    wrap.className = "tgrp" + (openGrps.has(gid) ? " open" : "");
    const row = document.createElement("div");
    row.className = "trow" + (hiddenG.has(gid) ? " hid" : "") +
                    (selGrp === gid ? " selp" : "");
    const tw = document.createElement("span");
    tw.className = "tw";
    tw.textContent = openGrps.has(gid) ? "▼" : "▶";
    const lbl = document.createElement("span");
    lbl.className = "lbl";
    lbl.textContent = ginfo.label + "  (" + ginfo.ids.length + ")";
    const eye = document.createElement("span");
    eye.className = "eye";
    eye.textContent = hiddenG.has(gid) ? "○" : "●";
    eye.title = "hide/show layer (F2)";
    row.append(tw, lbl, eye);
    row.onclick = () => { selGrp = gid; openGrps.add(gid); renderTree(); };
    row.ondblclick = () => drillTo(gid);
    eye.onclick = ev => { ev.stopPropagation(); toggleGrpHidden(gid); };
    wrap.appendChild(row);
    const kids = document.createElement("div");
    kids.className = "tkids";
    ginfo.ids.forEach(id => {
      const n = nodeById(id);
      if (!n) return;
      const k = document.createElement("div");
      k.className = "tkid";
      k.textContent = n.label;
      k.title = n.label;
      k.onclick = () => { selectNode(id);
        const r = document.getElementById("svg").getBoundingClientRect();
        zoomAt(r.left + n.x * view.k + view.x,
               r.top + n.y * view.k + view.y, 1); };
      kids.appendChild(k);
    });
    wrap.appendChild(kids);
    t.appendChild(wrap);
  });
}
// zoom the view so the layer fills the canvas ("take me through layers")
function drillTo(gid) {
  const ginfo = groupsOf().get(gid);
  if (!ginfo) return;
  selGrp = gid;
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
  let any = false;
  ginfo.ids.forEach(id => {
    const n = nodeById(id);
    if (!n) return;
    any = true;
    x0 = Math.min(x0, n.x - 60); y0 = Math.min(y0, n.y - 50);
    x1 = Math.max(x1, n.x + 60); y1 = Math.max(y1, n.y + 50);
  });
  if (!any) return;
  const r = svg.getBoundingClientRect();
  const k = Math.min(4, Math.max(0.05,
    Math.min(r.width / (x1 - x0), r.height / (y1 - y0))));
  view = { x: (r.width - (x1 - x0) * k) / 2 - x0 * k,
           y: (r.height - (y1 - y0) * k) / 2 - y0 * k, k };
  applyView();
  openGrps.add(gid);
  renderTree();
  status("layer: " + ginfo.label + " — double-click empty canvas to " +
         "zoom back out");
}
function toggleGrpHidden(gid) {
  if (hiddenG.has(gid)) hiddenG.delete(gid);
  else hiddenG.add(gid);
  sel = null; selSet.clear();     // selections may live in that layer
  loadIntoCanvas(serialize(), true);   // re-render with the layer gone
  renderTree();
  status(hiddenG.has(gid) ? "layer hidden (F2 brings it back)" :
         "layer shown");
}

// ---------------- subsystems (v3.3) ----------------
// AnimatLab-style containment: a node with a nested spec (node.sub)
// can be ENTERED -- its constituent nodes/edges become the canvas,
// fully editable; F4 or double-click on empty canvas writes the sub
// back into the host node and pops the parent context.
function materialize(spec) {
  const ns = [], es = [], byLabel = {};
  (spec.nodes || []).forEach(nd => {
    const n2 = { id: "n" + (uid++), type: nd.type, x: nd.x, y: nd.y,
                 label: nd.label, grp: nd.grp };
    if (nd.sub) n2.sub = nd.sub;
    if (nd.sub) {} // (nested subs ride along as raw specs)
    ns.push(n2);
    byLabel[nd.label] = n2;
  });
  (spec.edges || spec.synapses || []).forEach(sn => {
    if (byLabel[sn.from] && byLabel[sn.to]) {
      es.push({ id: "e" + (uid++), from: byLabel[sn.from].id,
                to: byLabel[sn.to].id, sign: sn.sign || "exc",
                gain: sn.gain === undefined ? 0.5 : sn.gain,
                tag: sn.tag || "",
                pts: (sn.pts || []).map(p => ({ x: +p.x, y: +p.y }))
                  .filter(p => isFinite(p.x) && isFinite(p.y)) });
    }
  });
  return { nodes: ns, edges: es };
}
function clearCanvasDom() {
  // remove ALL node/edge elements from the world group — the current
  // arrays may not reference everything (a just-exited subsystem's
  // elements are orphans; the 2026-09-25 tab-leak bug was exactly
  // this: only referenced els were removed, so entered-subs' DOM
  // floated over every other tab forever)
  world.querySelectorAll(".node, .edge").forEach(el => el.remove());
}
function updateCrumb() {
  const c = document.getElementById("crumb");
  if (!ctxStack.length) { c.textContent = ""; return; }
  c.innerHTML = "<button id='crumbBack' style='font-size:11px;" +
    "padding:1px 6px;cursor:pointer'>◀ back out (F4)</button> " +
    "<b>" + models[activeTab].name + "</b> ▸ " +
    ctxStack.map(f => "<span class='loc'>" + f.hostNode.label +
      "</span>").join(" ▸ ");
  document.getElementById("crumbBack").onclick = () => exitSub();
  [...c.querySelectorAll(".loc")].forEach((el, i) => {
    el.onclick = () => { while (ctxStack.length > i + 1) exitSub(); };
  });
}
function enterSub(node) {
  ctxStack.push({ hostNode: node, nodes, edges, view: { ...view },
                  undo: undoStack, redo: redoStack,
                  grpNames, hiddenG });
  const m = materialize(node.sub);
  nodes = m.nodes; edges = m.edges;
  sel = null; selSet.clear(); linkFrom = null;
  hiddenG = new Set();
  grpNames = node.sub.groups || {};
  undoStack = []; redoStack = [];
  clearCanvasDom();
  nodes.forEach(n => renderNode(n, false));
  redrawEdges(); renderTree(); updateCrumb(); fitView(); persist();
  status("ENTERED subsystem: " + node.label + " — edit its parts, " +
         "F4 to back out");
}
function exitSub() {
  if (!ctxStack.length) { status("not inside a subsystem"); return; }
  const f = ctxStack[ctxStack.length - 1];
  f.hostNode.sub = specOf(nodes, edges, grpNames, hiddenG);
  nodes = f.nodes; edges = f.edges; view = f.view;
  undoStack = f.undo; redoStack = f.redo;
  grpNames = f.grpNames; hiddenG = f.hiddenG;
  sel = null; selSet.clear(); linkFrom = null;
  ctxStack.pop();
  clearCanvasDom();
  nodes.forEach(n => renderNode(n, false));
  redrawEdges(); renderTree(); updateCrumb(); applyView(); persist();
  status("backed out of subsystem — edits saved into " + f.hostNode.label);
}
// pack the selected nodes (+ synapses between them) into one SUB node
function packSub() {
  if (selSet.size < 1) { status("select nodes to pack first"); return; }
  const ids = new Set(selSet);
  const packed = nodes.filter(n => ids.has(n.id));
  const inner = edges.filter(e => ids.has(e.from) && ids.has(e.to));
  if (packed.length < 1) return;
  const nm = prompt("Subsystem name:",
                    "SUB_" + (uid + 1));
  if (!nm) return;
  snapshot();
  const local = new Map(packed.map(n => [n.id, n]));
  const sub = specOfLocal(packed, inner);
  let cx = 0, cy = 0;
  packed.forEach(n => { cx += n.x; cy += n.y; });
  cx /= packed.length; cy /= packed.length;
  packed.forEach(n => { n.el && n.el.remove(); });
  const subNode = { id: "n" + (uid++), type: "SUB", x: cx, y: cy,
                    label: uniqueLabel(nm), sub };
  // boundary edges: RETARGET to the subsystem node (AnimatLab-style
  // ports) instead of dropping -- wires stay connected, none hang
  edges = edges.map(e => {
    if (ids.has(e.from) && !ids.has(e.to)) {
      const r = { ...e, from: subNode.id };
      return r;
    }
    if (ids.has(e.to) && !ids.has(e.from)) {
      const r = { ...e, to: subNode.id };
      return r;
    }
    return e;
  });
  edges = edges.filter(e => !(ids.has(e.from) && ids.has(e.to)));
  nodes = nodes.filter(n => !ids.has(n.id));
  nodes.push(subNode);
  renderNode(subNode);
  redrawEdges(); renderTree(); persist();
  sel = null; selSet = new Set([subNode.id]);
  applySelClasses(); updatePanel();
  status("packed " + packed.length + " node(s) into " + subNode.label +
         " — double-click to enter; boundary wires re-attach to it");
}
function unpackSub() {
  const node = sel && sel.kind === "node" && nodeById(sel.id);
  if (!node || !node.sub) { status("select a subsystem node to unpack"); return; }
  snapshot();
  const m = materialize(node.sub);
  let cx = 0, cy = 0;
  m.nodes.forEach(n => { cx += n.x; cy += n.y; });
  cx /= m.nodes.length || 1; cy /= m.nodes.length || 1;
  const ox = node.x - cx, oy = node.y - cy;
  m.nodes.forEach(n => { n.x += ox; n.y += oy; });
  node.el && node.el.remove();
  nodes.splice(nodes.indexOf(node), 1);
  edges = edges.filter(e => e.from !== node.id && e.to !== node.id);
  m.nodes.forEach(n => { nodes.push(n); renderNode(n, false); });
  m.edges.forEach(e => edges.push(e));
  redrawEdges(); renderTree(); persist();
  clearSel();
  status("unpacked " + node.label + " (" + m.nodes.length + " nodes)");
}
// specOf variant that resolves labels against a LOCAL array (pack):
function specOfLocal(ns, es) {
  const byId = new Map(ns.map(n => [n.id, n]));
  return JSON.parse(JSON.stringify({ nodes: ns.map(n => ({
    type: n.type, label: n.label, x: Math.round(n.x),
    y: Math.round(n.y), grp: n.grp || undefined,
    sub: n.sub || undefined })), edges: es.map(e => ({
    from: byId.get(e.from).label, to: byId.get(e.to).label,
    sign: e.sign, gain: e.gain, tag: e.tag,
    pts: (e.pts || []).map(p => ({
      x: Math.round(p.x * 10) / 10, y: Math.round(p.y * 10) / 10 })) })),
    groups: {} }));
}

// ---------------- interactions ----------------
svg.addEventListener("click", () => {
  if (suppressClick) { suppressClick = false; return; }
  if (!linkFrom) clearSel();
});
// Stale-click guard: a new press always resets suppressClick (e.g.
// a group drag released outside the canvas never fires a click).
window.addEventListener("mousedown",
  () => { suppressClick = false; }, true);
// Double-click EMPTY canvas = exit a subsystem (or zoom out at top)
svg.addEventListener("dblclick", ev => {
  if (ev.target !== svg) return;
  if (ctxStack.length) { exitSub(); return; }
  selGrp = null;
  fitView();
  renderTree();
  status("zoomed back out to the whole model");
});
// Drag from EMPTY canvas = rubber-band select every node the
// rectangle touches (replaces the current selection).
svg.addEventListener("mousedown", ev => {
  if (ev.target !== svg) return;          // node/edge hit: other paths
  if (ev.shiftKey || ev.altKey) return;   // modifiers stay link tools
  const p = svgPoint(ev);
  band = { x0: p.x, y0: p.y, x1: p.x, y1: p.y, moved: false };
  window.addEventListener("mousemove", onMove);
  window.addEventListener("mouseup", endDrag);
});

function startDrag(ev, id) {
  const nd = nodeById(id);
  if (ev.shiftKey || ev.altKey) {
    // pending until real movement: a no-move SHIFT press stays a
    // click (selection toggle), a no-move ALT press stays a click
    // (arms the two-click link) — neither creates anything alone.
    // v2.1 returned here WITHOUT its own listeners, so the link
    // gesture only completed via a stale linkDrag on the NEXT
    // gesture (and duplicated edges on the two-click flow); this
    // gesture is now self-contained.
    const p = svgPoint(ev);
    linkDrag = { from: id, x0: p.x, y0: p.y, moved: false };
    snapshot();
    status("drawing synapse from " + nd.label +
           " — release on the target");
    window.addEventListener("mousemove", onMove);
    window.addEventListener("mouseup", endDrag);
    return;
  }
  if (selSet.has(id) && selSet.size > 1) {
    // group drag: the grabbed node leads; every other selected
    // node keeps its offset from it. Positions mutate exactly like
    // the single-node drag (n.x/n.y), so localStorage, JSON
    // export and undo all see group moves.
    const p = svgPoint(ev);
    drag = { id, moved: false,
             dx: p.x - nd.x, dy: p.y - nd.y,
             group: [...selSet].map(gid => {
               const m = nodeById(gid);
               return m ? { id: gid, sx: m.x, sy: m.y } : null;
             }).filter(Boolean) };
  } else {
    const p = svgPoint(ev);
    drag = { id, moved: false,
             dx: p.x - nd.x, dy: p.y - nd.y };
  }
  window.addEventListener("mousemove", onMove);
  window.addEventListener("mouseup", endDrag);
}
let tempLine = null;
function onMove(ev) {
  if (bendDrag) {
    const p = svgPoint(ev);
    const cur = bendDrag.e.pts[bendDrag.k];
    if (!bendDrag.moved) {
      if (Math.hypot(p.x - cur.x, p.y - cur.y) <= 2) return;
      snapshot();               // pre-drag state, once per gesture
    }
    bendDrag.moved = true;
    bendDrag.e.pts[bendDrag.k] = { x: p.x, y: p.y };
    redrawEdges();
    return;
  }
  if (linkDrag) {
    if (!linkDrag.moved) {
      const p0 = svgPoint(ev);
      if (Math.hypot(p0.x - linkDrag.x0, p0.y - linkDrag.y0) <= 3)
        return;                      // still within click territory
      linkDrag.moved = true;
    }
    if (!tempLine) {
      tempLine = document.createElementNS(NS, "line");
      tempLine.setAttribute("stroke", "#B2182B");
      tempLine.setAttribute("stroke-width", "2");
      tempLine.setAttribute("stroke-dasharray", "6,4");
      world.appendChild(tempLine);
    }
    const a = nodeById(linkDrag.from);
    const tp = svgPoint(ev);
    tempLine.setAttribute("x1", a.x); tempLine.setAttribute("y1", a.y);
    tempLine.setAttribute("x2", tp.x);
    tempLine.setAttribute("y2", tp.y);
    return;
  }
  if (band) {
    const p = svgPoint(ev);
    band.x1 = p.x; band.y1 = p.y;
    if (Math.hypot(p.x - band.x0, p.y - band.y0) > 3)
      band.moved = true;
    if (!bandRect) {
      bandRect = document.createElementNS(NS, "rect");
      bandRect.setAttribute("fill", "rgba(205,52,181,0.08)");
      bandRect.setAttribute("stroke", "#CD34B5");
      bandRect.setAttribute("stroke-width", "1");
      bandRect.setAttribute("stroke-dasharray", "5,4");
      bandRect.setAttribute("pointer-events", "none");
      world.appendChild(bandRect);
    }
    bandRect.setAttribute("x", Math.min(band.x0, band.x1));
    bandRect.setAttribute("y", Math.min(band.y0, band.y1));
    bandRect.setAttribute("width", Math.abs(band.x1 - band.x0));
    bandRect.setAttribute("height", Math.abs(band.y1 - band.y0));
    return;
  }
  if (!drag) return;
  if (drag.group) {
    const p = svgPoint(ev);
    const nx = p.x - drag.dx, ny = p.y - drag.dy;
    const lead = drag.group.find(m => m.id === drag.id);
    const ddx = nx - lead.sx, ddy = ny - lead.sy;
    drag.group.forEach(m => {
      const n2 = nodeById(m.id);
      if (!n2) return;
      n2.x = m.sx + ddx;
      n2.y = m.sy + ddy;
      n2.el.remove(); renderNode(n2, true);
    });
    drag.moved = true;
    redrawEdges();
    return;
  }
  const nd = nodeById(drag.id);
  const dp = svgPoint(ev);
  nd.x = dp.x - drag.dx;
  nd.y = dp.y - drag.dy;
  drag.moved = true;
  nd.el.remove(); renderNode(nd, true); redrawEdges();
}
function endDrag(ev) {
  window.removeEventListener("mousemove", onMove);
  window.removeEventListener("mouseup", endDrag);
  if (linkDrag) {
    if (!linkDrag.moved && !tempLine) {
      // no movement: this was a click, not a drag. Do nothing here
      // and let the click event decide (SHIFT+click toggles the
      // selection, ALT+click arms the two-click link).
      undoStack.pop();   // no mutation happened
      linkDrag = null;
      return;
    }
    const hp = svgPoint(ev);
    const tgt = nodes.find(n =>
      n.id !== linkDrag.from &&
      Math.hypot(hp.x - n.x, hp.y - n.y) < 24);
    if (tempLine) { tempLine.remove(); tempLine = null; }
    if (tgt) {
      const e = { id: "e" + (uid++), from: linkDrag.from, to: tgt.id,
                  sign: "exc", gain: 0.5, tag: "" };
      edges.push(e);
      redrawEdges();
      selectEdge(e.id);
      status("synapse " + nodeById(linkDrag.from).label + " -> " +
             tgt.label);
    } else {
      linkFrom = linkDrag.from;
      status("release missed a node — click the target to finish");
      undoStack.pop();   // no mutation happened
    }
    linkDrag = null;
    return;
  }
  if (bendDrag) {
    // snapshot() already persisted the pre-drag state; save the final
    // routing. The follow-up click is left alone so the edge ends up
    // selected (panel shows the "remove bend points" button).
    if (bendDrag.moved) persist();
    bendDrag = null;
    return;
  }
  if (band) {
    if (bandRect) { bandRect.remove(); bandRect = null; }
    if (band.moved) {
      suppressClick = true;   // don't let the follow-up click clear it
      const rx0 = Math.min(band.x0, band.x1),
            ry0 = Math.min(band.y0, band.y1),
            rx1 = Math.max(band.x0, band.x1),
            ry1 = Math.max(band.y0, band.y1);
      const hits = nodes.filter(nd => {
        if (isHiddenNode(nd)) return false;   // hidden layers not selectable
        const d = TYPES[nd.type];
        const hw = d.shape === "ellipse" ? 30 :
                   d.shape === "rect" ? 26 : 18;
        const hh = d.shape === "circle" ? 18 : 13;
        return nd.x + hw > rx0 && nd.x - hw < rx1 &&
               nd.y + hh > ry0 && nd.y - hh < ry1;
      });
      if (ev.ctrlKey)          // CTRL+band ADDS to the selection
        hits.forEach(nd => selSet.add(nd.id));
      else
        selSet = new Set(hits.map(nd => nd.id));
      sel = selSet.size === 1 ? { kind: "node",
                                  id: [...selSet][0] } : null;
      applySelClasses();
      updatePanel();
      status(selSet.size ? selSet.size + " node(s) selected" :
             "no nodes in selection rectangle");
    }
    band = null;
    return;
  }
  if (drag && drag.moved) {
    snapshot();
    if (drag.group) suppressClick = true;   // keep group selected
  }
  drag = null;
}
function onClickNode(id, shift, alt) {
  if (alt && !linkFrom) {
    // ALT+click arms the two-click synapse flow (v2.1 used SHIFT
    // for this; SHIFT now toggles group selection)
    linkFrom = id;
    status("linking from " + nodeById(id).label +
           " — ALT-click the target");
    return;
  }
  if (linkFrom && linkFrom !== id) {
    const e = { id: "e" + (uid++), from: linkFrom, to: id,
                sign: "exc", gain: 0.5, tag: "" };
    edges.push(e);
    redrawEdges();
    status("synapse " + nodeById(linkFrom).label + " -> " +
           nodeById(id).label);
    linkFrom = null;
    selectEdge(e.id);
    snapshot();
    return;
  }
  if (shift && !linkFrom) {
    toggleNodeSel(id);
    return;
  }
  selectNode(id);
}

// ---------------- copy / paste (v3.0) ----------------
// Internal clipboard in the serialized (label-based) form, so a copy
// can paste into ANY tab. Paste relabels ("_2", "_3", ... first free
// suffix), keeps internal edges, drops edges that pointed outside the
// copied set, lands centered on the cursor, and selects the new group.
let clip = null;
function uniqueLabel(base) {
  const taken = new Set(nodes.map(n => n.label));
  if (!taken.has(base)) return base;
  const m = base.match(/^(.*?)(\d+)$/);
  let i = 2;
  let stem = base, n = 0;
  if (m && m[1]) { stem = m[1]; n = parseInt(m[2], 10) || 0; }
  let cand;
  do { cand = stem + (n + i); i++; } while (taken.has(cand));
  return cand;
}
function copySel() {
  if (!selSet.size) { status("nothing selected to copy"); return; }
  const ids = new Set(selSet);
  const ns = nodes.filter(n => ids.has(n.id))
    .map(n => ({ type: n.type, label: n.label, x: n.x, y: n.y,
                 grp: n.grp }));
  const es = edges.filter(e => ids.has(e.from) && ids.has(e.to))
    .map(e => ({ from: nodeById(e.from).label, to: nodeById(e.to).label,
                 sign: e.sign, gain: e.gain, tag: e.tag,
                 pts: (e.pts || []).map(p => ({ x: p.x, y: p.y })) }));
  clip = { nodes: ns, edges: es };
  status("copied " + ns.length + " node(s), " + es.length +
         " internal synapse(s) — Ctrl+V to paste");
}
function pasteClip() {
  if (!clip) { status("clipboard empty — select and Ctrl+C first"); return; }
  const map = {};
  let x0 = Infinity, y0 = Infinity;
  clip.nodes.forEach(n => { x0 = Math.min(x0, n.x); y0 = Math.min(y0, n.y); });
  const ox = lastMouseWorld.x - x0, oy = lastMouseWorld.y - y0;
  snapshot();
  const newIds = [];
  clip.nodes.forEach(n => {
    const nl = uniqueLabel(n.label);
    const nd = addNode(n.type, n.x + ox, n.y + oy, nl, true, n.grp,
                       n.sub);
    map[n.label] = nd.id;
    newIds.push(nd.id);
  });
  clip.edges.forEach(se => {
    if (map[se.from] && map[se.to]) {
      edges.push({ id: "e" + (uid++), from: map[se.from], to: map[se.to],
                   sign: se.sign, gain: se.gain, tag: se.tag,
                   pts: (se.pts || []).map(p => ({ x: p.x, y: p.y })) });
    }
  });
  redrawEdges();
  persist();
  sel = null;
  selSet = new Set(newIds);
  applySelClasses();
  updatePanel();
  renderTree();
  status("pasted " + clip.nodes.length + " node(s) at cursor — still on " +
         "the clipboard for another Ctrl+V");
}

// ---------------- delete / keyboard ----------------
function deleteSel() {
  if (sel && sel.kind === "edge") {
    snapshot();
    const e = edgeById(sel.id);
    e.el.remove();
    edges.splice(edges.indexOf(e), 1);
    redrawEdges();     // crossing hops elsewhere may disappear with it
    clearSel();
    return;
  }
  // node selection: the whole group (rubber-band / ctrl+click / paste)
  const targets = selSet.size ? [...selSet] : (sel ? [sel.id] : []);
  if (!targets.length) return;
  snapshot();
  const tset = new Set(targets);
  for (let i = edges.length - 1; i >= 0; i--) {
    if (tset.has(edges[i].from) || tset.has(edges[i].to)) {
      edges[i].el && edges[i].el.remove();
      edges.splice(i, 1);
    }
  }
  for (let i = nodes.length - 1; i >= 0; i--) {
    if (tset.has(nodes[i].id)) {
      nodes[i].el && nodes[i].el.remove();
      nodes.splice(i, 1);
    }
  }
  redrawEdges();
  renderTree();
  clearSel();
}
document.getElementById("bDel").onclick = deleteSel;
document.getElementById("bDel2").onclick = deleteSel;
window.addEventListener("keydown", ev => {
  const tag = (ev.target.tagName || "").toLowerCase();
  if (tag === "input" || tag === "select" || tag === "textarea") return;
  if (ev.key === "F2") {
    ev.preventDefault();
    if (!selGrp) {
      status("click a layer row (left tree) first, then F2 hides it");
      return;
    }
    toggleGrpHidden(selGrp);
  } else if (ev.key === "F3") {
    // show everything hidden again (the tree lists hidden layers with
    // a struck-through row, so you can see what F3 will bring back)
    ev.preventDefault();
    if (!hiddenG.size) { status("nothing hidden"); return; }
    hiddenG.clear();
    loadIntoCanvas(serialize(), true);
    renderTree();
    status("all hidden layers shown again");
  } else if (ev.key === "F4") {
    // back out of a drilled layer OR exit a subsystem
    ev.preventDefault();
    if (ctxStack.length) { exitSub(); return; }
    selGrp = null;
    fitView();
    renderTree();
    status("backed out to the whole model");
  } else if ((ev.ctrlKey || ev.metaKey) && ev.shiftKey &&
             ev.key.toLowerCase() === "g") {
    ev.preventDefault();
    unpackSub();
  } else if ((ev.ctrlKey || ev.metaKey) &&
             ev.key.toLowerCase() === "g") {
    ev.preventDefault();
    packSub();
  } else if (ev.key === "Delete" || ev.key === "Backspace") {
    ev.preventDefault();
    deleteSel();
  } else if ((ev.ctrlKey || ev.metaKey) &&
             ev.key.toLowerCase() === "c") {
    ev.preventDefault();
    copySel();
  } else if ((ev.ctrlKey || ev.metaKey) &&
             ev.key.toLowerCase() === "v") {
    ev.preventDefault();
    pasteClip();
  } else if (ev.key === "Escape") {
    // cancel any rubber band / armed link and clear the selection;
    // inside a subsystem Esc backs out of it too
    if (bandRect) { bandRect.remove(); bandRect = null; }
    band = null;
    linkFrom = null;
    if (ctxStack.length) { exitSub(); return; }
    clearSel();
    status("selection cleared (Esc)");
  } else if ((ev.ctrlKey || ev.metaKey) && !ev.shiftKey &&
             ev.key.toLowerCase() === "z") {
    ev.preventDefault(); undo();
  } else if ((ev.ctrlKey || ev.metaKey) &&
             (ev.key.toLowerCase() === "y" ||
              (ev.shiftKey && ev.key.toLowerCase() === "z"))) {
    ev.preventDefault(); redo();
  }
});
document.getElementById("bUndo").onclick = undo;
document.getElementById("bRedo").onclick = redo;

// ---------------- save / load ----------------
document.getElementById("bSave").onclick = () => {
  const data = modelSpec();
  const blob = new Blob([JSON.stringify(data, null, 1)],
                        { type: "application/json" });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = models[activeTab].name.replace(/\W+/g, "_") +
               "_connectome.json";
  a.click();
};
document.getElementById("bLoad").onclick = () =>
  document.getElementById("fileIn").click();
document.getElementById("fileIn").addEventListener("change", ev => {
  const f = ev.target.files[0];
  if (!f) return;
  f.text().then(txt => {
    snapshot();
    loadSpec(JSON.parse(txt));
  });
});
document.getElementById("bClear").onclick = () => {
  while (ctxStack.length) exitSub();
  loadSpec({ nodes: [], edges: [] });
};
function loadSpec(spec, keepUndo) {
  // TAB-LEVEL load: leaving a subsystem context first (the stack's
  // edits are written back by exitSub), then reload the top canvas.
  while (ctxStack.length) exitSub();
  reloadCurrent(spec, keepUndo);
}
function loadIntoCanvas(data, keepUndo) {
  // CONTEXT-AWARE reload: replaces the CURRENT context's content
  // (undo/redo/layer-hide inside a subsystem stay inside it).
  reloadCurrent(data, keepUndo);
}
function reloadCurrent(spec, keepUndo) {
  if (!keepUndo) snapshot();
  clearCanvasDom();   // ALL node/edge elements, incl. any orphans
  if (bandRect) { bandRect.remove(); bandRect = null; }
  band = null;
  nodes = []; edges = []; sel = null; selSet.clear(); linkFrom = null;
  hiddenG = new Set(spec.hidden || []);
  grpNames = spec.groups || {};
  selGrp = null;
  const byLabel = {};
  try {
    (spec.nodes || []).forEach(nd => {
      if (!TYPES[nd.type]) {
        status("UNKNOWN TYPE: " + nd.type + " (" + nd.label + ")");
        return;
      }
      const n = addNode(nd.type, nd.x, nd.y, nd.label, true, nd.grp,
                         nd.sub);
      byLabel[nd.label] = n.id;
    });
    // accept both key names (templates use `edges`, files use `synapses`)
    const syn = spec.synapses || spec.edges || [];
    syn.forEach(sn => {
      if (byLabel[sn.from] && byLabel[sn.to]) {
        edges.push({ id: "e" + (uid++), from: byLabel[sn.from],
                     to: byLabel[sn.to], sign: sn.sign,
                     gain: sn.gain, tag: sn.tag || "",
                     pts: (sn.pts || []).map(p => Array.isArray(p)
                       ? { x: +p[0], y: +p[1] } : { x: +p.x, y: +p.y })
                       .filter(p => isFinite(p.x) && isFinite(p.y)) });
      } else {
        status("edge skipped: " + sn.from + " -> " + sn.to);
      }
    });
  } catch (err) {
    status("ERROR loading: " + err.message);
  }
  redrawEdges();
  renderTree();
  persist();
}

// ---------------- TEMPLATES (full models, literature gains) --------
// One complete side of the literature model (circuit_dengstyle /
// circuit_literature equivalent): rhythm + lamination + PF +
// cross-reciprocal PF-INs + full reflex layer with INs + Renshaw +
// heel/toe + KINH + stance-Ib. Everything present, nothing direct.
function litSide(S) {
  const N = [], E = [];
  const nd = (t, l, x, y) => { N.push({ type: t, label: l, x, y });
                              return l; };
  const e = (f, t, sg, g, tag) =>
    E.push({ from: f, to: t, sign: sg, gain: g, tag: tag || S });
  const X = 400;
  // descending inputs
  nd("PORT-load", "DRIVE", 70, 150);
  nd("PORT-load", "POSTURE", 70, 240);
  // heel/toe cutaneous SN -> IN (R5)
  nd("SN-heel", "heel_SN_" + S, 170, 40);
  nd("SN-toe", "toe_SN_" + S, 170, 95);
  nd("IN-C", "heel_IN_" + S, 265, 60);
  // rhythm layer + lamination
  nd("HC-RG-E", "RG_E_" + S, 350, 170);
  nd("HC-RG-F", "RG_F_" + S, 350, 270);
  nd("IN-InE", "InE_" + S, 265, 190);
  nd("IN-InF", "InF_" + S, 265, 250);
  e("DRIVE", "RG_E_" + S, "exc", 1.0);
  e("DRIVE", "RG_F_" + S, "exc", 0.8);
  e("POSTURE", "RG_E_" + S, "exc", 0.6);
  e("RG_E_" + S, "InE_" + S, "exc", 1.0);
  e("InE_" + S, "RG_F_" + S, "inh", 1.0);
  e("RG_F_" + S, "InF_" + S, "exc", 1.0);
  e("InF_" + S, "RG_E_" + S, "inh", 1.0);
  e("heel_IN_" + S, "InE_" + S, "exc", 0.5, "heel_reset");
  e("heel_IN_" + S, "InF_" + S, "inh", 0.5, "heel_release");
  e("heel_SN_" + S, "heel_IN_" + S, "exc", 1.0);
  e("toe_SN_" + S, "heel_IN_" + S, "exc", 0.6, "toe_prolong");
  // pattern formation + cross-reciprocal PF-INs
  const pf = {};
  [["E1", 480, 130], ["E2", 480, 180], ["F1", 480, 240],
   ["F2", 480, 290]].forEach(([w, x, y]) =>
    pf[w] = nd("HC-PF-" + w[0], "PF_" + w + "_" + S, x, y));
  nd("IN-PF", "PF_IN_E_" + S, 575, 155);
  nd("IN-PF", "PF_IN_F_" + S, 575, 265);
  e("RG_E_" + S, "PF_E1_" + S, "exc", 1.0);
  e("RG_E_" + S, "PF_E2_" + S, "exc", 0.9);
  e("RG_F_" + S, "PF_F1_" + S, "exc", 1.0);
  e("RG_F_" + S, "PF_F2_" + S, "exc", 0.8);
  e("PF_E1_" + S, "PF_IN_E_" + S, "exc", 1.0);
  e("PF_E2_" + S, "PF_IN_E_" + S, "exc", 0.8);
  e("PF_F1_" + S, "PF_IN_F_" + S, "exc", 1.0);
  e("PF_F2_" + S, "PF_IN_F_" + S, "exc", 0.7);
  e("PF_IN_E_" + S, "PF_F1_" + S, "inh", 0.8);
  e("PF_IN_E_" + S, "PF_F2_" + S, "inh", 0.6);
  e("PF_IN_F_" + S, "PF_E1_" + S, "inh", 0.8);
  e("PF_IN_F_" + S, "PF_E2_" + S, "inh", 0.6);
  // motoneuron pools + muscles (extensor / flexor antagonist pair)
  nd("MN", "MN_ext_" + S, 720, 150);
  nd("MN", "MN_flex_" + S, 720, 260);
  nd("MUSCLE", "extensor_" + S, 850, 150);
  nd("MUSCLE", "flexor_" + S, 850, 260);
  e("PF_E1_" + S, "MN_ext_" + S, "exc", 1.0);
  e("PF_E2_" + S, "MN_ext_" + S, "exc", 0.9);
  e("PF_F1_" + S, "MN_flex_" + S, "exc", 1.0);
  e("PF_F2_" + S, "MN_flex_" + S, "exc", 0.7);
  e("MN_ext_" + S, "extensor_" + S, "exc", 1.0);
  e("MN_flex_" + S, "flexor_" + S, "exc", 1.0);
  // proprioceptive SNs
  nd("SN-Ia", "Ia_ext_" + S, 170, 350);
  nd("SN-Ia", "Ia_flex_" + S, 170, 470);
  nd("SN-II", "II_ext_" + S, 170, 395);
  nd("SN-II", "II_flex_" + S, 170, 515);
  nd("SN-Ib", "Ib_ext_" + S, 170, 440);
  nd("SN-Ib", "Ib_flex_" + S, 170, 560);
  // Ia: mono + reciprocal + IaIN mutual (rules 1)
  e("Ia_ext_" + S, "MN_ext_" + S, "exc", 1.0, "ia_mono");
  nd("IN-IaIN", "IaIN_ext_" + S, 640, 340);
  nd("IN-IaIN", "IaIN_flex_" + S, 640, 470);
  e("Ia_ext_" + S, "IaIN_ext_" + S, "exc", 1.0);
  e("Ia_flex_" + S, "IaIN_flex_" + S, "exc", 1.0);
  e("IaIN_ext_" + S, "MN_flex_" + S, "inh", 0.5, "ia_recip");
  e("IaIN_flex_" + S, "MN_ext_" + S, "inh", 0.5, "ia_recip");
  e("IaIN_ext_" + S, "IaIN_flex_" + S, "inh", 0.5, "iaIN_mutual");
  e("IaIN_flex_" + S, "IaIN_ext_" + S, "inh", 0.5, "iaIN_mutual");
  // II: two collateral INs (rule 2)
  nd("IN-IIe", "IIe_ext_" + S, 640, 395);
  nd("IN-IIi", "IIi_flex_" + S, 640, 515);
  e("II_ext_" + S, "IIe_ext_" + S, "exc", 1.0);
  e("IIe_ext_" + S, "MN_ext_" + S, "exc", 0.5, "ii_ag");
  e("II_flex_" + S, "IIi_flex_" + S, "exc", 1.0);
  e("IIi_flex_" + S, "MN_ext_" + S, "inh", 0.5, "ii_ant");
  // Ib: autogenic via IbIN + mutual (rule 3); reversal via IBEXC (4)
  nd("IN-IbIN", "IbIN_ext_" + S, 640, 440);
  nd("IN-IbIN", "IbIN_flex_" + S, 640, 560);
  e("Ib_ext_" + S, "IbIN_ext_" + S, "exc", 1.0);
  e("IbIN_ext_" + S, "MN_ext_" + S, "inh", 0.5, "ib_auto");
  e("IbIN_ext_" + S, "IbIN_flex_" + S, "inh", 0.5, "ibIN_mutual");
  e("IbIN_flex_" + S, "IbIN_ext_" + S, "inh", 0.5, "ibIN_mutual");
  nd("IN-Ib+", "IBEXC_" + S, 700, 200);
  e("Ib_ext_" + S, "IBEXC_" + S, "exc", 1.0);
  e("RG_E_" + S, "IBEXC_" + S, "exc", 1.0, "stance_gate");
  e("IBEXC_" + S, "MN_ext_" + S, "exc", 0.5, "ib_reversal");
  nd("IN-LBIN", "LBIN_" + S, 700, 245);
  e("IBEXC_" + S, "LBIN_" + S, "exc", 0.5);
  e("LBIN_" + S, "RG_E_" + S, "exc", 0.5, "stance_prolong");
  // Renshaw set (rule 5) + recurrent disinhibition
  nd("RC", "RC_ext_" + S, 800, 40);
  nd("RC", "RC_flex_" + S, 800, 330);
  e("MN_ext_" + S, "RC_ext_" + S, "exc", 1.0);
  e("RC_ext_" + S, "MN_ext_" + S, "inh", 0.5, "rc_own");
  e("MN_flex_" + S, "RC_flex_" + S, "exc", 1.0);
  e("RC_flex_" + S, "MN_flex_" + S, "inh", 0.5, "rc_own");
  e("RC_ext_" + S, "RC_flex_" + S, "inh", 0.5, "rc_mutual");
  e("RC_flex_" + S, "RC_ext_" + S, "inh", 0.5, "rc_mutual");
  e("RC_ext_" + S, "IaIN_ext_" + S, "inh", 0.5, "rc_disinhib");
  e("RC_flex_" + S, "IaIN_flex_" + S, "inh", 0.5, "rc_disinhib");
  // KINH swing suppression
  nd("IN-KINH", "KINH_" + S, 660, 90);
  e("PF_F1_" + S, "KINH_" + S, "exc", 1.0);
  e("KINH_" + S, "MN_ext_" + S, "inh", 0.5, "swing_suppress");
  return { nodes: N, edges: E };
}

// bilateral assembly: two mirrored sides + c1/V3 commissurals
function buildBilateral() {
  const R = litSide("R");
  const L = litSide("L");
  L.nodes.forEach(n => { n.x = 940 - n.x; });   // mirror layout
  // labels already carry the correct side suffix from litSide(S)
  const N = R.nodes.concat(L.nodes), E = R.edges.concat(L.edges);
  const c1r = { type: "IN-C", label: "c1_R", x: 520, y: 360 };
  const c1l = { type: "IN-C", label: "c1_L", x: 420, y: 360 };
  const v3r = { type: "IN-V3", label: "V3_R", x: 520, y: 30 };
  const v3l = { type: "IN-V3", label: "V3_L", x: 420, y: 30 };
  N.push(c1r, c1l, v3r, v3l);
  const x = (f, t, sg, g, tag) => E.push({ from: f, to: t, sign: sg,
    gain: g, tag: tag || "commissural" });
  x("RG_F_R", "c1_R", "exc", 1.0);
  x("c1_R", "RG_F_L", "inh", 0.1);
  x("RG_F_L", "c1_L", "exc", 1.0);
  x("c1_L", "RG_F_R", "inh", 0.1);
  x("RG_E_R", "V3_R", "exc", 1.0);
  x("V3_R", "RG_E_L", "exc", 0.1);
  x("RG_E_L", "V3_L", "exc", 1.0);
  x("V3_L", "RG_E_R", "exc", 0.1);
  return { nodes: N, edges: E };
}

// W2L contact-driven equivalent (from the documented working wiring)
function buildW2L() {
  return { nodes: [
    { type: "SN-heel", label: "ContactFoot_R", x: 100, y: 120 },
    { type: "SN-heel", label: "ContactFoot_L", x: 100, y: 260 },
    { type: "HC-RG-E", label: "RG_E_R", x: 300, y: 90 },
    { type: "HC-RG-F", label: "RG_F_R", x: 300, y: 160 },
    { type: "HC-RG-E", label: "RG_E_L", x: 300, y: 230 },
    { type: "HC-RG-F", label: "RG_F_L", x: 300, y: 300 },
    { type: "HC-PF-E", label: "PF_E_R", x: 480, y: 90 },
    { type: "HC-PF-F", label: "PF_F_R", x: 480, y: 160 },
    { type: "HC-PF-E", label: "PF_E_L", x: 480, y: 230 },
    { type: "HC-PF-F", label: "PF_F_L", x: 480, y: 300 },
    { type: "MN", label: "ext MN_R", x: 640, y: 100 },
    { type: "MN", label: "flx MN_R", x: 640, y: 170 },
    { type: "MN", label: "ext MN_L", x: 640, y: 240 },
    { type: "MN", label: "flx MN_L", x: 640, y: 310 },
    { type: "MUSCLE", label: "extensor_R", x: 780, y: 100 },
    { type: "MUSCLE", label: "flexor_R", x: 780, y: 170 },
    { type: "MUSCLE", label: "extensor_L", x: 780, y: 240 },
    { type: "MUSCLE", label: "flexor_L", x: 780, y: 310 } ],
    synapses: [
      { from: "ContactFoot_R", to: "RG_E_R", sign: "exc", gain: 6.0,
        tag: "w2l_contact" },
      { from: "ContactFoot_R", to: "RG_F_L", sign: "exc", gain: 6.0,
        tag: "w2l_contact_contra" },
      { from: "ContactFoot_L", to: "RG_E_L", sign: "exc", gain: 6.0,
        tag: "w2l_contact" },
      { from: "ContactFoot_L", to: "RG_F_R", sign: "exc", gain: 6.0,
        tag: "w2l_contact_contra" },
      { from: "RG_E_R", to: "PF_E_R", sign: "exc", gain: 0.5,
        tag: "w2l" },
      { from: "RG_F_R", to: "PF_F_R", sign: "exc", gain: 0.5,
        tag: "w2l" },
      { from: "RG_E_L", to: "PF_E_L", sign: "exc", gain: 0.5,
        tag: "w2l" },
      { from: "RG_F_L", to: "PF_F_L", sign: "exc", gain: 0.5,
        tag: "w2l" },
      { from: "PF_E_R", to: "ext MN_R", sign: "exc", gain: 0.5,
        tag: "w2l" },
      { from: "PF_F_R", to: "flx MN_R", sign: "exc", gain: 0.5,
        tag: "w2l" },
      { from: "PF_E_L", to: "ext MN_L", sign: "exc", gain: 0.5,
        tag: "w2l" },
      { from: "PF_F_L", to: "flx MN_L", sign: "exc", gain: 0.5,
        tag: "w2l" },
      { from: "ext MN_R", to: "extensor_R", sign: "exc", gain: 1.0,
        tag: "w2l" },
      { from: "flx MN_R", to: "flexor_R", sign: "exc", gain: 1.0,
        tag: "w2l" },
      { from: "ext MN_L", to: "extensor_L", sign: "exc", gain: 1.0,
        tag: "w2l" },
      { from: "flx MN_L", to: "flexor_L", sign: "exc", gain: 1.0,
        tag: "w2l" } ] };
}

// Deng A6 full pair: the reflex-layer template (agonist + antagonist,
// every IN present, mutual inhibitions, Renshaw)
function buildDengPair() {
  const s = litSide("R");
  // trim to the reflex core: keep SN/IN/MN/RC/IBEXC nodes only
  const keep = /^(Ia_|II_|Ib_|IaIN_|IbIN_|IIe_|IIi_|IBEXC_|RC_|MN_|extensor_|flexor_)/;
  return { nodes: s.nodes.filter(n => keep.test(n.label)),
           edges: s.edges.filter(e =>
             keep.test(e.from) && keep.test(e.to)) };
}

// layer stamping for the BUILT-IN templates — MIRROR of
// stamp_walker_grps() in make_editor_templates.py; keep rules identical
function stampWalkerGrps(spec) {
  const G = { drive: "Drive / posture inputs" };
  ["R", "L"].forEach(S => {
    G["rg_" + S] = "Rhythm RG + lamination (" + S + ")";
    G["pf_" + S] = "Pattern formation — 4 cells (" + S + ")";
    G["mn_" + S] = "Motoneurons, per muscle (" + S + ")";
    G["mus_" + S] = "Muscles (" + S + ")";
    G["aff_" + S] = "Afferents (" + S + ")";
    G["in_" + S] = "Reflex / phase INs (" + S + ")";
  });
  spec.nodes.forEach(n => {
    let l = n.label;
    let gid = null;
    if (l === "DRIVE" || l === "POSTURE" || l.startsWith("pf_gain=") ||
        l.startsWith("PM ")) gid = "drive";
    else {
      if (l.endsWith(" (pruned)")) l = l.slice(0, -9);  // pruned MNs
      const S = l.endsWith("_R") ? "R" : l.endsWith("_L") ? "L" : null;
      if (S) {
        if (l.startsWith("RG_") || l === "InE_" + S ||
            l === "InF_" + S || l.startsWith("c1_") ||
            l.startsWith("V3_") || l.startsWith("heel_IN_"))
          gid = "rg_" + S;
        else if (l.startsWith("PF_")) gid = "pf_" + S;
        else if (l.startsWith("MN_") || l.includes(" MN_"))
          gid = "mn_" + S;
        else if (/^(KINH_|PRESET_|IaIN_|LBIN_|IBEXC_|RC_|IIe_|IIi_|IbIN_)/
                 .test(l)) gid = "in_" + S;
        else if (/^(Ia_|II_|Ib_|heel_|toe_|ContactFoot_)/.test(l) ||
                 l.endsWith("_sig_" + S)) gid = "aff_" + S;
        else gid = "mus_" + S;
      }
    }
    if (gid) n.grp = gid;
  });
  spec.groups = G;
  return spec;
}

const TEMPLATES = {
  deng: { grp: "builtin", name: "Deng A6 full pair (ag + ant)",
          build: () => stampWalkerGrps(buildDengPair()) },
  lit:  { grp: "builtin", name: "Literature full model (one side)",
          build: () => stampWalkerGrps(litSide("R")) },
  bilateral: { grp: "builtin", name: "BILATERAL literature model",
               build: () => stampWalkerGrps(buildBilateral()) },
  w2l:  { grp: "builtin", name: "W2L contact-driven RG (coarse)",
          build: () => stampWalkerGrps(buildW2L()) },
  blank:{ grp: "builtin", name: "Blank", build: () => ({ nodes: [], edges: [] }) },
};
// w2l_equivalent_draft.json: loaded from disk via a synchronous
// XMLHttpRequest (works on file:// in most browsers)
TEMPLATES["w2lfile"] = {
  grp: "builtin", name: "W2L from w2l_equivalent_draft.json",
  build: () => {
    const xhr = new XMLHttpRequest();
    xhr.open("GET", "w2l_equivalent_draft.json", false);
    xhr.send();
    if (xhr.status === 200 || xhr.status === 0) {
      return stampWalkerGrps(JSON.parse(xhr.responseText));
    }
    throw new Error("cannot load w2l_equivalent_draft.json (status " +
                    xhr.status + ")");
  }
};

// ---- sidecar template library (connectome_templates.json, generated by
// make_editor_templates.py: AnimatLab mines, SNS walker versions, synergy,
// rules, replication connectomes). Loaded once at boot, sync, file://-safe.
const TPL_GROUPS = {
  animatlab: { el: "ogAnimat", label: "AnimatLab models (mined 2026-09-23)" },
  walker:    { el: "ogWalker", label: "SNS walker versions (winners)" },
  lit:       { el: "ogLit",    label: "Synergy / rules / literature" },
};
const TPL_JSON_GROUP = {
  biped2xcpg: "animatlab", w2laproj: "animatlab", bilateralrg: "animatlab",
  li: "animatlab",
  walker_v9: "walker", walker_v10: "walker", walker_s1: "walker",
  walker_s2: "walker", walker_s3: "walker", walker_s3k: "walker",
  w2lvar: "walker", syn6: "walker",
  walker_s3k_flat: "walker", w2lvar_flat: "walker", syn6_flat: "walker",
  walker_s3k_pruned: "walker", walker_s3k_pruned_flat: "walker",
  synergy6: "lit", rules: "lit", rules_motifs: "lit",
  shevtsova: "lit", shinohara: "lit", rybak: "lit",
  spiking_mirror: "walker",
};
const TPL_JSON_NAMES = {
  biped2xcpg: "Biped_2xCPG_wSubs (connectivity audit)",
  w2laproj: "Walker_2_Layer_CPG (aproj mine)",
  bilateralrg: "Walker_2_Layer_CPG_BilateralRG",
  li: "Li Model (walk tester rearranged)",
  walker_v9: "walker v9 winner", walker_v10: "walker v10 winner",
  walker_s1: "curriculum s1 (air, deaff)",
  walker_s2: "curriculum s2 (air + aff)",
  walker_s3: "curriculum s3 (ground)",
  walker_s3k: "s3k — CURRENT production (schematic)",
  w2lvar: "w2lvar — W2L-layout s3k variant (schematic)",
  syn6: "syn6 — 6-synergy walker (schematic)",
  walker_s3k_flat: "s3k FLAT (pre-schematic wiring)",
  walker_s3k_pruned: "s3k PRUNED — 4-way combo (09-30, schematic)",
  walker_s3k_pruned_flat: "s3k PRUNED FLAT (interleg/contact/Ib/rgweak cut)",
  w2lvar_flat: "w2lvar FLAT (pre-schematic wiring)",
  syn6_flat: "syn6 FLAT (pre-schematic wiring)",
  synergy6: "6 synergies per leg (rank-6 NMF)",
  rules: "Circuit rules — Ben (2026-09-24)",
  rules_motifs: "Rule motifs (auto-generated, old)",
  shevtsova: "Shevtsova 2026 laminar RG",
  shinohara: "Shinohara 2025 interlimb/load",
  rybak: "Rybak 2006a two-level RG/PF",
  spiking_mirror: "spiking mirror - hybrid LIF twin of the spinal net (2026-10-02)",
};
let TPL_SRC = "none";
(function loadTplJson() {
  // 1) embedded copy (works on file:// double-click — browsers block
  //    file XHR); 2) sidecar connectome_templates.json as fallback for
  //    regenerated libraries that predate an HTML re-embed.
  let data = window.TPL_DATA_EMBED || null;
  let src = "embedded";
  if (!data) {
    try {
      const xhr = new XMLHttpRequest();
      xhr.open("GET", "connectome_templates.json", false);
      xhr.send();
      if (xhr.status !== 200 && xhr.status !== 0) throw new Error(xhr.status);
      data = JSON.parse(xhr.responseText);
      src = "sidecar";
    } catch (e) {
      status("template library unavailable (no embedded copy, sidecar " +
             "missing): " + e.message);
      return;
    }
  }
  Object.keys(data).forEach(k => {
    const spec = data[k];
    TEMPLATES[k] = { grp: TPL_JSON_GROUP[k] || "lit",
                     name: TPL_JSON_NAMES[k] || k, build: () => spec };
  });
  TPL_SRC = src;
})();
(function renderTplOptions() {
  Object.entries(TPL_GROUPS).forEach(([gid, G]) => {
    const og = document.getElementById(G.el);
    if (!og) return;
    og.label = G.label;
    Object.entries(TEMPLATES).forEach(([k, T]) => {
      if ((TPL_JSON_GROUP[k] || (T.grp === "builtin" ? "builtin" : T.grp))
          !== gid) return;
      const o = document.createElement("option");
      o.value = k;
      o.textContent = T.name || k;
      og.appendChild(o);
    });
  });
})();

// Templates load into a NEW tab (never clobbers current work) and the
// view re-fits so the whole model is visible immediately.
document.getElementById("tplSel").addEventListener("change", ev => {
  const k = ev.target.value;
  ev.target.value = "";
  if (!k || !TEMPLATES[k]) { status("no template selected"); return; }
  let spec;
  try { spec = TEMPLATES[k].build(); }
  catch (err) { status("template ERROR: " + err.message); return; }
  const name = (TEMPLATES[k].name || k).replace(/[^\w ()+.-]/g, "").slice(0, 40);
  while (ctxStack.length) exitSub();
  models[activeTab].data = serialize();   // park current tab (saves it)
  models.push({ name: name, data: { nodes: [], edges: [] } });
  activeTab = models.length - 1;
  undoStack.length = 0; redoStack.length = 0;   // history is per-tab
  loadIntoCanvas(spec, true);
  persist();
  renderTabs();
  fitView();
  status("template loaded into new tab: " + name + " (" +
         nodes.length + " nodes, " + edges.length + " synapses)" +
         (spec && spec._note ? " — " + spec._note :
          " — (no provenance note)"));
});


// ---------------- boot ----------------
renderTabs();
loadIntoCanvas(cur(), true);
status("ready — v3.3: SUBSYSTEMS (double-click a ⊞ node = ENTER it, " +
       "edit its parts, F4 = back out; Ctrl+G = pack selection, " +
       "Ctrl+Shift+G = unpack); LAYERS tree (F2 = hide, F3 = show all); " +
       "ALT+click a LINE = bend point; crossing lines hop; templates " +
       "load into a NEW tab; tabs auto-saved");
