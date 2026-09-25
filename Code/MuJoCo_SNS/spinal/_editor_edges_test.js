// Behavioral test of the v3.1 edge machinery (wire hops, hit lines,
// bend points) against a minimal DOM stub. Extracts the edge section
// (const HOP_R .. applySelClasses) plus serialize() from the HTML and
// exercises the geometry + interaction logic under node.
const fs = require("fs");
const src = fs.readFileSync("connectome_block_editor.html", "utf8");
const js = src.slice(src.lastIndexOf("<script>") + 8,
                     src.lastIndexOf("</script>"));

// browser shims BEFORE eval (renderEdge touches document/world)
function FakeEl() {
  this.children = []; this.dataset = {}; this.textContent = "";
  this.setAttribute = () => {};
  this.appendChild = c => this.children.push(c);
  this.addEventListener = () => {};
  this.classList = { add() {}, toggle() {} };
  this.remove = () => {};
}
const document = { createElementNS: () => new FakeEl() };
const world = { appendChild: () => {} };
const NS = "svg";
const status = () => {};
const selectEdge = () => {};
const snapshot = () => {};
// v3.2 layers live outside the extracted region; tests run unhidden
const isHiddenNode = () => false;
const hiddenG = new Set();
const grpNames = {};

// graph scaffolding (same semantics as the editor's node helpers)
const nodeRadius = () => 18;
const anchor = (nd, tx, ty) => {
  const dx = tx - nd.x, dy = ty - nd.y;
  const L = Math.hypot(dx, dy) || 1;
  const r = nodeRadius() + 2;
  return { x: nd.x + dx / L * r, y: nd.y + dy / L * r };
};
let nodes = [], edges = [];
const nodeById = id => nodes.find(n => n.id === id);

// ---- extract the edge section + serialize ----
const region = js.slice(js.indexOf("const HOP_R"),
                        js.indexOf("function applySelClasses"));
const serFn = js.slice(js.indexOf("function serialize"),
                       js.indexOf("function snapshot"));
eval(region);
eval(serFn);

let pass = 0, fail = 0;
function T(name, cond) {
  if (cond) { pass++; console.log("PASS " + name); }
  else { fail++; console.log("FAIL " + name); }
}

// ---- segXint ----
T("segXint proper crossing",
  (() => { const p = segXint({ x1: 0, y1: 0, x2: 10, y2: 0 },
                             { x1: 5, y1: -5, x2: 5, y2: 5 });
           return p && Math.abs(p.x - 5) < 1e-9 &&
                  Math.abs(p.y) < 1e-9; })());
T("segXint parallel -> null",
  segXint({ x1: 0, y1: 0, x2: 10, y2: 0 },
          { x1: 0, y1: 5, x2: 10, y2: 5 }) === null);
T("segXint collinear overlap -> null",
  segXint({ x1: 0, y1: 0, x2: 10, y2: 0 },
          { x1: 5, y1: 0, x2: 15, y2: 0 }) === null);
T("segXint endpoint touch -> null",
  segXint({ x1: 0, y1: 0, x2: 10, y2: 0 },
          { x1: 10, y1: -5, x2: 10, y2: 5 }) === null);

// ---- edgeVerts ----
nodes = [{ id: "n1", x: 0, y: 0 }, { id: "n2", x: 100, y: 0 }];
edges = [{ id: "e1", from: "n1", to: "n2" }];
let G = edgeVerts(edges[0]);
T("edgeVerts straight: anchor at node rim",
  G.vs.length === 2 && Math.abs(G.vs[0].x - 20) < 1e-9 &&
  Math.abs(G.vs[1].x - 80) < 1e-9 && G.segs.length === 1);
edges[0].pts = [{ x: 50, y: 40 }];
G = edgeVerts(edges[0]);
T("edgeVerts with bend: 3 verts, 2 segs, bend in the middle",
  G.vs.length === 3 && G.segs.length === 2 &&
  Math.abs(G.vs[1].x - 50) < 1e-9 && Math.abs(G.vs[1].y - 40) < 1e-9);

// ---- pointAtFrac ----
let pf = pointAtFrac([{ x1: 20, y1: 0, x2: 80, y2: 0 }], 0.55);
T("pointAtFrac straight 55%",
  Math.abs(pf.x - 53) < 1e-9 && Math.abs(pf.y) < 1e-9);
pf = pointAtFrac([{ x1: 0, y1: 0, x2: 100, y2: 0 },
                  { x1: 100, y1: 0, x2: 100, y2: 100 }], 0.55);
T("pointAtFrac polyline 55% lands 10% into seg 2",
  Math.abs(pf.x - 100) < 1e-9 && Math.abs(pf.y - 10) < 1e-9);

// ---- computeHops: ownership, culling ----
nodes = [{ id: "n1", x: 0, y: 0 }, { id: "n2", x: 200, y: 0 },
         { id: "n3", x: 100, y: -100 }, { id: "n4", x: 100, y: 100 }];
edges = [{ id: "e1", from: "n1", to: "n2" },
         { id: "e2", from: "n3", to: "n4" }];
let hops = computeHops();
T("crossing: later edge (e2) is the hopper at t=0.5",
  hops.get("e2") && hops.get("e2")[0] &&
  Math.abs(hops.get("e2")[0][0].t - 0.5) < 1e-6);
T("crossing: earlier edge (e1) does NOT hop", !hops.get("e1"));
edges = [{ id: "e1", from: "n1", to: "n2" },
         { id: "e2", from: "n1", to: "n3" }];
hops = computeHops();
T("shared endpoint (fan-out): no hops",
  !hops.get("e1") && !hops.get("e2"));
edges = [{ id: "e1", from: "n1", to: "n2" },
         { id: "e2", from: "n1", to: "n2", pts: [{ x: 100, y: 150 }] }];
hops = computeHops();
T("AABB-disjoint bent edge: no hops",
  !hops.get("e1") && !hops.get("e2"));

// ---- addBend: creation + anti-spam cooldown ----
nodes = [{ id: "n1", x: 0, y: 0, label: "A" },
         { id: "n2", x: 200, y: 0, label: "B" }];
edges = [{ id: "e1", from: "n1", to: "n2" }];
addBend(edges[0], 0);
T("addBend: one midpoint at the clicked segment's middle",
  edges[0].pts.length === 1 &&
  Math.abs(edges[0].pts[0].x - 100) < 1e-9 &&
  Math.abs(edges[0].pts[0].y) < 1e-9);
addBend(edges[0], 0);          // held ALT / rapid second click
T("addBend: cooldown blocks rapid re-creation",
  edges[0].pts.length === 1);
setTimeout(() => {
  addBend(edges[0], 1);        // fresh ALT+click on the new segment
  T("addBend after cooldown on seg 2: second bend lands mid-seg-2",
    edges[0].pts.length === 2 &&
    Math.abs(edges[0].pts[1].x - 140) < 1e-9);   // seg2 = 100..180
  // renderEdge runs for real against the stub DOM (hops + handles)
  let ok = true;
  try { redrawEdges(); } catch (e) { ok = false; console.log(e); }
  T("redrawEdges with bends + stub DOM: no throw", ok);

  // ---- serialize keeps pts (undo / tabs / export) ----
  const s = serialize();
  T("serialize: edge carries pts for undo/export round-trip",
    s.edges[0].pts.length === 2 &&
    Math.abs(s.edges[0].pts[1].x - 140) < 1e-9 &&
    s.edges[0].from === "A");

  console.log(`\n${pass} pass, ${fail} fail`);
  process.exit(fail ? 1 : 0);
}, 300);
