// Functional test of the editor's template builders (pure functions —
// no DOM). Extracts litSide/buildBilateral/buildW2L/buildDengPair from
// the HTML and validates each output's shape against loadSpec's
// expectations: nodes[{type,label,x,y}] + synapses|edges[{from,to,
// sign,gain,tag}].
const fs = require("fs");
const src = fs.readFileSync("connectome_block_editor.html", "utf8");
// NB: since 2026-09-24 the file carries an EMBEDDED template <script>
// before the main one — slice the LAST script block, not the first
// closer (indexOf("</script>") lands inside the embed and empties js).
const js = src.slice(src.lastIndexOf("<script>") + 8,
                     src.lastIndexOf("</script>"));

// pull out only the template-builder region (from "function litSide"
// to the tplSel listener; excludes the w2lfile XHR template)
const start = js.indexOf("function litSide");
const end = js.indexOf('document.getElementById("tplSel")');
const region = js.slice(start, end);

const TYPES_NAMES = ["SN-Ia","SN-II","SN-Ib","SN-heel","SN-toe",
  "PORT-load","IN-V0D","IN-V0V","IN-V1","IN-V2a","IN-V2b","IN-V3",
  "IN-C","HC-RG-E","HC-RG-F","IN-InE","IN-InF","HC-PF-E","HC-PF-F",
  "IN-PF","IN-IaIN","IN-IbIN","IN-IIe","IN-IIi","IN-Ib+","IN-LBIN",
  "IN-KINH","RC","MN","MUSCLE"];
const TYPES = {};
TYPES_NAMES.forEach(t => TYPES[t] = {shape:"circle", color:"#000",
                                     w:2, grp:"X"});

// (v3.0 template-library IIFEs now ride along inside the region; give
// them just enough browser to no-op gracefully under node. MUST be
// defined before eval — the IIFEs execute during it.)
global.window = global.window || {};
global.status = global.status || (() => {});
global.document = global.document ||
  { getElementById: () => null };
global.XMLHttpRequest = global.XMLHttpRequest || function () {
  this.open = () => {}; this.send = () => { this.status = 404; };
};
eval(region);

const results = { deng: buildDengPair(), lit: litSide("R"),
                  bilateral: buildBilateral(), w2l: buildW2L() };
let fail = 0;
for (const [name, spec] of Object.entries(results)) {
  const syns = spec.synapses || spec.edges || [];
  const badTypes = spec.nodes.filter(n => !TYPES[n.type]);
  const labels = new Set(spec.nodes.map(n => n.label));
  const dangling = syns.filter(s => !labels.has(s.from) ||
                                    !labels.has(s.to));
  const badSign = syns.filter(s => s.sign !== "exc" &&
                                   s.sign !== "inh");
  const ok = spec.nodes.length > 0 && syns.length > 0 &&
             badTypes.length === 0 && dangling.length === 0 &&
             badSign.length === 0;
  console.log(`${ok ? "PASS" : "FAIL"} ${name}: ` +
    `${spec.nodes.length} nodes, ${syns.length} synapses` +
    (badTypes.length ? ` | BAD TYPES: ${badTypes.map(n=>n.type)}` : "") +
    (dangling.length ? ` | DANGLING: ${dangling.length}` : "") +
    (badSign.length ? ` | BAD SIGN: ${badSign.length}` : ""));
  if (!ok) fail++;
}
process.exit(fail);
