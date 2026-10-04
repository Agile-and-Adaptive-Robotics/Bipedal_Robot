// Node check that mirrors the editor's loadTplJson + renderTplOptions
// path for the new spiking_mirror entry: the embedded TPL_DATA_EMBED must
// carry it, TPL_JSON_GROUP/NAMES must resolve it into the walker group,
// and the spec must survive the same shape validation the editor's
// loadSpec path assumes (nodes[{type,label,x,y}], edges, groups, unique
// labels, no dangling edge endpoints).
const fs = require("fs");
const src = fs.readFileSync("connectome_block_editor.html", "utf8");

const embed = src.slice(src.indexOf("/*TPL-EMBED-START*/"),
                        src.indexOf("/*TPL-EMBED-END*/"));
const json = embed.slice(embed.indexOf("{"), embed.lastIndexOf(";"));
const data = JSON.parse(json);

const js = src.slice(src.lastIndexOf("<script>") + 8,
                     src.lastIndexOf("</script>"));
const grpSrc = js.slice(js.indexOf("const TPL_JSON_GROUP"),
                        js.indexOf("let TPL_SRC"))
  .replace("const TPL_JSON_GROUP", "globalThis.TPL_JSON_GROUP")
  .replace("const TPL_JSON_NAMES", "globalThis.TPL_JSON_NAMES");
eval(grpSrc);
const TPL_JSON_GROUP = globalThis.TPL_JSON_GROUP;
const TPL_JSON_NAMES = globalThis.TPL_JSON_NAMES;

let fail = 0;
const T = (name, cond) => {
  console.log((cond ? "PASS " : "FAIL ") + name);
  if (!cond) fail++;
};

T("embed carries spiking_mirror", !!data.spiking_mirror);
T("TPL_JSON_GROUP maps it to walker", TPL_JSON_GROUP.spiking_mirror === "walker");
T("TPL_JSON_NAMES has a label", typeof TPL_JSON_NAMES.spiking_mirror === "string"
  && TPL_JSON_NAMES.spiking_mirror.length > 0);
T("sidecar file matches embed key set",
  JSON.stringify(Object.keys(data).sort()) ===
  JSON.stringify(Object.keys(JSON.parse(
    fs.readFileSync("connectome_templates.json", "utf8"))).sort()));

const sm = data.spiking_mirror;
const labels = new Set(sm.nodes.map(n => n.label));
T("unique labels", labels.size === sm.nodes.length);
T("no dangling edges",
   sm.edges.every(e => labels.has(e.from) && labels.has(e.to)));
T("groups cover 19 layers", Object.keys(sm.groups).length === 19);
const grpSet = new Set(sm.nodes.map(n => n.grp).filter(Boolean));
T("every grp id exists in groups dict",
   [...grpSet].every(g => sm.groups[g]));
const nGrp = sm.nodes.filter(n => n.grp).length;
console.log(`spiking_mirror: ${sm.nodes.length} nodes (${nGrp} stamped), ` +
            `${sm.edges.length} edges, ${Object.keys(sm.groups).length} layers`);
process.exit(fail);
