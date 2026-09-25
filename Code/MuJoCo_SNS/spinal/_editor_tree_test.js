// Unit test of the v3.2 layer stamping rules (stampWalkerGrps) extracted
// from connectome_block_editor.html — mirrors stamp_walker_grps in
// make_editor_templates.py. Locks the label->layer mapping.
const fs = require("fs");
const src = fs.readFileSync("connectome_block_editor.html", "utf8");
const js = src.slice(src.lastIndexOf("<script>") + 8,
                     src.lastIndexOf("</script>"));
const s0 = js.indexOf("function stampWalkerGrps");
const s1 = js.indexOf("const TEMPLATES = {");
eval(js.slice(s0, s1));

let pass = 0, fail = 0;
const T = (name, cond) => {
  if (cond) { pass++; console.log("PASS " + name); }
  else { fail++; console.log("FAIL " + name); }
};
const stamp = labels => stampWalkerGrps({
  nodes: labels.map(l => ({ label: l })) });

// walker-style labels (incl. pruned MNs, sig nodes, heel_IN)
let r = stamp(["DRIVE", "POSTURE", "pf_gain=1.042", "PM T=0.31 ws=0.28",
  "RG_E_R", "InE_R", "c1_R", "V3_L", "heel_IN_R",
  "PF_E1_R", "PF_F2_L",
  "MN_glut_max1_R", "MN_quad_fem_R (pruned)", "ext MN_R",
  "glut_max1_R", "soleus_L",
  "Ia_hip_R", "II_knee_L", "Ib_ankle_R", "heel_R", "toe_L",
  "hip_ext_sig_R",
  "KINH_R", "PRESET_E_L", "IaIN_R", "LBIN_L", "IBEXC_R", "RC_L",
  "IIe_ext_R", "IIi_flex_L", "IbIN_ext_R"]);
const g = Object.fromEntries(r.nodes.map(n => [n.label, n.grp]));
T("ports/PM -> drive", g["DRIVE"] === "drive" && g["POSTURE"] === "drive" &&
  g["pf_gain=1.042"] === "drive" && g["PM T=0.31 ws=0.28"] === "drive");
T("RG/lamination/c1/V3/heel_IN -> rg_S",
  g["RG_E_R"] === "rg_R" && g["InE_R"] === "rg_R" &&
  g["c1_R"] === "rg_R" && g["V3_L"] === "rg_L" &&
  g["heel_IN_R"] === "rg_R");
T("PF -> pf_S", g["PF_E1_R"] === "pf_R" && g["PF_F2_L"] === "pf_L");
T("MNs (incl pruned) -> mn_S",
  g["MN_glut_max1_R"] === "mn_R" &&
  g["MN_quad_fem_R (pruned)"] === "mn_R" &&
  g["ext MN_R"] === "mn_R");
T("muscles -> mus_S", g["glut_max1_R"] === "mus_R" &&
  g["soleus_L"] === "mus_L");
T("afferents -> aff_S", g["Ia_hip_R"] === "aff_R" &&
  g["II_knee_L"] === "aff_L" && g["Ib_ankle_R"] === "aff_R" &&
  g["heel_R"] === "aff_R" && g["toe_L"] === "aff_L" &&
  g["hip_ext_sig_R"] === "aff_R");
T("reflex/phase INs -> in_S", g["KINH_R"] === "in_R" &&
  g["PRESET_E_L"] === "in_L" && g["IaIN_R"] === "in_R" &&
  g["LBIN_L"] === "in_L" && g["IBEXC_R"] === "in_R" &&
  g["RC_L"] === "in_L" && g["IIe_ext_R"] === "in_R" &&
  g["IIi_flex_L"] === "in_L" && g["IbIN_ext_R"] === "in_R");
T("groups dict carries all 13 layers",
  Object.keys(r.groups).length === 13 && !!r.groups.mn_R);

// ungrouped specs (replication/rules templates) stay untouched -> the
// editor falls back to type buckets. Side-suffixed unknowns are
// deliberately muscles (how new muscle names get layered); a label
// with NO side suffix is the fallback case.
let q = stamp(["V0C_interneuron"]);
T("non-walker label left unstamped (type-bucket fallback)",
  q.nodes[0].grp === undefined);

console.log(`\n${pass} pass, ${fail} fail`);
process.exit(fail ? 1 : 0);
