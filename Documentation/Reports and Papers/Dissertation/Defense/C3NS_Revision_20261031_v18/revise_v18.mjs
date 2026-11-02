import fs from "node:fs/promises";
import fsSync from "node:fs";
import path from "node:path";
import { FileBlob, PresentationFile } from "file:///C:/Users/Ben%20Bolen/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/@oai/artifact-tool/dist/artifact_tool.mjs";

const root = path.resolve(".");
const input = path.resolve("../C3NS_Revision_20261031_v17/output/Bolen_PhD_Defense_RoseGreen_20261031_v17.pptx");
const fig = (name) => path.join(root, "submitted_figs", name);
const build = path.join(root, "build");
await fs.mkdir(build, { recursive: true });

const deck = await PresentationFile.importPptx(await FileBlob.load(input));
const snapshot = await deck.inspect({
  select: "$.slides[*]",
  include: ["id", "index", "title", "elements.id", "elements.name", "elements.kind", "elements.text", "elements.position"],
});
const records = snapshot.records;
const byIndex = new Map(records.map((r) => [r.index, r]));

const DARK = "#174D24";
const LEAF = "#58D51B";
const INK = "#253543";
const MUTED = "#60706A";
const PALE = "#EAF4EC";
const LIGHT = "#F5F8F5";

function titleShape(index) {
  const record = byIndex.get(index);
  const candidate = record.elements.find((e) => e.kind === "textbox" && e.position?.top < 40 && e.position?.width > 650);
  if (!candidate) throw new Error(`No title found on slide ${index + 1}`);
  return deck.resolve(candidate.id);
}

function setTitle(index, text) {
  const t = titleShape(index);
  t.text = text;
  t.text.style = {
    fontFamily: "Arial",
    fontSize: text.length > 68 ? 30 : 33,
    bold: true,
    color: DARK,
  };
}

function clearBody(index) {
  const record = byIndex.get(index);
  const ids = record.elements
    .filter((e) => e.position && e.position.top >= 88 && e.position.top < 478)
    .filter((e) => ["textbox", "shape", "image", "chart", "table", "embeddedArtifact"].includes(e.kind))
    .map((e) => e.id);
  if (ids.length) deck.delete(ids);
}

function addImage(slide, imagePath, position, alt, options = {}) {
  const ext = path.extname(imagePath).toLowerCase();
  const contentType = ext === ".jpg" || ext === ".jpeg" ? "image/jpeg" : ext === ".webp" ? "image/webp" : "image/png";
  return slide.images.add({
    blob: fsSync.readFileSync(imagePath),
    contentType,
    alt,
    fit: options.fit ?? "contain",
    position,
    ...(options.crop ? { crop: options.crop } : {}),
    ...(options.geometry ? { geometry: options.geometry } : {}),
    ...(options.borderRadius ? { borderRadius: options.borderRadius } : {}),
  });
}

function addText(slide, text, position, options = {}) {
  const s = slide.shapes.add({
    geometry: options.geometry ?? "textbox",
    name: options.name ?? "v18 text",
    position,
    fill: options.fill ?? "none",
    line: options.line ?? { style: "solid", fill: "none", width: 0 },
    ...(options.borderRadius ? { borderRadius: options.borderRadius } : {}),
  });
  s.text = text;
  s.text.insets = options.insets ?? 4;
  s.text.style = {
    fontFamily: "Arial",
    fontSize: options.fontSize ?? 27,
    bold: options.bold ?? false,
    color: options.color ?? INK,
    alignment: options.alignment ?? "center",
  };
  return s;
}

function addLabel(slide, text, position, options = {}) {
  return addText(slide, text, position, {
    fontSize: options.fontSize ?? 25,
    bold: options.bold ?? true,
    color: options.color ?? DARK,
    fill: options.fill ?? "none",
    line: options.line ?? { style: "solid", fill: "none", width: 0 },
    borderRadius: options.borderRadius,
    insets: options.insets ?? 3,
    alignment: options.alignment ?? "center",
    name: options.name ?? "v18 label",
  });
}

function addCallout(slide, text, position, options = {}) {
  return addText(slide, text, position, {
    geometry: "roundRect",
    fontSize: options.fontSize ?? 27,
    bold: options.bold ?? true,
    color: options.color ?? DARK,
    fill: options.fill ?? PALE,
    line: options.line ?? { style: "solid", fill: options.lineColor ?? DARK, width: options.lineWidth ?? 1.5 },
    borderRadius: options.borderRadius ?? 12,
    insets: options.insets ?? 8,
    alignment: options.alignment ?? "center",
    name: options.name ?? "v18 callout",
  });
}

function setNotes(index, notes) {
  const slide = deck.slides.getItem(index);
  slide.speakerNotes.text = notes;
  slide.speakerNotes.setVisible(true);
}

function setSource(index, text) {
  const record = byIndex.get(index);
  const candidate = record.elements.find((e) => e.kind === "textbox" && e.position?.top >= 470 && e.position?.top < 510 && e.position?.width > 700);
  if (candidate) deck.resolve(candidate.id).text = text;
}

// Slide 3: the biological structure and the user's robot project are introduced visually.
{
  const i = 2, slide = deck.slides.getItem(i);
  clearBody(i);
  setTitle(i, "Human-scale musculoskeletal robots embody biological structure");
  addImage(slide, fig("kengoro_humanoid_comparison.png"), { left: 42, top: 112, width: 470, height: 305 }, "Kengoro human-mimetic humanoid in two poses");
  addImage(slide, fig("aarl_bpa_biped.png"), { left: 620, top: 112, width: 206, height: 305 }, "AARL braided-pneumatic-actuator biped prototype");
  addLabel(slide, "Human-mimetic tendon routing", { left: 72, top: 425, width: 410, height: 40 }, { fontSize: 25 });
  addLabel(slide, "AARL BPA biped", { left: 592, top: 425, width: 262, height: 40 }, { fontSize: 25 });
  setSource(i, "Submitted dissertation Ch. 2, Biomimetic Robotics");
  setNotes(i, "Human-scale musculoskeletal robots embody biological structure in a controllable model. Kengoro demonstrates the broader human-mimetic approach; the AARL BPA biped is the project platform that motivates this dissertation. The robot is not claimed to reproduce every biological degree of freedom or muscle property. Sources: submitted dissertation Background; Asano et al. (2017); AARL robot photograph.");
}

// Slide 4: use the submitted actuator comparison, cropped to the relevant mechanisms.
{
  const i = 3, slide = deck.slides.getItem(i);
  clearBody(i);
  setTitle(i, "Braided pneumatic actuators trade simple rotation for intrinsic compliance");
  addImage(slide, fig("actuator_technology_comparison.png"), { left: 55, top: 112, width: 850, height: 305 }, "Electric torque motor and braided pneumatic actuator mechanisms", { fit: "cover", crop: { left: 0, top: 0, right: 0, bottom: 0.50 } });
  addCallout(slide, "Pressure changes length and force", { left: 485, top: 421, width: 395, height: 44 }, { fontSize: 21 });
  setSource(i, "Submitted dissertation Ch. 2, Artificial Muscles");
  setNotes(i, "The electric motor converts electromagnetic torque to shaft rotation. A braided pneumatic actuator expands radially and contracts axially under pressure, so compliance and the force-length relation are part of the actuator itself. The dissertation therefore begins by characterizing force across resting length, pressure, and contraction. Sources: submitted dissertation Background; Suzumori et al. (2018); Liang et al. (2020).");
}

// Slide 6: make the replacement/ethics motivation concrete.
{
  const i = 5, slide = deck.slides.getItem(i);
  clearBody(i);
  setTitle(i, "Robotic models can replace invasive interventions with repeatable ones");
  addImage(slide, fig("decerebrate_cat_colorized.png"), { left: 38, top: 116, width: 555, height: 320 }, "Archival decerebrate-cat treadmill experiment", { fit: "cover", crop: { left: 0.02, top: 0.01, right: 0.04, bottom: 0.01 } });
  addImage(slide, fig("W2LNEW_TREADMILL_STILL.png"), { left: 640, top: 142, width: 250, height: 270 }, "Supported robotic walker model used for repeatable interventions", { fit: "cover", crop: { left: 0.12, top: 0.13, right: 0.12, bottom: 0.06 } });
  addLabel(slide, "Biology constrains the intervention", { left: 70, top: 439, width: 490, height: 38 }, { fontSize: 24, color: MUTED });
  addLabel(slide, "Model: repeat, remove, stimulate", { left: 600, top: 439, width: 330, height: 38 }, { fontSize: 24 });
  setSource(i, "Submitted dissertation Ch. 2, Ethical Case for Robotic Experiments");
  setNotes(i, "This archival preparation makes the ethical motivation concrete. Destructive pathway or muscle interventions cannot be performed as controlled experiments in healthy people, and animal experiments carry real moral cost. A validated robotic or computational model allows the same initial condition to be reused while one pathway is stimulated, removed, or changed. The archival frame is restored and machine-colorized; it is documentary context, not gait evidence. Source: submitted dissertation Background, ethical case; archival decerebrate-cat footage citation in the dissertation.");
}

// Slide 8: recurring design-loop overview from the submitted dissertation.
{
  const i = 7, slide = deck.slides.getItem(i);
  clearBody(i);
  setTitle(i, "Four contributions close one biomimetic design loop");
  const xs = [25, 255, 485, 715];
  const overview = [
    ["Gait2392_human_reference_oblique.png", "Human\nbenchmarks", { left: 0.22, top: 0.08, right: 0.18, bottom: 0.10 }],
    ["festo_bpas_inflated_deflated.jpg", "Actuator\nforce", undefined],
    ["testJigs_reviewed.png", "Joint\ntransmission", undefined],
    ["W2LNEW_TREADMILL_STILL.png", "Closed-loop\nexperiment", { left: 0.08, top: 0.13, right: 0.08, bottom: 0.05 }],
  ];
  for (let k = 0; k < 4; k++) {
    addImage(slide, fig(overview[k][0]), { left: xs[k], top: 122, width: 205, height: 190 }, overview[k][1].replace("\n", " "), { fit: "cover", crop: overview[k][2] });
    addLabel(slide, overview[k][1], { left: xs[k], top: 320, width: 205, height: 74 }, { fontSize: 27 });
    if (k < 3) addText(slide, "→", { left: xs[k] + 202, top: 218, width: 30, height: 44 }, { fontSize: 34, bold: true, color: LEAF });
  }
  addCallout(slide, "benchmark → actuator → joint → closed loop", { left: 190, top: 414, width: 580, height: 44 }, { fontSize: 22, fill: LIGHT, lineColor: LEAF });
  setSource(i, "Submitted dissertation overview and design workflow");
  setNotes(i, "This workflow is the recurring map for the talk. Human torque and range-of-motion targets guide actuator selection. The characterized BPA model supplies actuator force, and the corrected torque model supports route design. A route that meets the static targets can then move into closed-loop neuromechanical verification before hardware commitment. Source: submitted dissertation Discussion, design workflow.");
}

// Slide 11: independent gaps, not a misleading mechanics-to-database chain.
{
  const i = 10, slide = deck.slides.getItem(i);
  clearBody(i);
  setTitle(i, "Four independent gaps limit biomimetic robot design");
  const xs = [25, 255, 485, 715];
  const imgs = [
    ["festo_bpas_inflated_deflated.jpg", "Contracted and uncontracted BPAs"],
    ["testJigs_reviewed.png", "Pinned-knee and biomimetic-knee test configurations"],
    ["sadb_bubble.png", "Sensory Afferent Database citation map"],
    ["W2LNEW_TREADMILL_STILL.png", "Supported W2L treadmill model"],
  ];
  const labels = [
    "1  Force vs.\nresting length",
    "2  Loaded route\n→ joint torque",
    "3  Scattered\nsensory evidence",
    "4  Controller/body\nexperiments",
  ];
  for (let k = 0; k < 4; k++) {
    addImage(slide, fig(imgs[k][0]), { left: xs[k], top: 120, width: 205, height: 145 }, imgs[k][1], { fit: "cover", crop: k === 3 ? { left: 0, top: 0.12, right: 0, bottom: 0.02 } : undefined });
    addLabel(slide, labels[k], { left: xs[k], top: 278, width: 205, height: 82 }, { fontSize: 26 });
  }
  addCallout(slide, "Main talk: force → torque → neuromechanical experiments", { left: 155, top: 427, width: 650, height: 42 }, { fontSize: 21, fill: LIGHT, lineColor: LEAF, lineWidth: 1.2 });
  setSource(i, "Submitted dissertation research contributions");
  setNotes(i, "These are four independent knowledge or tool gaps, not a claim that SADb connects actuator mechanics to neural control. The main talk follows actuator force, loaded joint transmission, and neuromechanical experiments. SADb addresses the separate problem of dispersed sensory evidence and is retained at the end of the question slides. Source: submitted dissertation research contributions and Discussion.");
}

// Slide 44: submitted numerical corrections and a more informative assertion title.
{
  const i = 43;
  setTitle(i, "Route redesign raises mean torque coverage from 37.35% to 65.81%");
  const chart = deck.resolve("ch/t03uxwrq");
  chart.series.getItemAt(0).values = [37.35, 93.00];
  chart.series.getItemAt(1).values = [65.81, 99.31];
  chart.series.getItemAt(0).dataLabels = { showValue: true, position: "inEnd", numberFormatCode: '0.00"%"', numberFormatSourceLinked: false, textStyle: { fontSize: 14, fill: "#111111", bold: true } };
  chart.series.getItemAt(1).dataLabels = { showValue: true, position: "inEnd", numberFormatCode: '0.00"%"', numberFormatSourceLinked: false, textStyle: { fontSize: 14, fill: "#FFFFFF", bold: true } };
  setNotes(i, "The submitted dissertation corrects the whole-leg coverage summary. Across sampled joint-angle conditions, mean torque coverage rises from 37.35% to 65.81%. Mean active-force coverage rises from 93.00% to 99.31%. Positive actuator force alone is not the torque requirement; the route must also produce the required direction and moment arm. Source: submitted dissertation Chapter 4, Table 4.12, and Appendix C.");
}

// Slide 49: visual section opener.
{
  const i = 48, slide = deck.slides.getItem(i);
  clearBody(i);
  setTitle(i, "Neuromechanical Models and Tools");
  addImage(slide, fig("W2LNEW_TREADMILL_STILL.png"), { left: 135, top: 108, width: 690, height: 360 }, "Supported W2L robot model on a treadmill", { fit: "cover", crop: { left: 0, top: 0.105, right: 0, bottom: 0.03 } });
  addCallout(slide, "controller ↔ body ↔ feedback", { left: 270, top: 397, width: 420, height: 46 }, { fontSize: 24, fill: "#FFFFFF/90", lineColor: LEAF });
  setSource(i, "Submitted dissertation Ch. 6, Neuromechanical Models and Tools");
  setNotes(i, "The mechanical results supply a body that can be tested. The neuromechanical contribution asks how an inspectable neural controller, actuator law, sensory signals, mechanics, and support conditions can be assembled and changed in controlled experiments. The image is the newer supported W2L implementation; its support conditions are stated explicitly. Source: submitted dissertation Chapter 6.");
}

// Slide 53: the newer controller architecture, following the old-model videos.
{
  const i = 52, slide = deck.slides.getItem(i);
  clearBody(i);
  setTitle(i, "The updated W2L controller separates bilateral rhythm and sensory pathways");
  addImage(slide, fig("W2L108_MAIN_SUBNETWORKS.png"), { left: 55, top: 112, width: 850, height: 282 }, "Main subnetworks of the 108-neuron bilateral W2L controller", { fit: "cover", crop: { left: 0, top: 0, right: 0, bottom: 0.55 } });
  addCallout(slide, "108 neurons · 216 connections · contact + type-II feedback", { left: 130, top: 410, width: 700, height: 46 }, { fontSize: 22 });
  setSource(i, "Submitted dissertation Ch. 6, Bilateral W2L Controller");
  setNotes(i, "After the original shared-rhythm-generator transfer, I separated the left and right rhythm generators and linked them through commissural pathways. The revised controller retains bilateral pattern-formation and motor organization while adding distinct heel and toe contact channels and separate type-II relays for flexors and extensors at the hips, knees, and ankles. The full selected-subnetwork figure is retained in the question slides. Source: submitted dissertation Chapter 6, W2L controller.");
}

// Slide 54: the newer supported behavior, with support quantified.
{
  const i = 53, slide = deck.slides.getItem(i);
  clearBody(i);
  setTitle(i, "The updated W2L model alternates under quantified mechanical support");
  addImage(slide, fig("W2LNEW_TREADMILL_STILL.png"), { left: 42, top: 118, width: 635, height: 330 }, "Newer W2L model stepping on a MuJoCo treadmill", { fit: "cover", crop: { left: 0, top: 0.11, right: 0, bottom: 0.02 } });
  addCallout(slide, "9 L / 8 R", { left: 705, top: 145, width: 200, height: 72 }, { fontSize: 36 });
  addLabel(slide, "qualifying swings", { left: 690, top: 220, width: 230, height: 38 }, { fontSize: 22 });
  addCallout(slide, "40.5% BW", { left: 705, top: 280, width: 200, height: 72 }, { fontSize: 34, fill: "#F4F8E8", lineColor: LEAF });
  addLabel(slide, "mean net assistance", { left: 680, top: 356, width: 250, height: 38 }, { fontSize: 21 });
  addLabel(slide, "Passive rails · 30 mm guides", { left: 680, top: 410, width: 250, height: 40 }, { fontSize: 21, color: MUTED });
  setSource(i, "Submitted dissertation Ch. 6, Supported W2L Treadmill Test");
  setNotes(i, "The newer supported MuJoCo model produced nine left and eight right qualifying swings over ten seconds. Passive rails and 30 millimeter guides retained support, and mean net assistance was 40.5 percent of body weight. This is supported treadmill stepping, not harness-free walking. Source: submitted dissertation Chapter 6, supported W2L treadmill stepping.");
}

// Slide 55: refresh the paired human/robot model images from the submitted dissertation.
{
  const i = 54, slide = deck.slides.getItem(i);
  clearBody(i);
  setTitle(i, "A common anatomical frame supports human and BPA-driven model comparisons");
  addImage(slide, fig("Gait2392_human_reference_oblique.png"), { left: 110, top: 112, width: 285, height: 300 }, "Gait2392 human reference with muscle paths", { fit: "cover", crop: { left: 0.22, top: 0.03, right: 0.20, bottom: 0.02 } });
  addImage(slide, fig("mujoco_converted_robot.png"), { left: 560, top: 112, width: 285, height: 300 }, "Converted robot model in MuJoCo", { fit: "cover", crop: { left: 0.22, top: 0.03, right: 0.18, bottom: 0.02 } });
  addLabel(slide, "Human reference", { left: 110, top: 420, width: 285, height: 42 }, { fontSize: 25 });
  addLabel(slide, "Converted robot model", { left: 560, top: 420, width: 285, height: 42 }, { fontSize: 25 });
  setSource(i, "Submitted dissertation Ch. 6, Gait2392-Based Models");
  setNotes(i, "The Gait2392 human model and converted robot model use related anatomical coordinates and muscle-route descriptions. That common frame supports controlled comparisons, but conversion and route repairs require independent checks. It does not imply complete equivalence between engines or an unsupported gait result. Source: submitted dissertation Chapter 6 and Appendix D.");
}

// Slide 56: simplify the feedback ablation to the evidence that supports the claim.
{
  const i = 55, slide = deck.slides.getItem(i);
  clearBody(i);
  setTitle(i, "Feedback keeps bilateral foot loading near antiphase");
  addImage(slide, fig("w2l_feedback_short_20261031.png"), { left: 90, top: 145, width: 780, height: 205 }, "Heel-force traces with W2L feedback intact and disabled", { fit: "contain", crop: { left: 0, top: 0.78, right: 0, bottom: 0.01 } });
  addCallout(slide, "Intact\n0.498 cycle", { left: 160, top: 348, width: 245, height: 78 }, { fontSize: 26 });
  addCallout(slide, "Disabled\n0.223 cycle", { left: 555, top: 348, width: 245, height: 78 }, { fontSize: 26, fill: "#FFF3ED", lineColor: "#B34B24" });
  addText(slide, "Accepted interval: 0.35–0.65 cycle", { left: 275, top: 432, width: 410, height: 32 }, { fontSize: 20, color: MUTED });
  setSource(i, "Submitted dissertation Ch. 6, W2L Sensory-Feedback Comparison");
  setNotes(i, "With the saved controller, drive, and support held fixed, disabling all sensory feedback shifts median bilateral contact phase from 0.498 to 0.223 cycles. Cycling and forward progression persist, but the feedback-disabled condition falls outside the accepted 0.35 to 0.65 cycle interval. This isolates a contribution of the specified sensory pathways to bilateral contact coordination. Source: submitted dissertation Chapter 6, feedback comparison.");
}

// Slide 57: cross-implementation verification before coupled interpretation.
{
  const i = 56, slide = deck.slides.getItem(i);
  clearBody(i);
  setTitle(i, "NumPy and Simulink reproduce the same W2L rhythm core");
  addImage(slide, fig("sns_w2l_traces_short_20261031.png"), { left: 95, top: 112, width: 770, height: 300 }, "Matched NumPy and Simulink rhythm-generator voltage traces");
  addCallout(slide, "1.000 s period · all 10 bursts · r ≥ 0.99996", { left: 190, top: 416, width: 580, height: 46 }, { fontSize: 22, fill: "#FFFFFF/92", lineColor: LEAF });
  setSource(i, "Submitted dissertation Ch. 6, Matched-Input W2L Verification");
  setNotes(i, "Before interpreting coupled motion, the contact-driven W2L rhythm core was replayed with the same neuron definitions, connections, parameters, drive, and contact inputs. NumPy and Simulink matched the one-second burst period, reproduced all ten expected bursts, and achieved voltage-trace correlation of at least 0.99996. Source: submitted dissertation Chapter 6, matched-input W2L comparison.");
}

// Slide 58: crop to the readable functional catalog, not the paper-scale equation panel.
{
  const i = 57, slide = deck.slides.getItem(i);
  clearBody(i);
  setTitle(i, "Twelve editable blocks expose neural, sensory, and BPA definitions");
  addImage(slide, fig("sns_library_clean.png"), { left: 70, top: 120, width: 820, height: 280 }, "Twelve implemented SNS-Simscape blocks", { fit: "cover", crop: { left: 0, top: 0, right: 0, bottom: 0.60 } });
  addCallout(slide, "neurons + synapses + afferents + activation + BPA force", { left: 125, top: 414, width: 710, height: 46 }, { fontSize: 22 });
  setSource(i, "Submitted dissertation Ch. 6, SNS-Simscape Block Library");
  setNotes(i, "The graphical library exposes the parts that were previously buried inside one controller implementation. The twelve blocks cover graded and spiking neurons, synapses, summation, Ia and Ib afferents, activation, biological muscle, BPA force, and diameter-specific BPA variants. The underlying neuron formulation is retained from SNS-Toolbox; the contribution is the editable implementation and integration workflow. Source: submitted dissertation Chapter 6, SNS-Simscape library.");
}

// Slide 59: one result, large enough to read from the auditorium.
{
  const i = 58, slide = deck.slides.getItem(i);
  clearBody(i);
  setTitle(i, "Sensory feedback cuts loaded-arm sag from 7.85° to 2.00°");
  addImage(slide, fig("beverage_current_layout.png"), { left: 75, top: 112, width: 810, height: 350 }, "Loaded-arm feedback demonstration", { fit: "cover", crop: { left: 0, top: 0, right: 0, bottom: 0.40 } });
  setSource(i, "Submitted dissertation Ch. 6, Coupled Reflex Verification");
  setNotes(i, "The loaded-arm verification checks the full sensory-to-motor path. With feedback enabled, the forearm sags by 2.00 degrees as the cup is loaded; with feedback disabled, it sags by 7.85 degrees. This demonstrates a corrective effect in this specified coupled test, not a general balance or locomotion result. Source: submitted dissertation Chapter 6, loaded-arm demonstration.");
}

// Slide 60: replace generic boxes with the actual construction and interface workflow.
{
  const i = 59, slide = deck.slides.getItem(i);
  clearBody(i);
  setTitle(i, "One workflow connects CAD, OpenSim, Simscape, MuJoCo, and SNS tools");
  addImage(slide, fig("ah001_supported_control.png"), { left: 40, top: 125, width: 260, height: 260 }, "SOLIDWORKS-derived robot model in a supported test", { fit: "cover", crop: { left: 0.10, top: 0.04, right: 0.10, bottom: 0.02 } });
  addImage(slide, fig("w2l_selected_left_subnetwork.png"), { left: 350, top: 125, width: 260, height: 260 }, "Selected neural controller subnetwork", { fit: "cover", crop: { left: 0.12, top: 0.08, right: 0.10, bottom: 0.08 } });
  addImage(slide, fig("W2LNEW_TREADMILL_STILL.png"), { left: 660, top: 125, width: 260, height: 260 }, "Closed-loop walker experiment", { fit: "cover", crop: { left: 0.08, top: 0.13, right: 0.08, bottom: 0.05 } });
  addLabel(slide, "Construct the body", { left: 40, top: 395, width: 260, height: 50 }, { fontSize: 24 });
  addLabel(slide, "Construct the controller", { left: 350, top: 395, width: 260, height: 50 }, { fontSize: 24 });
  addLabel(slide, "Run the closed loop", { left: 660, top: 395, width: 260, height: 50 }, { fontSize: 24 });
  setSource(i, "Submitted dissertation Ch. 6, Model Construction and Controller Interfaces");
  setNotes(i, "The workflow separates model construction from controller experiments. SOLIDWORKS geometry feeds Simscape and separately assembled MuJoCo models. OpenSim follows the MyoConverter and repair route. The SNS-Simulink editor exchanges commands and sensory measurements with Simscape or with MuJoCo through the bridge, while SNS-Toolbox provides a scripted alternative for automated experiments. Source: submitted dissertation Chapter 6, model-construction workflow.");
}

// Slide 64: keep the controlled-intervention result; it is stronger than the prior spiking timing aside.
{
  const i = 63, slide = deck.slides.getItem(i);
  clearBody(i);
  setTitle(i, "Acceleration makes staged neuromechanical experiments practical");
  addCallout(slide, "9.21×", { left: 90, top: 145, width: 235, height: 120 }, { fontSize: 50, fill: PALE, lineColor: DARK });
  addLabel(slide, "Gait2392 coupled loop", { left: 70, top: 275, width: 275, height: 55 }, { fontSize: 25 });
  addCallout(slide, "17.75×", { left: 362, top: 145, width: 235, height: 120 }, { fontSize: 50, fill: "#F4F8E8", lineColor: LEAF });
  addLabel(slide, "W2L coupled loop", { left: 342, top: 275, width: 275, height: 55 }, { fontSize: 25 });
  addCallout(slide, "100.17×", { left: 635, top: 145, width: 235, height: 120 }, { fontSize: 50, fill: LIGHT, lineColor: MUTED });
  addLabel(slide, "1,027-neuron core", { left: 615, top: 275, width: 275, height: 55 }, { fontSize: 25 });
  addText(slide, "Matched inputs after warm-up · rendering and file output excluded", { left: 155, top: 350, width: 650, height: 40 }, { fontSize: 22, color: MUTED });
  addLabel(slide, "Complete W2L loop: 4.13× faster than real time", { left: 200, top: 414, width: 560, height: 44 }, { fontSize: 24 });
  setSource(i, "Submitted dissertation Ch. 6, Execution Speed");
  setNotes(i, "The acceleration work makes staged testing and controller searches practical. After warm-up, the complete Gait2392 balance loop accelerated 9.21 times, the W2L80 suspended loop 17.75 times, and a 1,027-neuron neural-only benchmark 100.17 times. The complete accelerated W2L loop ran about 4.13 times faster than real time. Loading, compilation, file output, and rendering are excluded. Source: submitted dissertation Chapter 6, runtime table.");
}

// Slide 65: return to the recurring overview and highlight the enabled verification step.
{
  const i = 64, slide = deck.slides.getItem(i);
  clearBody(i);
  setTitle(i, "The contributions shorten the path from benchmark to controlled experiment");
  addImage(slide, fig("design_workflow.png"), { left: 100, top: 108, width: 760, height: 350 }, "Design workflow enabled by the dissertation contributions", { fit: "contain", crop: { left: 0, top: 0, right: 0, bottom: 0.18 } });
  setSource(i, "Submitted dissertation Discussion, Design Workflow");
  setNotes(i, "The actuator work improves first-pass sizing, the torque corrections improve route prediction, and the neuromechanical tools make closed-loop verification and controlled interventions available before hardware commitment. SADb supports the separate evidence-selection problem and remains in the question slides. The highlighted box marks the step enabled by the Chapter 6 tools. Source: submitted dissertation Discussion, design workflow.");
}

// Slide 69: use the submitted simulation-to-hardware workflow instead of abstract boxes.
{
  const i = 68, slide = deck.slides.getItem(i);
  clearBody(i);
  setTitle(i, "The implemented tools form a simulation-to-hardware path");
  addImage(slide, fig("pipeline.png"), { left: 80, top: 110, width: 800, height: 360 }, "Model conversion, closed-loop simulation, and proposed hardware deployment");
  setSource(i, "Submitted dissertation Discussion, Simulation-to-Hardware Workflow");
  setNotes(i, "The implemented path covers model conversion and closed-loop simulation. The future hardware branch would connect neural evaluation on the embedded computer to joint pressure and valve control, then return encoder, pressure, force, contact, and inertial measurements. Computer-in-the-loop testing should precede treadmill deployment. The dashed branch is proposed future work, not a completed hardware controller. Source: submitted dissertation Discussion, simulation-to-hardware workflow.");
}

// Export focused previews, a montage, and the editable draft.
const changed = [2, 3, 5, 7, 10, 43, 48, 52, 53, 54, 55, 56, 57, 58, 59, 63, 64, 68];
for (const index of changed) {
  const preview = await deck.render({ slide: index, scale: 1 });
  await fs.writeFile(path.join(build, `artifact-slide-${String(index + 1).padStart(3, "0")}.png`), new Uint8Array(await preview.arrayBuffer()));
  await fs.writeFile(path.join(build, `artifact-slide-${String(index + 1).padStart(3, "0")}.layout.json`), await preview.layout.text(), "utf8");
}
const montage = await deck.export({ format: "webp", montage: { columns: 5, slideWidth: 300, padding: 12, gap: 12 } });
await fs.writeFile(path.join(build, "v18-artifact-montage.webp"), new Uint8Array(await montage.arrayBuffer()));
const pptx = await PresentationFile.exportPptx(deck);
await pptx.save(path.join(build, "Bolen_PhD_Defense_RoseGreen_20261031_v18_artifact_draft.pptx"));
console.log(`Edited ${changed.length} slides; exported ${deck.slides.items.length} slides.`);
