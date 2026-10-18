// Isometric torque paper bridge — dissertation chapter 30-torque + Testing_Data archives
// bridged into the Frontiers manuscript at Knee_Torque_Test/ProofFinal.
// Every subagent reads tmp/paper_bridge/BRIEF.md first (rules, paths, figure spec).

interface MinerResult {
  /** Path of the digest markdown this miner wrote. */
  digestPath: string;
  /** Headline findings the planner must not miss. */
  highlights: string[];
  /** Questions only Ben can answer; surfaced in the final report. */
  openQuestions: string[];
}

interface BridgePlan {
  planPath: string;
  /** Paper sections the plan maps (e.g. "Methods", "Results"). */
  sections: string[];
  /** Planned tables (name + planned .xlsx filename). */
  tables: string[];
  /** Stale result figures the plan regenerates (filenames). */
  figuresToRegenerate: string[];
  notes: string[];
}

interface PlanReview {
  approved: boolean;
  findings: { severity: "low" | "medium" | "high"; what: string }[];
}

interface NumberSource {
  /** The load-bearing number as it appears in the fragment. */
  value: string;
  /** Ledger row / chapter location / citation that sources it. */
  source: string;
}

interface SectionDraft {
  section: string;
  fragmentPath: string;
  /** Citations the section needs that are absent from BolenFrontiers22.bib. */
  missingCitations: { key: string; claim: string }[];
  /** The 10 most load-bearing numbers used, each with its named source. */
  numberSources: NumberSource[];
  notes: string[];
}

interface TablesOutcome {
  fragments: string[];
  xlsxFiles: string[];
  numberSources: NumberSource[];
  notes: string[];
}

interface FiguresOutcome {
  regenerated: string[];
  /** Figures skipped because data source or content was ambiguous (Ben decides). */
  skipped: string[];
  scriptPaths: string[];
  notes: string[];
}

interface GateFinding {
  /** Which gate produced it (e.g. "Gate 2 slop"). */
  gate: string;
  /** Location in the manuscript, with a line when it applies. */
  where: string;
  /** One sentence: what failed the gate. */
  what: string;
  /** How much it matters. Reserve "high" for plagiarism, wrong numbers, or directive violations. */
  severity: "low" | "medium" | "high";
}

interface GateVerdict {
  verdict: "GO" | "BLOCK";
  findings: GateFinding[];
}

interface ColdRead {
  unclear: string[];
  unsupported: string[];
  reviewerQuestions: string[];
}

const BRIEF = "tmp/paper_bridge/BRIEF.md";
const MINING = "tmp/paper_bridge/mining";
const SECTIONS = "tmp/paper_bridge/sections";
const MAINTEX = "Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Main.tex";
const SUPPMAT = "Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/frontiers_SupplementaryMaterial.tex";
const MAINAUX = "Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Main.aux";

artifact.table("progress", {
  title: "Bridge progress",
  columns: [
    { field: "item", label: "Workstream" },
    { field: "status", label: "Status" },
    { field: "detail", label: "Detail" },
  ],
  key: "item",
});

const READ_ONLY = " Read-only mission: the only files you may write are inside tmp/paper_bridge/.";
const BRIEF_LINE = `Read ${BRIEF} first and follow it exactly (folder rules, editorial contract, figure spec).`;

phase("Mine the chapter, the data archives, and the current draft");
log("Four miners: dissertation chapter, archived results, Actuators-paper register, current manuscript state.");
const miners = await Promise.all([
  agent("dissertation miner").ask<MinerResult>(
    `${BRIEF_LINE}${READ_ONLY}
Mission: produce the full bridge digest of the committed dissertation torque chapter
Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/30-torque.tex (read the whole file),
its sibling Figures/30-torque/ inventory, and any torque-relevant tables in chapters/92-AppendixA.tex,
93-AppendixB.tex, 94-AppendixC.tex (settled-Xi tables live in an appendix). Record:
(1) section map, marking the publication-note overlap sections (bolen_isometric_2025);
(2) every equation (label + one-line meaning), especially the fictive-strain framing via epsilon_620;
(3) every figure (file, caption gist, what it evidences); (4) every table; (5) numbers of record with
chapter locations; (6) how the Discussion carries the follow-up identification (Xi0-Xi3) and
route-redesign results. Write ${MINING}/chapter_digest.md and return the result.`,
  ),
  agent("results archivist").ask<MinerResult>(
    `${BRIEF_LINE}${READ_ONLY}
Mission: build the numbers ledger for the paper. Inventory Testing_Data/2022_02_Festo (focus on
Dig_out/ result mats and logs and the minimize* evaluator/driver family), read the summary tabs of
Testing_Data/2026_06_Festo/Results_table_10mm_pinned.xlsx and Results_table_20mm.xlsx (use
D:/Anaconda/python.exe with openpyxl, READ-ONLY), and inventory Code/Matlab/Mesh_Optimization/Results
(route-redesign outcomes). For each candidate number: value, unit, source file/tab/sheet, the paper
claim it supports, and confidence (dissertation-record vs single-run artifact). Prefer values as
recorded in the dissertation chapter; flag disagreements instead of averaging. Write
${MINING}/results_ledger.md as a markdown table and return the result.`,
  ),
  agent("actuators-paper analyst").ask<MinerResult>(
    `${BRIEF_LINE}${READ_ONLY}
Mission: distill the register and structure of the published Actuators force paper (bolen_2026) from
Documentation/Reports and Papers/Isometric_Force_Test/ (find and read the manuscript .tex or PDF).
Capture: section skeleton and approximate lengths; how Methods introduces equations; Results voice;
how tables and figures are captioned and referenced; citation density and style; abstract shape.
Write ${MINING}/actuators_style.md as a style card the section writers calibrate against.`,
  ),
  agent("manuscript analyst").ask<MinerResult>(
    `${BRIEF_LINE}${READ_ONLY}
Mission: map the current manuscript state at Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/
(Main.tex whole file; skim frontiers_SupplementaryMaterial.tex; Figures/ inventory; BolenFrontiers22.bib key
list; latexmkrc). Record: sections present vs missing; stale elements (title block, abstract register,
introduction gaps); every \\includegraphics target vs existing files; every \\cite key vs the bib (list
missing ones); class/preamble constraints (what may and may not be added); TODO markers. Write
${MINING}/paper_state.md.`,
  ),
]);
report({ item: "Mining", status: "done", detail: `${miners.length} digests in ${MINING}/` }, "progress");

phase("Draft the bridge plan");
const planner = agent("bridge planner");
const plan = await planner.ask<BridgePlan>(
  `${BRIEF_LINE}
Read the four digests in ${MINING}/. Write the bridge plan at tmp/paper_bridge/plan.md:
(1) section-by-section mapping chapter to paper (what bridges, what is rewritten to journal register,
what is newly written); (2) the long-Discussion outline carrying the follow-up identification and
route-redesign results as ADDITIONAL RESULTS, fictive-strain framing throughout; (3) abstract, title,
and keywords strategy (working title per BRIEF); (4) tables list (name, columns, ledger source rows,
planned .xlsx filename into Tables/); (5) figures list: keep / regenerate / new, with filenames;
(6) citation-gap candidates. Be decisive; where Ben must decide something, put it in an OPEN QUESTIONS
section of the plan rather than guessing.`,
);
report({ item: "Plan", status: "drafted", detail: `${plan.sections.length} sections, ${plan.tables.length} tables, ${plan.figuresToRegenerate.length} figures to regenerate` }, "progress");

phase("Independent review of the bridge plan");
const planReviewer = agent("plan reviewer", {
  system: "You review a bridge plan you did not write. Ask what would break it, not whether it looks fine. Escalate rather than guess.",
});
let planFeedback = "none";
let planApproved = false;
for (let round = 1; round <= 2; round++) {
  const review = await planReviewer.ask<PlanReview>(
    `${BRIEF_LINE}${READ_ONLY}
Read tmp/paper_bridge/plan.md, then spot-check it against the evidence: the dissertation chapter
(chapters/30-torque.tex), the digests in ${MINING}/, and the current ProofFinal/Main.tex. Find what
would break it: missing sections, mis-mapped content, tables or figures without sources, register
risks, anything violating the BRIEF rules. Start your reply with the plan restated in your own words
in two sentences. Previous critique already addressed: ${planFeedback}. Approve only if complete and faithful.`,
  );
  const blocking = review.findings.filter((f) => f.severity === "high");
  if (review.approved || blocking.length === 0) {
    planApproved = true;
    break;
  }
  planFeedback = JSON.stringify(blocking);
  await planner.ask(`Revise tmp/paper_bridge/plan.md to address these blocking findings (round ${round}): ${planFeedback}`);
}
report({ item: "Plan", status: planApproved ? "approved" : "best-effort", detail: "independent plan review complete" }, "progress");

phase("Write the paper sections in parallel");
const writerAsk = (section: string, fragment: string, extra: string) =>
  `${BRIEF_LINE}
Read tmp/paper_bridge/plan.md fully plus the digests you need in ${MINING}/. Write the ${section} of the
paper as a LaTeX fragment to ${SECTIONS}/${fragment} (body content only: \\section/\\subsection, text,
equations, figure and table floats; NO preamble, NO \\begin{document}). Journal register per the BRIEF
and the actuators style card. Cite only keys that exist in ProofFinal/BolenFrontiers22.bib; report any
missing citation instead of inventing it. Every load-bearing number comes from the results ledger or
chapter digest. Return your result with the 10 most load-bearing numbers and their sources.${extra}`;

// Join: the assembler needs every fragment before merging, so all six land together.
const [intro, methods, results, discussion, tables, figures] = await Promise.all([
  agent("introduction writer").ask<SectionDraft>(
    writerAsk("Abstract, keywords, and Introduction", "10_abstract_intro.tex", " The abstract states the design-factor contribution, not a diary."),
  ),
  agent("methods writer").ask<SectionDraft>(
    writerAsk("Methods and Materials", "20_methods.tex", " Reproduce the actuator force model only as briefly as the Actuators paper reproduced its own background; the fictive-strain quantities are defined here."),
  ),
  agent("results writer").ask<SectionDraft>(
    writerAsk("Results", "30_results.tex", " Results are reported, not interpreted; interpretation belongs to the Discussion."),
  ),
  agent("discussion writer").ask<SectionDraft>(
    writerAsk("long Discussion and Conclusion", "40_discussion_conclusion.tex", " The Discussion carries the follow-up identification (Xi0-Xi3) and the route-redesign results as ADDITIONAL RESULTS, per the plan."),
  ),
  agent("tables builder").ask<TablesOutcome>(
    `${BRIEF_LINE}
Read tmp/paper_bridge/plan.md (tables list) and ${MINING}/results_ledger.md. Build each planned table
twice: (a) a LaTeX table fragment per table into ${SECTIONS}/tables_main.tex (one float per table,
booktabs style, Frontiers-friendly); (b) an .xlsx export per table into
Documentation/Reports and Papers/Knee_Torque_Test/Tables/ via D:/Anaconda/python.exe + openpyxl
(export script in tmp/paper_bridge/tabscripts/). ONLY .xlsx files may be created in Tables/. Every
number comes from the ledger or chapter digest. Return fragments, xlsx paths, and number sources.`,
  ),
  agent("figures regenerator").ask<FiguresOutcome>(
    `${BRIEF_LINE}
Read tmp/paper_bridge/plan.md (figures list) and follow the BRIEF figure specification exactly
(FlxPin_group.fig exemplar conventions, Plot_KneeFlx* script settings, PlotPrompt palette,
figure_qc rules). For each figure marked regenerate: write a MATLAB script in
tmp/paper_bridge/figscripts/, run it headless with matlab -batch (addpath Code/Matlab/Functions and
Code/Matlab/Robot_Data as the existing plot scripts require; data from Testing_Data/2022_02_Festo),
and export EPS to Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Figures/ under the
SAME filename it replaces. Do NOT edit any .tex file. If a figure's data source or content is
ambiguous, SKIP it and report why (Ben decides) rather than guessing. Every regenerated figure's
caption must be reminded to state that negative values denote knee flexion; include that reminder
list in your notes.`,
  ),
]);
report({ item: "Sections", status: "drafted", detail: `${[intro, methods, results, discussion].length} fragments; ${tables.xlsxFiles.length} xlsx; ${figures.regenerated.length} figures regenerated, ${figures.skipped.length} skipped` }, "progress");

phase("Assemble the manuscript and compile");
const assembler = agent("manuscript assembler");
const assemblerResult = await assembler.ask<SectionDraft>(
  `${BRIEF_LINE}
Read ${MINING}/paper_state.md and tmp/paper_bridge/plan.md. Assemble the final manuscript by editing
ONLY Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Main.tex: merge the fragments from
${SECTIONS}/ in order, set the working title (long and short) per the BRIEF, place the table floats at
their planned positions, keep the Frontiers class, the author block, and the supplementary material
untouched. Preamble additions allowed only when needed (e.g. mathabx for \\diameter). Every
\\includegraphics must resolve to an existing file in ProofFinal/Figures/; if a figure was skipped by
the figures regenerator, keep the old figure and say so in your notes. Every \\cite key must exist in
ProofFinal/BolenFrontiers22.bib; list any missing keys in missingCitations. Report section-by-section
what you changed and anything preserved-but-stale. Return the load-bearing numbers you placed, with sources.`,
);

const missingCits = [intro, methods, results, discussion, assemblerResult]
  .flatMap((d) => d.missingCitations.map((c) => `${c.key}: ${c.claim}`));
if (missingCits.length > 0) {
  await agent("bibliographer").ask(
    `${BRIEF_LINE}
Verify and add missing bibliography entries for these claims:
${missingCits.join("\n")}
Confirm each reference exists (follow the DOI, publisher page, or arXiv listing) BEFORE writing it;
deliver only verified BibTeX keyed to match the \\cite commands. Append the entries ONLY to
Documentation/Reports and Papers/Knee_Torque_Test/Bibliography/BolenFrontiers22.bib and mirror the
same entries into ProofFinal/BolenFrontiers22.bib. Anything you cannot verify: report it, do not add it.`,
  );
}
report({ item: "Assembly", status: "assembled", detail: missingCits.length > 0 ? `${missingCits.length} bib entries verified/added` : "no bib gaps" }, "progress");

async function compileAll(): Promise<{ ok: boolean; pending: boolean; tail: string }> {
  try {
    await world.run("D:/MiKTeX/miktex/bin/x64/pdflatex.exe", ["-cd", "-interaction=nonstopmode", SUPPMAT], { timeoutMs: 600000 });
    const r1 = await world.run("D:/MiKTeX/miktex/bin/x64/pdflatex.exe", ["-cd", "-interaction=nonstopmode", MAINTEX], { timeoutMs: 900000 });
    if (r1.exitCode !== 0) {
      return { ok: false, pending: false, tail: (r1.stdout + r1.stderr).slice(-4000) };
    }
    await world.run("D:/MiKTeX/miktex/bin/x64/bibtex.exe", [MAINAUX], { timeoutMs: 300000 });
    await world.run("D:/MiKTeX/miktex/bin/x64/pdflatex.exe", ["-cd", "-interaction=nonstopmode", MAINTEX], { timeoutMs: 600000 });
    const r3 = await world.run("D:/MiKTeX/miktex/bin/x64/pdflatex.exe", ["-cd", "-interaction=nonstopmode", MAINTEX], { timeoutMs: 600000 });
    return { ok: r3.exitCode === 0, pending: false, tail: (r3.stdout + r3.stderr).slice(-4000) };
  } catch {
    return { ok: false, pending: true, tail: "" };
  }
}

let compile = await compileAll();
let compileRounds = 0;
while (!compile.ok && !compile.pending && compile.tail !== "" && compileRounds < 2) {
  compileRounds++;
  await assembler.ask(`The local pdflatex compile failed (round ${compileRounds}). Tail of the output:
${compile.tail}
Fix the manuscript per the BRIEF rules (edit only ProofFinal files) and report what you fixed.`);
  compile = await compileAll();
}
report(
  { item: "Compile", status: compile.pending ? "pending TeX install" : compile.ok ? "passed" : "failed", detail: compile.pending ? "MiKTeX not ready; main agent will compile after install" : compile.tail.slice(-200) },
  "progress",
);

phase("Four-gate review and revision");
const gateReviewer = agent("gate reviewer", {
  system: "You are the dissertation-style supervisor adapted to a journal manuscript. Judge against the rubric and the BRIEF; GO only when every gate passes. Escalate scope questions instead of reinterpreting them.",
});
const gateAsk = (attempt: string) =>
  `${BRIEF_LINE}${READ_ONLY}
Apply the four gates of Documentation/Reports and Papers/Dissertation/Notes/revision_campaign_20261015/supervisor-rubric.md,
adapted to this journal manuscript (${attempt}) — read ${MAINTEX} fully:
Gate 1 self-plagiarism/adaptation: proper adaptation of Ben's dissertation text into his own
unpublished manuscript (bolen_isometric_2025); no verbatim reuse of the published Actuators paper
(bolen_2026) beyond short quotation norms. Gate 2: slop regexes over the new prose. Gate 3: sample 10
load-bearing numbers and verify each against the ledger or chapter digest. Gate 4 directive intent:
no campaign dates, file names, or pick numbers; fictive-strain framing present; the long Discussion
carries the follow-up identification and route-redesign results; author block untouched. Verdict
GO or BLOCK with file:line findings.`;
let gate = await gateReviewer.ask<GateVerdict>(gateAsk("first pass"));
let gateRounds = 0;
const reviser = agent("reviser");
while (gate.verdict === "BLOCK" && gateRounds < 2) {
  gateRounds++;
  await reviser.ask(
    `${BRIEF_LINE}
Fix these gate findings in the manuscript (edit ONLY Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/
files), minimum churn per finding: ${JSON.stringify(gate.findings)}. Report what you changed.`,
  );
  if (!compile.pending) compile = await compileAll();
  gate = await gateReviewer.ask<GateVerdict>(gateAsk(`revision ${gateRounds}`));
}
report({ item: "Gates", status: gate.verdict === "GO" ? "GO" : "BLOCK", detail: `${gate.findings.length} findings after ${gateRounds} revision rounds` }, "progress");

const coldRead = await agent("cold reader").ask<ColdRead>(
  `${BRIEF_LINE}${READ_ONLY}
Read ${MAINTEX} the way a Frontiers reviewer would, from the text alone (checks are someone else's
job). Report: what is unclear; which claims lack support; what a reviewer will ask first.`,
);
if (coldRead.unsupported.length > 0) {
  await reviser.ask(
    `${BRIEF_LINE}
The cold reader flagged these unsupported or unclear points; fix the top ones with minimum churn
(edit only ProofFinal files), and report which you fixed vs deferred: ${JSON.stringify({
      unclear: coldRead.unclear,
      unsupported: coldRead.unsupported,
    })}`,
  );
  if (!compile.pending) compile = await compileAll();
}

await artifact.markdown(
  "bridge-report",
  [
    "# Isometric torque paper — bridge report",
    "",
    `- Manuscript: ${MAINTEX}`,
    `- Plan: tmp/paper_bridge/plan.md; digests in ${MINING}/`,
    `- Tables (xlsx): ${tables.xlsxFiles.join(", ") || "none"}`,
    `- Figures regenerated: ${figures.regenerated.join(", ") || "none"}; skipped (Ben decides): ${figures.skipped.join(", ") || "none"}`,
    `- Gates: ${gate.verdict} with ${gate.findings.length} findings after ${gateRounds} revision rounds`,
    `- Compile: ${compile.pending ? "PENDING (MiKTeX still installing)" : compile.ok ? "passed" : "failed"}`,
    "",
    "## Open questions for Ben",
    ...miners.flatMap((m) => m.openQuestions).map((q) => `- ${q}`),
    ...plan.notes.map((q) => `- (plan) ${q}`),
    "",
    "## Cold-read findings (reported, not all fixed)",
    ...coldRead.reviewerQuestions.map((q) => `- ${q}`),
  ].join("\n"),
  { title: "Paper bridge report", description: "What was bridged, verified, and what still needs Ben.", primary: true },
);
if (compile.ok) {
  try {
    await artifact.file("manuscript-pdf", "Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Main.pdf", {
      title: "Compiled manuscript (local MiKTeX)",
    });
  } catch {
    log("PDF artifact publish skipped (file missing or over cap).");
  }
}

interface ReportFinding {
  where: string;
  what: string;
  evidence: string;
  status: "verified" | "unconfirmed";
  severity: "low" | "medium" | "high";
}
interface WorkflowReport {
  conclusion: string;
  findings: ReportFinding[];
  verified: string[];
  notCovered: string[];
}
const result: WorkflowReport = {
  conclusion:
    `Bridged the committed dissertation torque chapter and archived test results into the Frontiers manuscript ` +
    `(${plan.sections.length} sections, ${tables.xlsxFiles.length} tables with xlsx, ${figures.regenerated.length} figures regenerated). ` +
    `Gates: ${gate.verdict}. Compile: ${compile.pending ? "pending MiKTeX install (main agent will compile)" : compile.ok ? "passed locally" : "FAILED — see findings"}.`,
  findings: [
    ...gate.findings.map<ReportFinding>((f) => ({
      where: f.where,
      what: `${f.gate}: ${f.what}`,
      evidence: "gate reviewer over Main.tex against supervisor-rubric + ledger",
      status: "unconfirmed",
      severity: f.severity,
    })),
    ...figures.skipped.map<ReportFinding>((f) => ({
      where: f,
      what: "Figure regeneration skipped: ambiguous data source or content; Ben decides.",
      evidence: "figures regenerator report",
      status: "unconfirmed",
      severity: "medium",
    })),
  ],
  verified: [
    compile.pending ? "compile pending: MiKTeX was still installing" : `local pdflatex+bibtex compile ${compile.ok ? "passed" : "ran with errors"}`,
    "gate 3 sampled 10 load-bearing numbers against the results ledger / chapter digest",
    "plan reviewed by an independent reviewer before drafting",
  ],
  notCovered: [
    "Overleaf upload of the bridged source (Ben's click)",
    "author-list discrepancy: Main.tex lists 4 authors vs 7 contributors in the dissertation publication note (Ben's call)",
    "Connor-comments PDF (Sep 2025 review round) not incorporated — not requested",
    "dissertation chapter continues to be revised until 2026-10-21; reconcile again before submission if it changes",
  ],
};
return result;
