"""Generate MODEL_XREF.md — the SADb → neuromechanical-modeling cross-reference.

RULE FOR SESSIONS: before building/tuning any neuromechanical model in
AnimatLab, SNS Simulink/SNS_Simscape, or the MuJoCo-SNS toolbox, READ THIS
FILE (regenerate after curation batches: myo python build_model_xref.py).

Contents are GENERATED from the corpus (feedback rules + key models + key
studies stay current automatically); the platform cross-reference section is
hand-authored from facts already recorded in AGENTS.md.
"""
import json, os
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
EXPORT = os.path.join(HERE, "export")
records = json.load(open(os.path.join(EXPORT, "sadb_export.json"), encoding="utf-8"))
layout = json.load(open(os.path.join(EXPORT, "sadb_layout.json"), encoding="utf-8"))

# ---- feedback rule index: rule -> papers demonstrating it (live first) ----
rule_papers = defaultdict(list)
for r in records:
    for f in r["feedback"]:
        rule_papers[f].append(r)

# ---- key models: papers that ARE models (Is-the-model-paper links) ----
model_papers = defaultdict(list)
for r in records:
    for m in (r.get("models2") or []):
        model_papers[m].append(r)

top = sorted(records, key=lambda r: -(layout.get(r["id"], {}).get("c", 0) or 0))

out = ["# MODEL_XREF — SADb rules/studies/models for neuromechanical modeling",
       "",
       "**Standing rule (Ben, 2026-10-02):** any session building, tuning, or",
       "debugging a neuromechanical model in **AnimatLab**, **SNS Simulink /",
       "SNS_Simscape**, or the **MuJoCo-SNS toolbox** — on any machine — reads",
       "this file FIRST and grounds pathway wiring/gains in the cited studies.",
       "Regenerate after every curation batch:",
       "`myo python SADb_audit/build_model_xref.py`. Deeper browsing:",
       "`SADb_audit/knowledge_base/INDEX.md` and the SADb Explorer app",
       "(`SADb_audit/app/sadb_app.html`, 943 papers incl. archived).",
       "",
       f"Corpus snapshot: {len(records)} papers "
       f"({sum(1 for r in records if r.get('archived'))} app-only/archived).",
       "",
       "## 1. Feedback-pathway RULES (the wiring vocabulary)",
       "",
       "Each rule below is demonstrated by the linked papers (author year,",
       "OpenAlex citations; [A] = archived from Airtable, still in the app).",
       "When you implement a pathway, cite the demonstrating paper in comments.",
       ""]
for rule in sorted(rule_papers, key=lambda k: -len(rule_papers[k])):
    papers = sorted(rule_papers[rule], key=lambda r: -(layout.get(r["id"], {}).get("c", 0) or 0))
    refs = ", ".join(f"{r['primary'] or '?'} {r['year'] or ''}"
                     f"{'[A]' if r.get('archived') else ''}"
                     for r in papers[:4])
    out.append(f"- **{rule}** — {len(papers)} paper(s): {refs}")

out += ["",
        "## 2. Key MODELS (papers that ARE models — 'Is the model paper' links)",
        "",
        "NOTE: this field has known multi-link pollution — treat names as",
        "leads, verify against the paper before citing.",
        ""]
for m in sorted(model_papers, key=lambda k: -len(model_papers[k]))[:25]:
    ps = model_papers[m]
    if len(ps) == 1 and (ps[0]["primary"] or "").lower() in m.lower():
        r = ps[0]
        out.append(f"- **{m}** — {r['title'][:80]} ({r['year'] or '?'}) · "
                   f"{layout.get(r['id'], {}).get('c', 0)} cites"
                   f"{' · [A]' if r.get('archived') else ''}")

out += ["", "## 3. Key STUDIES (top-cited; cite these for mechanisms)", ""]
for r in top[:30]:
    c = layout.get(r["id"], {}).get("c", 0) or 0
    note = (r.get("notes") or "").strip().split(". ")[0][:160]
    out.append(f"- **{r['primary'] or '?'} {r['year'] or ''}** ({c} cites"
               f"{' · [A]' if r.get('archived') else ''}) — {r['title'][:85]}"
               + (f"  \n  {note}." if note else ""))

out += ["",
        "## 4. Afferent cheat sheet (encoding conventions used across the stacks)",
        "",
        "- **Ia** (spindle primary): length + velocity. MuJoCo-SNS encoder:",
        "  `(L−Lmid)/Lhalf` and `L̇/0.6` → afferent current (spinal/ build_network).",
        "- **Ib** (Golgi tendon organ): force. Encoder: `F/Fmax`.",
        "- **II** (spindle secondary): static length; often stance-gated to avoid",
        "  a global co-contraction floor (the ungated-i0 lesson in DESIGN.md).",
        "- **group III/IV**: metabolic/fatigue-nociceptive — 'group III/IV fatigue' rule.",
        "- **Cutaneous / mechanosensory**: heel/toe contact sensors → HEEL_IN/TOE_IN",
        "  ports; insect hair plates/campaniform/chordotonal = Mechanosensory.",
        "- **Phase → afferent gain** (Akazawa 1982 direction): reflex gain is",
        "  locomotor-phase dependent — the RG/PF state gates every afferent pathway.",
        "",
        "## 5. Platform cross-reference (where rules live in each toolchain)",
        "",
        "| Rule family | MuJoCo-SNS (spinal/) | SNS Simulink (SNS_Simscape) | AnimatLab (Neuromechanical_Models) |",
        "|---|---|---|---|",
        "| Rhythm generation (RG half-centers, Brown 1911 lineage) | per-leg RG; NaP/tau-h traps in AGENTS | SNS_Library NonSpikingNeuron; SNS_SpinalNetwork | Biped_2xCPG_wSubs RG (contact-driven; W2L) |",
        "| Ia reciprocal inhibition | params ia_in / IaIN population | KneeReflexDemo reciprocal Ia pair | Deng-style Ia-IN chains in aproj |",
        "| Ib autogenic + stance-gated reversal (IBEXC) | full_rules Ib branch; ib_rge | — | W2L Ib autogenic excitation |",
        "| Ib/Ia disynaptic excitation & inhibition (Angel/Jankowska) | full_rules group-I INs | — | Ia relay INs in W2L |",
        "| Heel/toe contact reset at PF layer | heel_pf_layer, toe_df_inh (Ben's rules JSON) | — | contact neurons → RG (W2L reference) |",
        "| Renshaw recurrent inhibition | --renshaw G | — | Renshaw chains in aproj |",
        "| Stance-gated phase reset / PRESET transients | PRESET_E/F high-pass onset | — | — |",
        "| II afferents stance-gated | ii gating in build_network | — | — |",
        "| Knee convention traps | flexion-negative knee (audit_signs) | — | — |",
        "",
        "Platform facts above are summaries of AGENTS.md's toolchain sections —",
        "read the full sections there before implementation work.",
        "",
        "## 6. Non-negotiables when porting rules across platforms",
        "",
        "- Tau-h semantics differ (sns_toolbox tau_h(V) quenches NaP half-centers;",
        "  fixed tau_h works — spinal/_tau_h_check.py, deng_cpg_ode.py).",
        "- Ia/Ib/II afferents must be pure-signal (resting tone → constant",
        "  reciprocal inhibition — the i0 lesson).",
        "- Summation-order chaos flips marginal limit cycles between backends",
        "  (basin gate: basin_gate.py).",
        "- SNS_Library vs sns_toolbox use OPPOSITE synapse-saturation conventions.",
        ""]

path = os.path.join(HERE, "MODEL_XREF.md")
open(path, "w", encoding="utf-8").write("\n".join(out) + "\n")
print(f"wrote {path}: {len(rule_papers)} rules, {len(model_papers)} model names, "
      f"top-{30} studies, afferent sheet, platform xref")
