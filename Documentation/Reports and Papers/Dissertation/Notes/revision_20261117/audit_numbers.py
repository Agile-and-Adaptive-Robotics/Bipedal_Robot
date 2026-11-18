"""Compare every numeric-bearing edited TeX block against the independently captured live source."""
import argparse
import difflib
import json
import re
from pathlib import Path

notes=Path(__file__).resolve().parent
dis=notes.parents[1]
parser=argparse.ArgumentParser()
parser.add_argument("baseline",type=Path,help="Directory containing the 23 individual captured live source files")
args=parser.parse_args()
active=dis/"Final/chapters"
number=re.compile(r"[-+−]?(?:\d+(?:,\d{3})*(?:\.\d+)?|\.\d+)(?:[eE][-+]?\d+)?")
def study(ch,n):
    if ch in ("92","93","94","96"):
        return {"92":"A.1–A.4, signed torque exports, bracket fits, and cable implementation",
                "93":"B.1–B.3, spiking and graded controller recordings and input ports",
                "94":"C.1–C.3, saved route and fold metadata, polar arc, and forward pose grid",
                "96":"E.1–E.3, SADb export and primary experiments"}[ch], "APPENDIX_AC_E_EDITS.md"
    if ch=="70":
        if n<40: return "Published resting-length work and proposed BPA pressure series or boundary conditions", "FIGURE_PROVENANCE.json; REVISION_REPORT.md"
        if n<51: return "Joint mechanical power proposal and author ruling on a fifth factor", "REVISION_REPORT.md"
        if n<62: return "SADb Table view and McCrea experimental preparations","INDEPENDENT_REVIEW.md"
        return "Supported human and W2L experiments or proposed follow-up","DEFENSE_EVIDENCE.md; APPENDIX_D_EDITS.md"
    if ch=="95": return "Section-specific configuration and precisely defined scoring intervals by revised line","APPENDIX_D_EDITS.md; DEFENSE_SOURCE_CHECKS.json; TECHNICAL_READER.md"
    if n<=148: return "Course project, tools, model identities, and experiment definitions","TECHNICAL_READER.md; DEFENSE_EVIDENCE.md"
    if n<=165: return "Matched poses, bridge, circuit and supported W2L rhythm-core tests","APPENDIX_D_EDITS.md; TECHNICAL_READER.md"
    if n<=177: return "Knee pulse from t=30 s or cup ramp t=10 s, 0.625 kg","APPENDIX_D_EDITS.md; TECHNICAL_READER.md"
    if n<=202: return "BiW2Lvar+ 16 s air runs, cadence at 5 ≤ t < 14 s; 5–16 s spring support; phase-locked stimuli","DEFENSE_SOURCE_CHECKS.json; APPENDIX_D_EDITS.md"
    if n<=214: return "Historical W2L intervention, 30 s protocol with 8–8.15 s common disturbance; syn6 timing checkpoints","APPENDIX_D_EDITS.md; TECHNICAL_READER.md"
    if n<=221: return "Original and repaired syn6 PF checks, offline or in-body interval per trace","DEFENSE_SOURCE_CHECKS.json; APPENDIX_D_EDITS.md"
    if n<=237: return "W2L108 and CAD supported runs, printed post-startup or whole-trial windows","DEFENSE_SOURCE_CHECKS.json; APPENDIX_D_EDITS.md"
    if n<=258: return "Timing workload, clock and parity for each named neural backend","DEFENSE_EVIDENCE.md; APPENDIX_D_EDITS.md"
    return "Coupled controller discussion and supported drive-weaning conditions","TECHNICAL_READER.md; APPENDIX_D_EDITS.md"

rows=[]
for prefix in ["60","70","92","93","94","95","96"]:
    end=next(active.glob(prefix+"-*.tex"))
    old=(args.baseline/end.name).read_text(encoding="utf-8-sig").splitlines()
    new=end.read_text(encoding="utf-8-sig").splitlines()
    for tag,i,j,k,l in difflib.SequenceMatcher(None,old,new,autojunk=False).get_opcodes():
        if tag=="equal": continue
        a,b="\n".join(old[i:j]),"\n".join(new[k:l])
        oa,nb=number.findall(a.replace("--"," – ")),number.findall(b.replace("--"," – "))
        if not (oa or nb): continue
        config,proof=study(prefix,i+1)
        rows.append({"file":end.name,"original_lines":[i+1,j],"revised_lines":[k+1,l],"operation":tag,
                     "original_exact":a,"revised_exact":b,"old_numeric_tokens":oa,"new_numeric_tokens":nb,
                     "numeric_token_order_or_values_changed":oa!=nb,"context":config,"evidence":proof})
out={"baseline":"Live Overleaf text, 2026-11-18 00:29:14 UTC, independently matched against local project ZIP.",
     "method":"Every nonidentical line block containing a digit. Includes numeric labels, file names, version identifiers, math subscripts, formatting, and values as printed. Empty side denotes insertion or deletion. Full wording and interval definitions disambiguate reordered values. See the human ledger for scientific meaning.","rows":rows}
(notes/"NUMERICAL_SOURCE_DIFF.json").write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
for row in rows:
    if row["file"][:2] in ("60","70") and row["numeric_token_order_or_values_changed"]:
        print(row["file"][:2],row["original_lines"]," ".join(row["old_numeric_tokens"])[:95]," -> "," ".join(row["new_numeric_tokens"])[:110])
