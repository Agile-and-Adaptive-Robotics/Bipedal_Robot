"""Retrieve historical advisor examples for prose and verify recorded local judgments."""
import argparse
import calendar
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path


HERE = Path(__file__).resolve().parent
NOTES = HERE.parent
CHAPTERS = NOTES.parent / "Final" / "chapters"
DEFAULT_CHAPTERS = [CHAPTERS / "60-sns.tex", CHAPTERS / "70-overall.tex"]
RECENT_TAGS = {
    "context": ["reader_context", "gap_and_motivation"],
    "interface": ["mechanism", "definitions_and_terminology"],
    "intro_actuator_model": ["mechanism", "evidence_and_validation"],
    "intro_torque_mechanism": ["mechanism", "reader_context"],
    "intro_torque_prediction": ["evidence_and_validation", "mechanism"],
    "intro_route_design": ["result_first", "evidence_and_validation"],
    "intro_neural_motivation": ["gap_and_motivation", "reader_context"],
    "intro_sadb_contribution": ["contribution_and_credit", "reader_use"],
}
QUESTIONS = {
    "reader_context": "Can a reader identify the system, its parts, and the question before the details arrive?",
    "gap_and_motivation": "What could the previous method not tell the researcher?",
    "mechanism": "What acts on what, through which physical or neural quantity?",
    "evidence_and_validation": "What was changed, what was held fixed, what was measured, and what comparison answers the question?",
    "result_first": "What happened in the test, and what design decision follows?",
    "contribution_and_credit": "Who measured, built, or adapted it? What was already available?",
    "reader_use": "What can a researcher select, run, or measure with the existing tools? Which experiments would extend the present evidence?",
    "ownership_and_narration": "Does first person identify the author of a contribution or measurement? Can a tool or experiment state its own capability or result? Does the paragraph merely announce what the writer will describe?",
    "definitions_and_terminology": "At its first appearance in the printed chapter, can the reader say what a named model, controller, program, or signal is and what it does? For a variant, what changes, and which body does it control?",
    "narrative_structure": "Why does this paragraph follow the previous one?",
    "concision_and_repetition": "Does a repeated statement supply a new comparison or interpretation?",
    "scope_and_uncertainty": "Do the preparation, support, condition, and measurement justify the inference?",
}
STOP = set("a an and are as at be been by can did do does for from had has have how i in into is it its may of on one or our that the their these they this to was we were what when where which who will with would".split())
MONTHS = {name.lower(): i for i, name in enumerate(calendar.month_abbr) if i}


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def digest(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def ownership_cues(passage):
    """Locate phrases for a reader to inspect. A first-person pronoun alone carries no score."""
    patterns = {
        "roadmap_or_intention": r"\b(?:I (?:begin|show|describe|next describe|wanted|needed|would)|(?:show|describe|explain) how I)\b",
        "attributed_contribution": r"\b(?:I|we) (?:developed|built|adapted|implemented|connected|measured|verified|compared|tested|timed|added|assembled)\b",
        "named_capability": r"\b(?:the (?:new )?(?:system|editor|interface|script|tools?)|a researcher|researchers) can\b",
    }
    return {kind: found for kind, pattern in patterns.items()
            if (found := re.findall(pattern, passage, flags=re.IGNORECASE))}


def dated_signal(row, batch):
    if row.get("date"):
        return row["date"]
    day, month = row["time"].split(",", 1)[0].split()[:2]
    return f"{batch[:4]}-{MONTHS[month[:3].lower()]:02}-{int(day):02}"


def signals():
    history = read_json(HERE / "HISTORICAL_SIGNALS.json")
    known = {}
    for row in history["signals"]:
        if row["id"] in known:
            raise ValueError(f"Duplicate historical event {row['id']}")
        known[row["id"]] = row
    message_ids = {r["platform_id"]: r["id"] for r in known.values()
                   if r["type"] == "comment" and r.get("platform_id")}
    for path in sorted(NOTES.glob("advisor_scoring_*/CLASSIFIER_RESULTS.json")):
        doc = read_json(path)
        for row in doc["advisor_signals"]:
            if row["author"] != "ajh26":
                continue
            kind = "tracked_edit" if row.get("kind") in {"replace", "insert", "delete"} else "comment"
            edit = {k: row.get(k, "") for k in ("original", "proposed")} if kind == "tracked_edit" else None
            if edit is not None:
                edit["operation"] = row["kind"]
            entry = {
                "id": row["id"], "type": kind, "project": "dissertation", "author": row["author"],
                "date": dated_signal(row, path.parent.name.split("_")[-1]), "time_local": row["time"],
                "file": row["file"], "platform_id": row.get("comment_id"), "aliases": [],
                "comment": row.get("text"), "edit": edit, "score": row["original_score"],
                "tags": RECENT_TAGS.get(row["issue"], row.get("tags", [])),
                "rationale": row["reason"], "issue_id": row["issue"],
                "anchor": row.get("live_tex_excerpt") or row.get("anchor"),
                "sources": [path.relative_to(NOTES).as_posix() + "#" + row["id"]],
            }
            if not entry["tags"]:
                raise ValueError(f"Assign historical presentation tags to {entry['id']} before using this capture")
            platform_id = entry.get("platform_id")
            if kind == "comment" and platform_id:
                alias = message_ids.setdefault(platform_id, entry["id"])
                if alias != entry["id"]:
                    raise ValueError(f"Reconcile two IDs for the same Overleaf message: {alias}, {entry['id']}")
            prior = known.get(entry["id"])
            if prior:
                # The dated ledgers may repeat an earlier capture. An inconsistent
                # identity must be reconciled from the source before being scored.
                for key in ("author", "file", "date", "type", "score", "comment", "edit"):
                    if prior[key] != entry[key]:
                        raise ValueError(f"Inconsistent recapture of {entry['id']}, field {key}")
            else:
                known[entry["id"]] = entry
    return list(known.values())


def words(text):
    text = re.sub(r"\\(?:ref|label|cite\w*)\s*(?:\[[^\]]*\])?\{[^}]*\}", " ", text)
    text = re.sub(r"\\[a-zA-Z]+", " ", text.casefold())
    return [w for w in re.findall(r"[a-z][a-z0-9]+", text) if w not in STOP]


def features(text):
    tokens = words(text)
    return Counter(tokens + [a + " " + b for a, b in zip(tokens, tokens[1:])])


def as_text(value):
    if isinstance(value, str):
        return value
    if isinstance(value, dict):
        return " ".join(as_text(v) for v in value.values())
    if isinstance(value, list):
        return " ".join(as_text(v) for v in value)
    return ""


def refers_to_paragraph(signal, paragraph):
    """An exact text/context match reopens local review. A shared topic or line number alone does not."""
    anchor = signal.get("anchor")
    if isinstance(anchor, dict):
        candidates = [anchor.get("text"), anchor.get("context")]
    else:
        candidates = [anchor]
    paragraph = " ".join(paragraph.split())
    for candidate in candidates:
        if not isinstance(candidate, str):
            continue
        selected = " ".join(candidate.split())
        if len(selected) >= 24 and (selected in paragraph or paragraph in selected):
            return True
    return False


def recurring_questions(rows):
    """Give the semantic reviewer direct examples across dates, including distant topics."""
    result = {}
    for tag, question in QUESTIONS.items():
        relevant = [r for r in rows if r["score"] < 0 and tag in r["tags"]]
        if not relevant:
            continue
        dates = sorted({r["date"] for r in relevant})
        examples = []
        for month in sorted({date[:7] for date in dates}):
            choices = [r for r in relevant if r["date"].startswith(month)]
            # Prefer an actual specific reader question over a punctuation edit.
            chosen = min(choices, key=lambda r: (not bool(r.get("comment")),r["score"],
                         -min(len(r.get("comment") or ""),240),r["id"]))
            examples.append({"id":chosen["id"],"date":chosen["date"],"file":chosen["file"],
                             "comment":chosen.get("comment"),"edit":chosen.get("edit"),"issue":chosen["issue_id"]})
        result[tag] = {"question":question,"events":len(relevant),"event_dates":dates,
                       "distinct_issues":len({r["issue_id"] or r["id"] for r in relevant}),
                       "historical_examples_by_month":examples}
    return result


class EvidenceIndex:
    """Cosine retrieval over distinct historical issues. Examples never assign a grade."""
    def __init__(self, rows):
        self.rows = [r for r in rows if r["score"] < 0]
        counts = []
        frequency = Counter()
        for row in self.rows:
            content = " ".join(as_text(row.get(k)) for k in ("anchor", "comment", "edit", "rationale"))
            terms = features(content)
            frequency.update(terms.keys())
            counts.append(terms)
        self.idf = {key: 1 + math.log((1 + len(counts)) / (1 + count)) for key, count in frequency.items()}
        self.vectors = [self.vector(terms) for terms in counts]

    def vector(self, counts):
        weights = {k: (1 + math.log(v)) * self.idf[k] for k, v in counts.items() if k in self.idf}
        length = math.sqrt(sum(v*v for v in weights.values()))
        return {k: v/length for k, v in weights.items()}

    def nearest(self, paragraph, limit=3):
        query = self.vector(features(paragraph))
        ranks = sorted(
            ((sum(weight * vector.get(key, 0) for key, weight in query.items()), row)
             for row, vector in zip(self.rows, self.vectors)),
            key=lambda pair: (-pair[0], pair[1]["id"]),
        )
        matches, issues = [], set()
        for similarity, row in ranks:
            if similarity <= 0 or len(matches) == limit:
                break
            identity = row["issue_id"] or row["id"]
            if identity in issues:
                continue
            issues.add(identity)
            matches.append({"id": row["id"], "date": row["date"], "issue": identity, "tags": row["tags"],
                            "cosine": round(similarity, 4)})
        return matches


def passages(path):
    """Blank-delimited narrative plus distinct caption/table lines, retaining TeX anchors."""
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    start, parts = None, []
    in_table = False
    environments_to_skip = {"align", "align*", "equation", "equation*", "verbatim", "lstlisting"}
    skip = None
    for number, line in enumerate(lines + [""], 1):
        s = line.strip()
        if skip:
            if r"\end{" + skip + "}" in s:
                skip = None
            continue
        opening = re.match(r"\\begin\{([^}]+)\}", s)
        if opening and opening[1] in environments_to_skip:
            skip = opening[1]
        cell = in_table and "&" in s
        caption = s.startswith(r"\caption")
        wrapper = bool(re.match(r"^\\noindent\\begin\{minipage\}\{[^}]+\}", s)) and len(words(s)) > 4
        source_prose = bool(re.match(r"^[A-Za-z({$\\]", s)) and not s.startswith((r"\toprule", r"\midrule", r"\bottomrule"))
        structural = s.startswith("\\") and not (caption or cell or wrapper or re.match(r"^\\(?:cite[tp]|emph|textit|textbf|qty)[{[]", s))
        narrative = source_prose and not structural and not in_table and not skip and not caption
        if not narrative and parts:
            yield start, number - 1, "prose", "\n".join(parts)
            start, parts = None, []
        if caption or cell:
            yield number, number, "caption" if caption else "table", line
        elif narrative:
            start = start or number
            parts.append(line)
        if s.startswith((r"\begin{tabular", r"\begin{longtable")):
            in_table = True
        if s.startswith((r"\end{tabular", r"\end{longtable")):
            in_table = False


def make_review(paths, rows, decisions=None):
    index = EvidenceIndex(rows)
    saved = defaultdict(list)
    for case in (decisions or {}).get("passages", []):
        saved[case["file"]].append(case)
    prior_ids = set((decisions or {}).get("advisor_event_ids_at_review", ()))
    unseen = [r for r in rows if prior_ids and r["id"] not in prior_ids and r["score"] < 0]
    outputs, found, reopened, newly_anchored = [], set(), [], set()
    for path in paths:
        for start, end, kind, content in passages(path):
            sha = digest(content)
            case = next((c for c in saved[path.name] if sha == c["after_sha256"] and c["id"] not in found), None)
            output = {"file": path.name, "start": start, "end": end, "kind": kind, "source": content,
                      "sha256": sha, "historical_examples": index.nearest(content),
                      "status": "unlabeled", "local_judgment": None}
            if cues := ownership_cues(content):
                output["ownership_phrases_for_contextual_reading"] = cues
            if case:
                found.add(case["id"])
                output.update(status="locally reviewed", local_judgment=case["id"])
                later = [r["id"] for r in unseen if r["file"] == path.name and refers_to_paragraph(r, content)]
                if later:
                    reopened.append(case["id"])
                    newly_anchored.update(later)
                    output.update(status="reopened by previously unseen advisor feedback",
                                  new_advisor_signals=later, local_judgment=None, prior_judgment=case["id"])
            outputs.append(output)
    assessed_files = {p.name for p in paths}
    stale = [c["id"] for c in (decisions or {}).get("passages", [])
             if c["file"] in assessed_files and c["id"] not in found]
    months = defaultdict(list)
    for row in rows:
        months[row["date"][:7]].append(row["score"])
    return {
        "method": "Sublinear TF-IDF with word unigrams/bigrams and cosine similarity ranks prior advisor examples. A language-model reader makes the separate semantic judgment in the dated passage record. Similarity, silence, and a local score of 3 do not establish advisor acceptance.",
        "questions_for_every_pass": QUESTIONS,
        "evidence_by_recurring_reader_question":recurring_questions(rows),
        "teaching_examples_from_prior_revisions":read_json(HERE/"EDITORIAL_EXAMPLES.json")["examples"],
        "advisor_signal_totals_by_event_month":
            {k: {"signals": len(v), "score": sum(v)} for k, v in sorted(months.items())},
        "observed_advisor_signals": len(rows), "observed_advisor_score": sum(r["score"] for r in rows),
        "saved_judgments_stale_or_absent": stale, "reopened_judgments":reopened,
        "additional_same_chapter_feedback_to_check":[r["id"] for r in unseen
            if r["file"] in assessed_files and r["id"] not in newly_anchored], "passages": outputs,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("chapters", nargs="*", type=Path, help="TeX files, defaults to local Chapters 6 and 7")
    parser.add_argument("--before", help="Use only historical signals strictly earlier than this ISO event date")
    parser.add_argument("--decisions", type=Path, help="Verify a dated semantic review against the current text")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = [r for r in signals() if not args.before or r["date"] < args.before]
    if not rows:
        parser.error("No earlier advisor examples in the selected period")
    result = make_review(args.chapters or DEFAULT_CHAPTERS, rows, read_json(args.decisions) if args.decisions else None)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"{len(rows)} advisor signals; {len(result['passages'])} passage/caption/table units; "
          f"{len(result['saved_judgments_stale_or_absent'])} saved judgments missing or stale; "
          f"{len(result['reopened_judgments'])} reopened.")
    if result["saved_judgments_stale_or_absent"] or result["reopened_judgments"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
