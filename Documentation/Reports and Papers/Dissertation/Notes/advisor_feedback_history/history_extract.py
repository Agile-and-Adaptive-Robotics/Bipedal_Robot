"""Rebuild the dated historical dissertation-advisor evidence from immutable captures.

Run from any directory within this repository with Python 3.10+. No packages are
required. The EB475WS4 captures are read-only. An existing historical JSON file
or the original hash-checked lab files can reproduce those archived inputs.
"""
import argparse
import base64
from collections import Counter, defaultdict
from datetime import datetime
from hashlib import sha256
from html.parser import HTMLParser
from io import BytesIO
import json
from pathlib import Path
import re
import subprocess
import unicodedata
from zipfile import ZipFile
import xml.etree.ElementTree as ET
import zlib

HERE = Path(__file__).resolve().parent
ROOT = Path(
    subprocess.check_output(["git", "-C", str(HERE), "rev-parse", "--show-toplevel"], text=True).strip()
)
DIS = "Documentation/Reports and Papers/Dissertation/"
MIRROR = "a35f0e36f268be3fc324de2da8b4cc6b433e7397"
EARLY = "1d0cfc2b1f98193a5ab72e4f2d3caaec5d3c5288"
JULSEP = "5ba9fff684e4565fb63758fe0532f6bcb71c30b6"
OCT26 = "d192bc766fb71a3543efaac58f67b73f924fa17e"
PROJECT = "644dc82bcd6a5481e3c3a8bb"
DOM = DIS + "Notes/advisor_review_20261023/"
EARLIER = "Notes/advisor_review_20261022/"
DAY = re.compile(r"(\d{1,2} [A-Za-z]+)(?:, (\d{1,2}:\d{2} [ap]m))?", re.I)
TAGS = {
    "context": (
        "reader_context",
        "Introduce the system, condition, or antecedent needed to follow the claim.",
    ),
    "gap": ("gap_and_motivation", "Name the missing relationship or evidence and why it motivated the work."),
    "mech": ("mechanism", "Explain the physical, neural, or computational operation and its causal link."),
    "test": ("evidence_and_validation", "Specify the test, measured observable, comparison, and result."),
    "result": ("result_first", "Lead with the scientific result so that the figure or detail has a purpose."),
    "credit": ("contribution_and_credit", "Identify who developed, adapted, supplied, or tested the work."),
    "use": ("reader_use", "Show what the reader can do and how an artifact or example supports that action."),
    "term": (
        "definitions_and_terminology",
        "Define the term, symbol, or referenced quantity in this context.",
    ),
    "short": ("concision_and_repetition", "Condense a repetition, excessive detail, or a low-value aside."),
    "flow": (
        "narrative_structure",
        "Put material in the section and sequence that support the scientific argument.",
    ),
    "fig": (
        "visual_communication",
        "Make a figure, table, caption, cross-reference, or page technically readable.",
    ),
    "unit": ("notation_and_units", "Correct or explain notation, signs, formatting of quantities, or units."),
    "num": (
        "numerical_accuracy",
        "Correct a number, metric, input, sign, category, or mathematical interpretation.",
    ),
    "cite": (
        "sources_and_citations",
        "Identify the primary source or ensure that a claim is supported by its citation.",
    ),
    "copy": ("copyediting", "Repair wording, syntax, spelling, or punctuation."),
    "scope": ("scope_and_uncertainty", "Bound the claim by what the experiment or data actually support."),
    "praise": (
        "praise",
        "Positive assessment without a request for change, or a praise component of a mixed message.",
    ),
    "admin": (
        "editorial_administration",
        "A review action or presentation requirement without a technical mechanism.",
    ),
}
SOURCES = {}
RAW = {}
EVENTS = []
GAPS = []
EXCLUDED = Counter()


def git(*args):
    return subprocess.check_output(["git", "-C", str(ROOT), *args])


def digest(v):
    return sha256(v).hexdigest()


def source(sid, path, commit=MIRROR):
    if sid not in RAW:
        spec = f"{commit}:{path}"
        data = git("show", spec)
        SOURCES[sid] = dict(
            id=sid,
            path=path,
            commit=commit,
            git_blob=git("rev-parse", spec).decode().strip(),
            sha256=digest(data),
            origin="Git",
            immutable_ref=spec,
        )
        RAW[sid] = data
    return RAW[sid]


def txt(sid, path, commit=MIRROR):
    return source(sid, path, commit).decode("utf-8-sig")


def jfile(sid, path, commit=MIRROR):
    return json.loads(txt(sid, path, commit))


def norm(s):
    return " ".join(
        unicodedata.normalize("NFC", str(s)).replace("Â\u00a0", " ").replace("\u00a0", " ").split()
    ).casefold()


def pdate(d):
    m = DAY.search(d)
    if not m:
        raise ValueError("unrecognized local date " + d)
    return datetime.strptime(m[1] + " 2026", "%d %B %Y").date().isoformat(), (m[2] or "").lower() or None


def anchor(status="unavailable", text=None, context=None, line=None, offset=None):
    return dict(status=status, text=text, context=context, line=line, offset=offset)


def add(
    sid,
    loc,
    uid,
    file,
    author,
    rawdate,
    kind,
    body=None,
    edit=None,
    anch=None,
    alias=None,
    pid=None,
    thread=None,
    quality="captured primary review",
    ordno=0,
):
    if author != "ajh26":
        EXCLUDED[(sid, author, kind)] += 1
        return
    date, time = pdate(rawdate)
    ref = dict(source_id=sid, locator=loc)
    if anch:
        anch["source_id"] = sid
    e = dict(
        file=file,
        author=author,
        date=date,
        time_local=time,
        type=kind,
        comment=body,
        edit=edit,
        anchor=anch or anchor(),
        aliases=[alias] if alias else [],
        platform_id=pid,
        thread_id=thread,
        evidence_quality=quality,
        sources=[ref],
        event_uid=uid,
        occurrence=ordno,
    )
    EVENTS.append(e)


def ce(op, a, b):
    return {"operation": op, "original": a, "proposed": b}


class Node:
    def __init__(self, tag, attrs):
        self.tag, self.attrs, self.children = tag, dict(attrs), []

    def text(self):
        return "".join(n if isinstance(n, str) else n.text() for n in self.children)

    def walk(self):
        for n in self.children:
            if not isinstance(n, str):
                yield n
                yield from n.walk()

    def cls(self, c):
        return [n for n in self.walk() if c in n.attrs.get("class", "").split()]


class DOMParser(HTMLParser):
    def __init__(self, markup):
        super().__init__(convert_charrefs=True)
        self.root = Node("", [])
        self.stack = [self.root]
        self.feed(markup)

    def handle_starttag(self, tag, attrs):
        n = Node(tag, attrs)
        self.stack[-1].children.append(n)
        if tag == "br":
            n.children.append("\n")
        if tag not in {"br", "hr", "img", "input", "path", "wbr"}:
            self.stack.append(n)

    def handle_endtag(self, tag):
        for i in range(len(self.stack) - 1, 0, -1):
            if self.stack[i].tag == tag:
                del self.stack[i:]
                return

    def handle_data(self, v):
        self.stack[-1].children.append(v)


def add23():
    matrices = jfile("23_matrix", DOM + "review_matrix.json")
    mm = defaultdict(list)
    for mr in matrices:
        mm[(mr["file"], mr["offset_utf16"], mr["kind"].endswith("tracked edit"))].append(mr)
    for fname in [
        "entries_final_20261023",
        "entries_20261023_1141",
        "entries_20261023_1111",
        "entries_manual_20261023_1110",
        "entries_20261023_1011",
        "entries_live",
        "entries_latest",
        "matrix_entries",
    ]:
        sid = "23_" + fname
        rows = jfile(sid, DOM + fname + ".json")
        for i, r in enumerate(rows):
            tree = DOMParser(r["html"]).root
            base = f"/{i}/html"
            f = r["file"]
            choices = mm.get((f, int(r["pos"]), r["type"] == "edit"), [])
            when = tree.cls("review-panel-entry-time")[0].text().strip()
            candidates = [m for m in choices if m["date"] == when]
            if fname == "entries_final_20261023":
                m = matrices[i]
                assert m in candidates, (i, f, when)
            else:
                m = candidates[0] if len(candidates) == 1 else None
            a = anchor("source context at captured offset", offset=r["pos"])
            alias = None
            if m:
                alias = m["id"]
                a = anchor(m["anchor_status"], context=m["context"], line=m["line"], offset=m["offset_utf16"])
                if (m.get("highlight") or "").startswith("Historical exact highlight"):
                    a["text"] = m["highlight"].split(": ", 1)[-1]
                    a["status"] = "historical exact selection, current source context recorded separately"
                if m.get("historical"):
                    a["historical_reference"] = m["historical"]
            if r["type"] == "comment":
                thread = None
                for k, n in enumerate(tree.cls("review-panel-comment")):
                    who = n.cls("review-panel-entry-user")[0].text().strip()
                    time = n.cls("review-panel-entry-time")[0].text().strip()
                    content = n.cls("review-panel-comment-body")[0].text()
                    ids = [
                        x.attrs["id"].split("-")[-1]
                        for x in n.walk()
                        if x.attrs.get("id", "").startswith("review-panel-comment-options-btn-")
                    ]
                    assert len(ids) == 1, (sid, i, k)
                    if k == 0:
                        thread = ids[0]
                    if "show more" in n.text().lower():
                        continue
                    add(
                        sid,
                        base + f"::comment[{k}]",
                        ids[0],
                        f,
                        who,
                        time,
                        "comment",
                        body=content,
                        anch=a.copy(),
                        alias=alias,
                        pid=ids[0],
                        thread=thread,
                    )
            else:
                node = tree.cls("review-panel-change-body")[0]
                who = tree.cls("review-panel-entry-user")[0].text().strip()
                time = tree.cls("review-panel-entry-time")[0].text().strip()
                if "show more" in node.text().lower():
                    continue
                vals = [x.text() for x in node.cls("review-panel-expandable-content")]
                icon = node.cls("material-symbols")[0].text().strip()
                if icon == "edit" and len(vals) == 2:
                    ed = ce("replace", *vals)
                elif icon == "add_circle" and len(vals) == 1:
                    ed = ce("insert", "", vals[0])
                elif icon == "delete" and len(vals) == 1:
                    ed = ce("delete", vals[0], "")
                else:
                    raise ValueError((sid, i, icon, vals))
                add(sid, base, None, f, who, time, "tracked_edit", edit=ed, anch=a.copy(), alias=alias)


def plain_edit(s):
    s = re.sub(r"\nshow less(?=(?:\n|$))", "", s)
    if "show more" in s:
        return None
    if "Changed:" in s:
        v = s.split("Changed:", 1)[1]
        v = v[1:] if v.startswith((" ", "\u00a0")) else v
        vals = re.split(r"\nto[ \u00a0]", v, maxsplit=1)
        return ce("replace", *vals) if len(vals) == 2 else None
    for head, op in [("Added:", "insert"), ("Deleted:", "delete")]:
        if head in s:
            v = s.split(head, 1)[1]
            v = v[1:] if v.startswith((" ", "\u00a0")) else v
            return ce(op, "", v) if op == "insert" else ce(op, v, "")
    return None


def add22():
    hl = {}
    for f in ["highlights", "highlights_final_additions"]:
        for row in jfile("22_" + f, EARLIER + f + ".json"):
            hl[(row["file"], row["position"])] = row.get("highlight")
    for f in ["advisor_inventory", "review_anchored", "review_final", "review_snapshot"]:
        sid = "22_" + f
        for ri, group in enumerate(jfile(sid, EARLIER + f + ".json")):
            offsets = Counter(x["position"] for x in group["entries"])
            for ei, entry in enumerate(group["entries"]):
                path = f"/{ri}/entries/{ei}"
                fn = group["file"]
                pos = entry["position"]
                h = entry.get("highlight") or (hl.get((fn, pos)) if offsets[pos] == 1 else None)
                a = anchor(
                    "exact historical selection" if h else "offset or nearby recorded context",
                    h,
                    entry.get("context"),
                    entry.get("line"),
                    pos,
                )
                if entry["type"] == "comment":
                    messages = entry["messages"]
                    tid = messages[0]["id"].split("-")[-1] if messages else None
                    for k, m in enumerate(messages):
                        pid = m["id"].split("-")[-1] if m.get("id") else None
                        if m["text"].rstrip().endswith("...") and f != "advisor_inventory":
                            continue
                        add(
                            sid,
                            path + f"/messages/{k}",
                            pid,
                            fn,
                            m["author"],
                            m["time"],
                            "comment",
                            body=m["text"],
                            anch=a.copy(),
                            alias=entry.get("report_id"),
                            pid=pid,
                            thread=tid,
                        )
                else:
                    e = plain_edit(entry.get("change", ""))
                    if e:
                        add(
                            sid,
                            path + "/change",
                            None,
                            fn,
                            entry["author"],
                            entry["time"],
                            "tracked_edit",
                            edit=e,
                            anch=a,
                            alias=entry.get("report_id"),
                        )


def fenced(block):
    return re.findall(r"^" + chr(96) * 3 + r"[^\n]*\n(.*?)\n" + chr(96) * 3 + r"(?=\n|$)", block, re.M | re.S)


def add15():
    sid = "15_report"
    s = txt(sid, "CHATGPT_REPORT.md", EARLY)
    matches = list(re.finditer(r"^#{3,4} ([CRT]\d+) [—-] (.+)$", s, re.M))
    current = None
    for i, m in enumerate(matches):
        tail = matches[i + 1].start() if i + 1 < len(matches) else len(s)
        block = s[m.end() : tail]
        alias, part = m.groups()
        fpart = re.findall(r"^### ([\w-]+\.tex)$", s[: m.start()], re.M)
        fn = part.split(" — ")[0] if alias[0] in "TR" else (fpart[-1] if fpart else None)
        fline = re.search(r"\*\*Source line:\*\* (\d+)", block)
        high = block.split("**Comment", 1)[0].split("**Author", 1)[0]
        h = fenced(high)
        a = anchor(
            "captured highlight" if h else "historical line or unresolved selection",
            h[0] if h else None,
            line=int(fline[1]) if fline else None,
        )
        if alias.startswith("T"):
            au = re.search(r"\*\*Author:\*\* (\w+) — ([^\n]+)", block)
            fs = fenced(block)
            op = part.split(" — ")[-1]
            ed = (
                ce("replace", *fs[:2])
                if op == "Changed"
                else (ce("insert", "", fs[0]) if op == "Added" else ce("delete", fs[0], ""))
            )
            add(
                sid,
                f"line {s.count(chr(10),0,m.start())+1}",
                None,
                fn,
                au[1],
                au[2],
                "tracked_edit",
                edit=ed,
                anch=a,
                alias=alias,
            )
        else:
            items = list(
                re.finditer(r"\*\*(Comment|Reply) — (\w+) — ([^:\n]+:\d{2} [ap]m):\*\*", block, re.I)
            )
            for j, k in enumerate(items):
                sub = block[k.end() :]
                fs = fenced(sub)
                if not fs:
                    raise ValueError((alias, j))
                add(
                    sid,
                    f"line {s.count(chr(10),0,m.end()+k.start())+1}",
                    None,
                    fn,
                    k[2],
                    k[3],
                    "comment",
                    body=fs[0],
                    anch=a.copy(),
                    alias=alias,
                    thread="report:" + alias,
                )


def add_early():
    sid = "early_inventory"
    s = txt(sid, DIS + "summary_inventory_ajh26.md", JULSEP)
    parts = list(re.finditer(r"^(## [^\n]+|### [^\n]+)", s, re.M))
    fn = ""
    for i, m in enumerate(parts):
        title = m[0]
        if title.startswith("## "):
            fn = title[3:]
            continue
        block = s[m.end() : parts[i + 1].start() if i + 1 < len(parts) else len(s)]
        date = re.search(r"Date: ([^\n]+)", block)[1].strip()
        ln = re.search(r"Lines?: ([^\n]+)", block)
        fs = fenced(block)
        a = anchor(
            "transcribed source selection" if fs else "no exact selection",
            fs[0] if title == "### Comment" and fs else None,
            line=ln[1].strip() if ln else None,
        )
        loc = f"line {s.count(chr(10),0,m.start())+1}"
        alias = "early:" + str(i)
        if title == "### Comment":
            c = re.search(r'ajh26 commented: "(.*)"', block, re.S)
            if not c:
                raise ValueError((i, title, block))
            body = c[1].split('"\n', 1)[0]
            add(
                sid,
                loc,
                None,
                fn,
                "ajh26",
                date,
                "comment",
                body=body,
                anch=a,
                alias=alias,
                quality="dated advisor inventory transcript",
            )
        else:
            if "replacement" in title:
                ed = ce("replace", *fs[:2])
            elif "addition" in title:
                ed = ce("insert", "", fs[0])
            else:
                ed = ce("delete", fs[0], "")
            add(
                sid,
                loc,
                None,
                fn,
                "ajh26",
                date,
                "tracked_edit",
                edit=ed,
                anch=a,
                alias=alias,
                quality="dated advisor inventory transcript",
            )


REMOTE = {
    "ch4_review.txt": "2fc717c593d6fb917f452460918ad332fe7ae9c8344ac9ec2dddb2f7a72b9616",
    "ch5_review.txt": "2c9faf6827a829288e4c75c9b73cd371abe5bcd317434991de37bb0a1861c7be",
    "comment_cleanup_final_20261028.json": "681c64c1afc99183c3084efcafcb0243b37c176ac219a76b5e7a69604059e0ac",
    "comment_resolutions_20261028.json": "284dbf0929438d1c3dd151d24b56f002d3f536d8937e95c2135890adde2a6437",
    "dissertation_torque_comments.json": "7c9f26af7215962bc89a28be3b543112a9c3cb94d0b379a299657a26593a59af",
    "resolved_comments_view.txt": "bb98083289840a505826257753f9ab47b473d1ce4f9395ebbf709a3a6ea71cad",
    "remaining_comments_20261028.txt": "0fc9eb1b42da386a787d8138cbf588e3bedf4d493f9b7b6325a45f0d8f7db042",
    "dissertation_review_initial.txt": "8721b0cd0d24317be9a579a3e6951f0d7ead69eef60eb6c6d2ba70775cd91e13",
}


def remote_inputs(old):
    for fname, sha in REMOTE.items():
        folder = "ch45_polish_20261029" if fname.startswith(("ch4", "ch5")) else "reconciliation_20261028"
        path = DIS + "Archive/" + folder + "/" + fname
        sid = "eb_" + fname.split(".")[0]
        cache = ROOT / "tmp/history_review_20261118/reader_raw" / fname
        if (ROOT / path).is_file():
            data = (ROOT / path).read_bytes()
        elif cache.is_file():
            data = cache.read_bytes()
        elif sid in old:
            data = zlib.decompress(base64.b64decode(old[sid]["embedded_zlib_base64"]))
        else:
            raise RuntimeError(
                "Restore the hash-checked source or keep the distributed JSON beside this generator: " + path
            )
        assert digest(data) == sha, fname
        RAW[sid] = data
        SOURCES[sid] = dict(
            id=sid,
            path=path,
            origin="EB475WS4 archive; raw bytes embedded",
            commit=None,
            git_blob=None,
            sha256=sha,
            immutable_ref="sha256:" + sha,
            embedded_zlib_base64=base64.b64encode(zlib.compress(data, 9)).decode(),
        )


def scalar(v):
    if len(v) > 1 and v[0] == v[-1] == '"':
        try:
            return json.loads(v)
        except json.JSONDecodeError:
            return v[1:-1]
    if len(v) > 1 and v[0] == v[-1] == "'":
        return v[1:-1].replace("''", "'")
    return v


def paragraphs(data):
    with ZipFile(BytesIO(data)) as z:
        r = ET.fromstring(z.read("word/document.xml"))
    w = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
    return ["".join(t.text or "" for t in p.iter(w + "t")) for p in r.iter(w + "p")]


def add29():
    sid = "29_comment_transcript"
    ps = paragraphs(
        source(
            sid,
            DIS + "Notes/committee_revision_20261027/CH45_COMMENTS_REPLIES_AND_FIXES_REVISED_20261029.docx",
        )
    )
    ids = jfile("30_live_after", DIS + "Notes/second_pdf_comments_20261030/LIVE_REVIEW_AFTER.json")[
        "comments"
    ]
    indices = [i for i, p in enumerate(ps) if re.fullmatch(r"Comment [45] \d+", p)]
    for idx, start in enumerate(indices):
        block = ps[start : indices[idx + 1] if idx + 1 < len(indices) else len(ps)]
        card = block[0].replace("Comment ", "").replace(" ", ".")
        f = "30-torque.tex" if card.startswith("4.") else "50-database.tex"
        authorline = next((p for p in block if p.startswith("Advisor comment ")), None)
        if not authorline:
            EXCLUDED[(sid, "bbolen83", "comment")] += 1
            continue
        m = re.match(r"Advisor comment\s+(\d+ \w+, \d+:\d+ [ap]m)\s+(.*)", authorline)
        assert m, authorline
        line = next(
            (
                re.search(r"Source: [^,]*, line (\d+)", p)[1]
                for p in block[:8]
                if re.search(r"Source: [^,]*, line (\d+)", p)
            ),
            None,
        )
        ctx = next((p for p in block[:7] if p.startswith("Ben’s saved document instruction")), None)
        a = anchor(
            "transcript's historical local revision reference"
            if line
            else "review card; original selection unavailable",
            context=ctx,
            line=int(line) if line else None,
        )
        pid = None
        for r in ids:
            lines = r["text"].splitlines()
            if lines[1] != m[1] or not m[2].startswith(lines[6].removesuffix("...")):
                continue
            pid = r["thread"].split("-")[-1]
        add(
            sid,
            f"word/document.xml paragraph {start+block.index(authorline)+1} (1-based); card {card}",
            pid,
            f,
            "ajh26",
            m[1],
            "comment",
            body=m[2],
            anch=a.copy(),
            alias="card-" + card,
            pid=pid,
            thread=pid or ("card-29-" + card),
            quality="expanded advisor-comment transcript",
        )
        if pid:
            EVENTS[-1]["sources"].append(
                dict(
                    source_id="30_live_after",
                    locator="/comments/"
                    + str(next(i for i, v in enumerate(ids) if v["thread"].endswith(pid))),
                )
            )
        for j, t in enumerate(block):
            if not t.startswith("Advisor follow up "):
                continue
            follow = re.fullmatch(r"Advisor follow up\s+(\d+ \w+, \d+:\d+ [ap]m)\s+(.*)", t)
            assert follow, t
            add(
                sid,
                f"word/document.xml paragraph {start+j+1} (1-based); card {card} follow-up",
                None,
                f,
                "ajh26",
                follow[1],
                "comment",
                body=follow[2],
                anch=a.copy(),
                alias="card-" + card + "-follow-up",
                thread=pid or ("card-29-" + card),
                quality="expanded advisor-comment transcript",
            )
    assert len(indices) == 42
    assert sum(v.startswith("Advisor follow up ") for v in ps) == 2


def aria(sid):
    raw = RAW[sid].decode("utf-8-sig").splitlines()
    cards = [i for i, v in enumerate(raw) if re.match(r"^      - '?button (?:\"|\[active\])", v)]
    f = "30-torque.tex" if sid == "eb_ch4_review" else "50-database.tex" if sid == "eb_ch5_review" else None
    for j, start in enumerate(cards):
        title = raw[start]
        block = raw[start + 1 : cards[j + 1] if j + 1 < len(cards) else len(raw)]
        stop = next((i for i, b in enumerate(block) if b.strip() and len(b) - len(b.lstrip()) < 8), None)
        if stop is not None:
            block = block[:stop]
        filename = re.search(r'button "([\w-]+\.tex) \d+', title)
        if filename:
            f = filename[1]
            continue
        m = re.search(r'button "(\w+) (\d+ \w+, \d+:\d+ [ap]m) (Accept change|Resolve comment)', title)
        if not m and title.strip() == "- button [active]:" and len(block) >= 3:
            subject = block[0].removeprefix("        - generic: ")
            when = block[1].removeprefix("        - generic: ")
            action = next(
                (a for a in ("Accept change", "Resolve comment") if f'        - button "{a}":' in block), None
            )
            if action:
                m = re.fullmatch(
                    r"(\w+) (\d+ \w+, \d+:\d+ [ap]m) (Accept change|Resolve comment)",
                    f"{subject} {when} {action}",
                )
        if not m or not f:
            continue
        loc = f"lines {start+1}-{start+len(block)+1}"
        if m[3] == "Accept change":
            if m[1] != "ajh26":
                EXCLUDED[(sid, m[1], "tracked_edit")] += 1
                continue
            op = next(
                (
                    act
                    for act, tag in [
                        ("replace", '"Changed:"'),
                        ("insert", '"Added:"'),
                        ("delete", '"Deleted:"'),
                    ]
                    if any(b.strip() == "- text: " + tag for b in block)
                ),
                None,
            )
            if op is None:
                continue

            def side(lines):
                vals = [
                    scalar(q[1])
                    for b in lines
                    if (q := re.match(r"^\s+- (?:deletion|insertion|generic|text): (.+)$", b))
                    and q[1] not in ["show more", "show less"]
                ]
                return "".join(vals)

            k = next(i for i, b in enumerate(block) if re.search(r'- text: "(?:Changed|Added|Deleted):"', b))
            sides = block[k + 1 :]
            to = next(
                (
                    i
                    for i, q in enumerate(sides)
                    if q == "        - text: to" or q.startswith("        - text: to ")
                ),
                None,
            )
            sides1 = sides[:to] if to is not None else sides
            sides2 = sides[to + 1 :] if to is not None else []
            a1 = side(sides1)
            a2 = side(sides2)
            if to is not None and sides[to].startswith("        - text: to "):
                a2 = scalar(sides[to][len("        - text: to ") :]) + a2
            ed = ce(op, "" if op == "insert" else a1, a1 if op == "insert" else a2)
            partial = any('button "show more"' in v for v in sides)
            if partial:
                ed["displayed_fragments"] = {k: v for k, v in ed.items() if k in {"original", "proposed"}}
                if any('button "show more"' in v for v in sides1):
                    ed["proposed" if op == "insert" else "original"] = None
                if any('button "show more"' in v for v in sides2):
                    ed["proposed"] = None
            add(
                sid,
                loc,
                None,
                f,
                m[1],
                m[2],
                "tracked_edit",
                edit=ed,
                quality="distinct edit card; one or more sides collapsed"
                if partial
                else "raw expanded or short Overleaf review card",
            )
        else:
            who = [
                i
                for i, b in enumerate(block)
                if re.fullmatch(r"        - generic: (?:ajh26|bbolen83|mcneal|leow)", b)
            ]
            for n, idx in enumerate(who):
                msg = block[idx : who[n + 1] if n + 1 < len(who) else len(block)]
                author = msg[0].split(": ", 1)[1]
                if author != "ajh26":
                    EXCLUDED[(sid, author, "comment")] += 1
                    continue
                time = next(
                    (
                        re.match(r"        - generic: (\d+ \w+, \d+:\d+ [ap]m)$", t)[1]
                        for t in msg
                        if re.match(r"        - generic: (\d+ \w+, \d+:\d+ [ap]m)$", t)
                    ),
                    None,
                )
                if not time:
                    continue
                opts = [
                    scalar(v[1])
                    for b in msg[2:]
                    if (v := re.fullmatch(r"        - (?:generic|text): (.+)", b))
                    and v[1] not in ("show more", "show less", "Reply")
                ]
                if not opts:
                    continue
                body = "\n".join(opts)
                partial = any('button "show more"' in b for b in msg)
                add(
                    sid,
                    loc + f" message {n+1}",
                    None,
                    f,
                    author,
                    time,
                    "comment",
                    body=body,
                    quality="collapsed comment preview"
                    if partial
                    else "raw expanded or short Overleaf review card",
                )


def svalue(r):
    if r["type"] == "comment":
        return norm(r["comment"])
    e = r["edit"]
    return "|".join(
        (
            e["operation"],
            norm(e["original"] if e["original"] is not None else e["displayed_fragments"]["original"]),
            norm(e["proposed"] if e["proposed"] is not None else e["displayed_fragments"]["proposed"]),
        )
    )


def is_preview(r):
    return "collapsed" in r["evidence_quality"] and r["type"] == "comment"


def idbasis(r):
    return (r["date"], r["file"], r["time_local"] or "", r["type"], svalue(r))


def combine(a, b):
    a["sources"].extend(x for x in b["sources"] if x not in a["sources"])
    a["aliases"].extend(x for x in b["aliases"] if x not in a["aliases"])
    for attr in ("platform_id", "thread_id", "time_local"):
        if not a[attr]:
            a[attr] = b[attr]
    for k, v in b["anchor"].items():
        if v is not None and a["anchor"].get(k) is None:
            a["anchor"][k] = v
            a["anchor"].setdefault("additional_field_sources", {})[k] = b["anchor"].get("source_id")
    if (b["anchor"]["text"] or b["anchor"]["context"]) and a["anchor"]["status"] == "unavailable":
        a["anchor"]["status"] = b["anchor"]["status"]


def dedup():
    counts = Counter()
    allrows = []
    ids = {}
    pending = []
    for r in EVENTS:
        sid = r["sources"][0]["source_id"]
        if (
            (sid == "early_inventory" and r["aliases"][0] in {"early:56", "early:57"})
            or (
                r["date"] == "2026-10-14"
                and r["type"] == "tracked_edit"
                and r["edit"]["operation"] == "delete"
                and norm(r["edit"]["original"]) == "these tools were"
            )
            or (
                sid == "eb_ch5_review"
                and r["type"] == "tracked_edit"
                and r["date"] == "2026-10-23"
                and "displayed_fragments" in r["edit"]
            )
        ):
            why = (
                "The original inventory has two identical question marks on two distinct selections. Two separately timed messages C31/C32 survive without their original selections. The correspondence of the two original selections cannot be established."
                if sid == "early_inventory"
                else "Later capture of an altered or collapsed edit overlaps the original review at the identical file, author, minute and wording. Separate origin is unverified."
            )
            GAPS.append(
                dict(
                    type="unresolved_version_or_identity_overlap",
                    reason=why,
                    date=r["date"],
                    file=r["file"],
                    time_local=r["time_local"],
                    anchor=r["anchor"],
                    edit=r["edit"],
                    comment=r["comment"],
                    sources=r["sources"],
                    aliases=r["aliases"],
                )
            )
            continue
        base = idbasis(r)
        ordinal = counts[(sid, base)]
        counts[(sid, base)] += 1
        r["occurrence"] = ordinal
        if is_preview(r):
            pending.append(r)
            continue
        choices = [
            e
            for e in allrows
            if e["type"] == r["type"]
            and e["date"] == r["date"]
            and e["occurrence"] == ordinal
            and svalue(e) == svalue(r)
            and (not e["time_local"] or not r["time_local"] or e["time_local"] == r["time_local"])
            and (not e["platform_id"] or not r["platform_id"] or e["platform_id"] == r["platform_id"])
        ]
        samefile = [e for e in choices if e["file"] == r["file"]]
        if samefile:
            choices = samefile
        elif r["type"] == "tracked_edit" or not r["time_local"]:
            choices = []
        if r["platform_id"] in ids:
            existing = ids[r["platform_id"]]
        elif len(choices) == 1:
            existing = choices[0]
        else:
            existing = None
        if existing is None:
            allrows.append(r)
            existing = r
        else:
            combine(existing, r)
        if r["platform_id"]:
            ids[r["platform_id"]] = existing
    for r in pending:
        snippet = norm(r["comment"]).removesuffix("...").removesuffix("…")
        cand = [
            e
            for e in allrows
            if e["type"] == "comment"
            and e["date"] == r["date"]
            and e["file"] == r["file"]
            and (not e["time_local"] or e["time_local"] == r["time_local"])
            and (svalue(e).startswith(snippet) or snippet.startswith(svalue(e)))
        ]
        if len(cand) == 1:
            combine(cand[0], r)
        else:
            GAPS.append(
                dict(
                    type="collapsed_comment",
                    reason="Full message or unique match unavailable.",
                    file=r["file"],
                    date=r["date"],
                    time_local=r["time_local"],
                    preview=r["comment"],
                    sources=r["sources"],
                )
            )
    for e in allrows:
        month = {"07": "JUL", "09": "SEP", "10": "OCT"}[e["date"][5:7]]
        short = "C" if e["type"] == "comment" else "T"
        sig = json.dumps([*idbasis(e), e["occurrence"]], ensure_ascii=False, separators=(",", ":"))
        e["id"] = f"{month}-{short}-" + (e["platform_id"] or digest(sig.encode())[:16])
        e["identity_basis"] = (
            "Overleaf message ID"
            if e["platform_id"]
            else "date, file, visible local time when available, exact event wording, and within-capture occurrence"
        )
        e["project"] = "Bolen-Dissertation"
        e.pop("event_uid")

    def time_sort(t):
        return datetime.strptime(t, "%I:%M %p").strftime("%H:%M") if t else ""

    return sorted(
        allrows, key=lambda r: (r["date"], r["file"], time_sort(r["time_local"]), r["type"], r["id"])
    )


def extract():
    destination = HERE / "HISTORICAL_SIGNALS.json"
    old = json.loads(destination.read_text(encoding="utf-8"))["sources"] if destination.exists() else {}
    remote_inputs(old)
    add23()
    add22()
    add15()
    add_early()
    add29()
    for sid in ["eb_ch4_review", "eb_ch5_review", "eb_resolved_comments_view"]:
        aria(sid)
    return dedup()


def inspection(rows):
    for i, e in enumerate(rows, 1):
        val = e["comment"] or (str(e["edit"]["original"]) + " -> " + str(e["edit"]["proposed"]))
        if e["edit"] and "displayed_fragments" in e["edit"]:
            val = str(e["edit"]["displayed_fragments"])
        print(
            f"{i:03d} {e['date'][-5:]} {e['file'][:3]} {e['time_local'] or '-':8s} {e['type'][0]} {' '.join(val.split())}"[
                :640
            ]
        )


# Each line labels the ten consecutively numbered original inputs indicated in
# its comment. The binding below verifies ordered canonical IDs before use.
# An integer is a comment score magnitude. Every tracked edit has weight −1.
LABELS = """
2gap+flow 2context+flow 2flow 1short short+scope 1short+flow short scope test scope
term scope copy 2gap 2flow+credit 1copy 2flow 2flow 2context+fig 3credit
3unit 3num+term 1term 0praise 1fig scope scope term short 2fig+cite
2credit+gap 1flow 1copy 2credit+gap 1flow 3credit 2cite+credit 1term 2credit 2short+flow
2scope+use 2use+test scope+short short+flow scope+short copy scope 3num+test 3num+context 1term
2test+context 1test+praise 2test 2mech+context 3mech+num 2mech+gap 1flow 2scope+mech 2mech 1flow
2context+mech 1context 1context 0praise 2context+mech 1flow 1flow+copy 2use+credit 2gap+praise 2gap+flow
scope 1flow 1credit 1term 1fig 2fig 1copy credit credit 1flow
1flow copy 2scope+cite 2flow+result scope 2scope 2gap+test  mech+scope 2gap+credit+flow 1term
1term 1mech+term 1mech 1mech 1short+copy 2copy+flow 2context+term 2scope+flow 1context 2gap+context
3scope 3scope short 1flow 2gap+credit 2flow 3credit+cite 1term 1term 1term
3mech 1term 1context 1flow 1term flow 2flow+test 2use+result+short 1term 2context+use+flow
1credit 2test 2mech+test 2credit 1context 2gap+credit 1context 1term copy term
2gap+use 2gap+flow 1context 1fig 1context short+credit 1short 2gap+context 2cite+test  scope+flow
scope term+use 1term term copy scope 2context 2num+test copy 1copy
term 2flow+gap 2credit+context 1term 1term 1term 1term copy use copy
1short  copy 1context use credit 1term test test+short 1term copy
use+gap+short term copy short 3fig+num 2fig+short 1term 1copy term 3test+scope
1context 1copy 1copy 1context copy 2mech+term 1term 1copy 1flow flow
term+credit fig 1flow short context 1term 1mech copy 1short 1short
2flow 2test 1mech short fig+short short+flow term copy scope term
copy copy short+scope copy 2flow+scope 2scope+flow 1flow copy context copy
1flow short+fig gap term scope term test use+gap gap+use 2test
2scope+gap 1context copy short short copy short short+scope flow flow
flow cite+term flow credit+mech mech+use copy mech term copy copy
copy context 1flow 1context 1context scope+short gap+credit use copy 1context
2mech+term 2flow+context copy term term copy use+term short term short+credit
scope scope cite+credit term 2context short 1context 1context 1term short
mech short 2mech+term 1term context flow flow 2scope+use credit+context 1term
1credit 1credit 1short 2result+test 1context 1context credit+use mech+test+credit term credit
copy term test+credit 1context copy credit use short copy 1term
1copy 1copy 1context 1term 2test+mech 2term+flow 1fig 2result+fig 3test+num 1unit
3test+scope 1short short 2fig+test short+scope 2result+fig flow 2context+gap 1term 1use+term
1use+term term 2use+result mech context term copy term term copy
term copy scope credit 1term 2mech+test scope scope 1context 2mech
1short 1context 1term 2short+scope 1term 2context+use 1context num+scope 1term 3scope+test
scope 1context 1copy 2flow 2flow+use 2fig+use 1use+fig 1context 1credit 1term
use+term context 3num 1mech+term short+fig fig 3cite+credit copy copy copy
test+flow short 1result+praise 2use+flow 2flow 2flow 2use+mech
"""
ORDER_HASH = "63a5b417af15045a91f7d6984c0b9d2d8e36cafbe3af80d1e084db339f160393"
NOTES = {
    1: "Give the dissertation a scientific motivation and connect the chapter map to that purpose.",
    14: "The reader needs an explicit sentence identifying what remains unknown.",
    20: "The lab credited for the cited work is factually wrong.",
    24: "The advisor endorses this particular passage without requesting a change.",
    39: "Assign the technical advance to its actual developer.",
    41: "Separate completed work from proposed experiments and name the present benefit.",
    52: "The positive assessment is qualified by a request for another person to check the equations.",
    55: "The stated independent variables have a dependency that changes how the model should be represented.",
    64: "The alternative introduction receives praise without a requested change.",
    68: "Describe the tool developed, who built it, and what it lets a reader do.",
    69: "Praise for one stated gap is followed by a request to identify the others.",
    83: "The proposed ethics argument needs the dependency on animal experiments.",
    84: "A methods-only sequence makes the reader carry too much. The advisor asks for topic chapters with their own outcomes.",
    87: "The reader cannot reconstruct the rationale, design sequence, or scientific work from the current account.",
    89: "Introduce what is being built and why, then distinguish background from dissertation contributions.",
    92: "Replace an ambiguous implementation verb with the actual operation and language.",
    97: "Give the mechanical-engineering reader sufficient neuroscience context or explain the models in accessible language.",
    101: "The described alternatives are incomplete, limiting the apparent scope of the database field.",
    102: "The comment disputes the apparent completeness of the classification.",
    105: "State the missing knowledge then identify what the author contributed to fill it.",
    107: "The account of prior presentations appears to omit a separate version.",
    118: "Lead with the finished resource and retain the methods needed to recreate it.",
    120: "The reader needs the database contents and organization before the curation details.",
    124: "Replace the indefinite actor with the author's own role.",
    126: "Technical detail is obscuring both the research question and what the author built.",
    132: "The advisor requests the gap and the work that fills it, with implications moved to Discussion.",
    138: "Explain what prior tools do well enough to establish the specific missing capability.",
    139: "The reader asks for a direct comparison and attribution to the existing libraries.",
    152: "Present the knowledge gaps and topic chapters in matching order.",
    153: "Name the existing platform and the technical advance this dissertation contributed.",
    167: "The edit replaces a demonstration label with the validation purpose.",
    168: "The edit replaces a list of exercises with a concrete assembly and validation sequence.",
    171: "The short replacement explains how the four topic chapters support the larger objective.",
    175: "The advisor identifies a visibly misplaced center in a technical schematic.",
    180: "Identify whose data were collected before interpreting their relevance.",
    199: "File-saving detail does not explain the engineering.",
    202: "The passage claims a comparison without giving the compared entities.",
    223: "Frame the scientific problem as a missing answer that can be supplied.",
    227: "The tracked replacement makes the stated purpose validation.",
    228: "The addition identifies the researcher's difficulty using the prior material.",
    229: "Replace a negative list of operations with the missing complete model and the hypotheses it would let a researcher test.",
    257: "State the historical gap from the author's perspective.",
    258: "Name the design activity a reader can carry out with the resource.",
    270: "The deletion removes maintenance logistics from the scientific explanation.",
    271: "An insertion and deletion of the same scoped sentence are separately visible actions. They may represent a move.",
    272: "An insertion and deletion of the same scoped sentence are separately visible actions. They may represent a move.",
    283: "Distinguish the direction of a citation from the direction of the cited-by relationship.",
    288: "Name the unfinished relational links between pathways, reviews, and papers.",
    297: "The advisor replaces a missing-combinations inventory with the author's purpose and first-person contributions.",
    298: "The new text names the source model, conversion, connection, and controller test.",
    303: "Identify what was tested and introduce the next system with the author as actor.",
    307: "State the capability the tools enable.",
    311: "The one-character message is treated as a minimal editorial correction. Its intended change is indeterminate.",
    315: "Report a concrete defect, how it was detected, and how it was corrected if one occurred. The advisor's examples are hypothetical.",
    316: "Introduce the BPA interface and connect it to the preceding model before discussing its effect.",
    318: "State the engineering result before directing the reader to the figure.",
    319: "The advisor questions the interpretation and suggests a zero offset. The suggestion does not establish that explanation.",
    321: "The advisor questions whether a fit metric diagnoses out-of-sample model tests. Ben's response and decision are separate evidence.",
    326: "Lead with the result, then point to the supporting display.",
    333: "Name the reader-facing Explorer as the principal deliverable at this location.",
    344: "The inserted wording explicitly identifies the advisor's contribution.",
    350: "The reader cannot tell why records entered the collection and were then removed.",
    354: "The advisor objects to an unnecessary statement about missing information and describes its effect on reader confidence. The message contains no finding of misconduct.",
    360: "Coverage of pathway-level classification is materially smaller than the whole database. Qualify the claim about the finished resource.",
    364: "Give the paragraph a topic sentence that distinguishes all views from the table view.",
    365: "Group the map and table explanations by view so that the reader can follow the example.",
    367: "Explain whether the interface changes axes or toggles between views.",
    373: "The advisor flags an apparently inconsistent count. The historical record does not confirm a corrected value.",
    377: "The historical account needs both presentations and their source citations.",
    383: "The advisor praises the features while requesting a more positive opening.",
    384: "Introduce the database before explaining its use to design testable control models.",
    387: "State a specific research use, then walk the reader through how the database supports it.",
}

GROUPS = {
    "OCT-I-PASSIVE-TO-AUTHOR-1014": (78, 79),
    "OCT-I-GAPS-AND-CHAPTER-ORDER-1022": (152, 154, 155, 156, 157),
    "OCT-I-VALIDATION-BEFORE-ASSEMBLY-1022": (167, 168),
    "OCT-I-CPG-PARAGRAPH-PLACEMENT-1022": (215, 216, 217),
    "OCT-I-SADB-ORIGINAL-GAP-1023": (256, 257),
    "OCT-I-SADB-CURATION-CLAIM-1023": (271, 272),
    "OCT-I-SNS-PURPOSE-AND-DELIVERABLES-1023": (297, 298, 299, 300, 301, 302, 303, 305, 306, 307),
    "OCT-I-MODEL-CONVERSION-CHECKS-1023": (313, 314, 315),
    "OCT-I-HYBRID-POINTS-1029": (324, 325),
    "OCT-I-EXPLORER-VIEWS-1029": (355, 356),
    "OCT-I-DATABASE-RESEARCH-USE-1029": (385, 387),
}


def classify(rows):
    extra = [
        r
        for r in rows
        if r["date"] == "2026-10-29"
        and r["time_local"] == "10:13 am"
        and "card-5.23-follow-up" in r["aliases"]
    ]
    assert len(extra) == 1, "Check October 29 chapter 5 follow-up identity"
    base = [r for r in rows if r is not extra[0]]
    assert (
        digest("\n".join(r["id"] for r in base).encode()) == ORDER_HASH
    ), "Input identities changed. Reconcile the classification against the original inputs."
    labels = LABELS.split()
    assert len(labels) == len(base)
    lookup = {i: group for group, indices in GROUPS.items() for i in indices}
    groups = {}
    for i, (r, lab) in enumerate(zip(base, labels), 1):
        m = re.fullmatch(r"([0-3]?)([a-z+]+)", lab)
        assert m and bool(m[1]) == (r["type"] == "comment"), (i, r["type"], lab)
        codes = m[2].split("+")
        r["tags"] = [TAGS[c][0] for c in codes]
        r["score"] = -int(m[1]) if m[1] else -1
        vague = r["type"] == "comment" and norm(r["comment"]) in {
            "?",
            "huh?",
            "*",
            "huh",
            "i'm lost.",
            "unclear",
            "i don't know what this means",
            "i don't really know what this means",
            "which?",
        }
        r["interpretation_confidence"] = (
            "bounded"
            if vague or "collapsed" in r["evidence_quality"] or r["anchor"]["status"] == "unavailable"
            else "direct or corroborated"
        )
        if i in NOTES:
            r["rationale"] = NOTES[i]
        elif vague:
            r[
                "rationale"
            ] = "The message records a reader difficulty. Its wording alone does not identify the exact technical defect."
        elif r["type"] == "tracked_edit":
            r["rationale"] = (
                "The "
                + r["edit"]["operation"]
                + " action addresses "
                + r["tags"][0].replace("_", " ")
                + ". "
                + TAGS[codes[0]][1]
            )
        else:
            r["rationale"] = TAGS[codes[0]][1]
        issue = (
            lookup.get(i)
            or (groups.get(r["thread_id"]) if r["thread_id"] else None)
            or r["id"].replace("-C-", "-I-").replace("-T-", "-I-")
        )
        r["issue_id"] = issue
        if r["thread_id"]:
            groups[r["thread_id"]] = issue
        r["occurrence_in_capture"] = r.pop("occurrence") + 1
    r = extra[0]
    r.update(
        score=-2,
        tags=[TAGS[c][0] for c in ("mech", "scope", "use")],
        rationale="The advisor sees apparently removed or irrelevant records in the curated app and asks for an explanation of which records the reader can see.",
        interpretation_confidence="direct or corroborated",
        issue_id=groups[r["thread_id"]],
        occurrence_in_capture=r.pop("occurrence") + 1,
    )


def attach_sources(rows):
    led = jfile("15_ledger", DIS + "Notes/revision_campaign_20261015/comment_ledger.json")
    for i, item in enumerate(led):
        for r in rows:
            if item["id"] in r["aliases"] and item["author"] == r["author"]:
                r["sources"].append(
                    dict(source_id="15_ledger", locator=f"/{i}; derived campaign cross-reference only")
                )
    sms = jfile(
        "26_sns_replies",
        DIS + "Notes/advisor_revision_20261023/Revision_E/review/overleaf_latest_sns_replies.json",
        OCT26,
    )
    original = sms["new_reply"]
    rr = [
        r
        for r in rows
        if r["file"] == "60-sns.tex"
        and r["comment"] == original["comment"]
        and (r["date"], r["time_local"]) == pdate(original["advisor_time"])
    ]
    assert len(rr) == 1 and original["advisor"] == "ajh26"
    rr[0]["sources"].append(
        dict(source_id="26_sns_replies", locator="/new_reply; October 26 reply belongs to bbolen83")
    )
    cleanup = json.loads(RAW["eb_comment_cleanup_final_20261028"].decode("utf-8-sig"))
    a = Counter()
    for i, x in enumerate(cleanup["newResolutions"]):
        if x["document"] != "dissertation":
            a["other_project"] += 1
            continue
        q = x["thread"]
        if q.startswith("bbolen83 "):
            a["author_ben"] += 1
            continue
        if q == "23 October 11:28 huh?":
            c = [
                r
                for r in rows
                if r["file"] == "60-sns.tex"
                and r["time_local"] == "11:28 am"
                and r["comment"] == "huh?"
                and r["date"] == "2026-10-23"
            ]
        else:
            m = re.match(r"ajh26 (\d+ \w+, \d+:\d+ [ap]m) Resolve comment More options (.*)", q)
            if m:
                d, t = pdate(m[1])
                partial = re.split(
                    r"(?:\.\.\. show more|\s+(?:bbolen83|ajh26|mcneal) \d+ \w+,| Reply\s*$)", m[2]
                )[0]
                c = [
                    r
                    for r in rows
                    if r["date"] == d
                    and r["time_local"] == t
                    and r["type"] == "comment"
                    and norm(r["comment"]).startswith(norm(partial))
                ]
            else:
                c = []
        ref = dict(source_id="eb_comment_cleanup_final_20261028", locator=f"/newResolutions/{i}")
        if q == "This belongs in the other section then.":
            hits = [
                j
                for j, u in enumerate(cleanup["newResolutions"])
                if u["document"] == "dissertation"
                and u["thread"].startswith("bbolen83 ")
                and (q + " Reply") in u["thread"]
            ]
            assert len(hits) == 1
            GAPS.append(
                dict(
                    type="excluded_non_advisor_resolution",
                    reason="This authorless export repeats Ben's attributed comment in the same archive.",
                    author="bbolen83",
                    preview=q,
                    sources=[
                        ref,
                        dict(
                            source_id="eb_comment_cleanup_final_20261028",
                            locator=f"/newResolutions/{hits[0]}",
                        ),
                    ],
                )
            )
            a["duplicate_ben_excerpt"] += 1
            continue
        if len(c) == 1:
            c[0]["sources"].append(ref)
            a["uniquely_corroborated"] += 1
        else:
            why = (
                "Ambiguous original review event."
                if len(c) > 1
                else "Resolution excerpt lacks verified original message identity or text."
            )
            GAPS.append(
                dict(
                    type="resolution_control",
                    reason=why,
                    preview=q,
                    possible_ids=[o["id"] for o in c],
                    sources=[ref],
                )
            )
            a["unresolved_mapping"] += 1
    more = json.loads(RAW["eb_comment_resolutions_20261028"].decode("utf-8-sig"))
    ben = [
        (i + 1, norm(s))
        for i, s in enumerate(RAW["eb_resolved_comments_view"].decode("utf-8-sig").splitlines())
        if re.match(r'^      - \'?button "bbolen83 \d+ \w+, \d+:\d+ [ap]m Resolve comment', s)
    ]
    for i, q in enumerate(more["dissertation"]):
        c = [r for r in rows if r["type"] == "comment" and norm(q) == svalue(r)]
        ref = dict(source_id="eb_comment_resolutions_20261028", locator=f"/dissertation/{i}")
        if q.startswith("bbolen83\n"):
            identifier = norm(" ".join(q.splitlines()[:2]))
            body = norm(q.splitlines()[-1])
            candidates = [
                n
                for n, display in ben
                if identifier in display and ("more options " + body + " reply") in display
            ]
        else:
            candidates = [n for n, display in ben if norm(q)[:70] in display]
        if len(candidates) == 1:
            GAPS.append(
                dict(
                    type="excluded_non_advisor_resolution",
                    reason="The archived dissertation Review panel identifies Ben as the author of the resolved comment.",
                    author="bbolen83",
                    preview=q,
                    sources=[
                        ref,
                        dict(source_id="eb_resolved_comments_view", locator=f"line {candidates[0]}"),
                    ],
                )
            )
        elif len(c) == 1:
            c[0]["sources"].append(ref)
        else:
            GAPS.append(
                dict(
                    type="unattributed_resolution",
                    reason="The resolution export gives text without an original advisor identifier. A unique original advisor message could not be confirmed.",
                    preview=q,
                    sources=[ref],
                )
            )
    for sid, path in {
        "chair": "ProofFinal/main.tex",
        "campaign_rules": "Notes/revision_campaign_20261015/RULINGS_20261016.md",
        "oct22_teaching": "summary.md",
        "oct22_rounds": "Notes/advisor_revision_20261022/ROUND4_NOTES.md",
        "oct29_annotation": "Notes/committee_revision_20261027/CH45_ANNOTATED_REVIEW_20261029.md",
        "oct27_committee": "Notes/committee_revision_20261027/REVIEW_REPORT.md",
        "oct30_pdf_user": "Notes/second_pdf_comments_20261030/USER_COMMENTS_19.json",
        "oct30_paper": "Notes/second_pdf_comments_20261030/PAPER_LIVE_REVIEW_AFTER.json",
        "oct30_defense": "Defense/AJH_Review_20261030/comments.json",
        "oct30_weighted_rubric": "Notes/publication_audit_20261030/REVIEW_RUBRIC.md",
    }.items():
        source(sid, DIS + path)
    source("23_campaign_report", "CHATGPT_REPORT.md", MIRROR)
    source(
        "26_torque_entries",
        "Documentation/Reports and Papers/Knee_Torque_Test/review_audit_20261026/dissertation_torque_review_entries.json",
    )
    for r in rows:
        for alias in r["aliases"]:
            if re.fullmatch("M-[0-9]{3}", alias):
                r["sources"].append(
                    dict(
                        source_id="23_matrix",
                        locator=f"/{int(alias[2:])-1} ({alias}); historical context/readiness only",
                    )
                )
    return dict(a)


def tally(rows):
    return dict(
        events=len(rows),
        comments=sum(r["type"] == "comment" for r in rows),
        tracked_edits=sum(r["type"] == "tracked_edit" for r in rows),
        score=sum(r["score"] for r in rows),
        score_histogram={str(n): sum(r["score"] == n for r in rows) for n in [0, -1, -2, -3]},
        platform_message_ids=sum(bool(r["platform_id"]) for r in rows),
        tracked_edits_with_incomplete_text=sum(
            r["type"] == "tracked_edit" and "displayed_fragments" in r["edit"] for r in rows
        ),
    )


def make_output(rows):
    classify(rows)
    controls = attach_sources(rows)
    grouped = {}
    for g in GAPS:
        values = [
            g.get("date"),
            g.get("file"),
            g.get("time_local"),
            g.get("type"),
            g.get("preview"),
            g.get("edit"),
            g.get("comment"),
            g.get("aliases"),
        ]
        version_edit = (
            g["type"] == "unresolved_version_or_identity_overlap"
            and g.get("date") == "2026-10-14"
            and g.get("edit") is not None
        )
        if version_edit:
            values = [g["date"], g["file"], g["time_local"], g["type"], norm(g["edit"]["original"])]
            g["display_variants"] = [
                dict(sources=g["sources"].copy(), raw_edit=g["edit"], aliases=g["aliases"])
            ]
        key = digest(json.dumps(values, ensure_ascii=False, sort_keys=True).encode())[:16]
        prefix = {"07": "JUL", "09": "SEP", "10": "OCT"}.get((g.get("date") or "")[5:7], "OCT")
        g["id"] = prefix + "-G-" + key
        if key in grouped:
            grouped[key]["sources"].extend(ref for ref in g["sources"] if ref not in grouped[key]["sources"])
            grouped[key].setdefault("display_variants", []).extend(g.get("display_variants", []))
        else:
            grouped[key] = g
    for r in rows:
        assert r["score"] in (0, -1, -2, -3) and r["author"] == "ajh26" and r["sources"]
        if r["type"] == "tracked_edit":
            assert r["score"] == -1 and r["edit"]
        for s in r["sources"]:
            assert s["source_id"] in SOURCES
    assert len({r["id"] for r in rows}) == len(rows)
    source_counts = {}
    for sid in [
        "15_report",
        "early_inventory",
        "22_advisor_inventory",
        "23_entries_final_20261023",
        "29_comment_transcript",
        "eb_ch4_review",
        "eb_ch5_review",
        "eb_resolved_comments_view",
    ]:
        linked = [r for r in rows if any(s["source_id"] == sid for s in r["sources"])]
        source_counts[sid] = dict(
            comments_or_replies=sum(r["type"] == "comment" for r in linked),
            tracked_edits=sum(r["type"] == "tracked_edit" for r in linked),
        )
    assert source_counts["15_report"] == dict(comments_or_replies=68, tracked_edits=6)
    assert source_counts["23_entries_final_20261023"] == dict(comments_or_replies=67, tracked_edits=63)
    assert source_counts["29_comment_transcript"] == dict(comments_or_replies=44, tracked_edits=0)
    return dict(
        schema_version="1.0",
        scope=dict(
            project="Bolen-Dissertation",
            project_id=PROJECT,
            original_event_dates=sorted(set(r["date"] for r in rows)),
            reconstruction_date="2026-11-18",
            unit="One verified advisor comment message or one distinct proposed tracked edit in the dissertation. No November events included.",
            attribution=dict(
                account="ajh26",
                name="Alexander (Alex) Hunt",
                role="dissertation advisor and committee chair",
                corroboration=[
                    "chair:93 identifies Alexander Hunt as chair.",
                    "23_matrix attributes ajh26 tracked changes to Advisor.",
                    "eb_ch5_review:471–482 records ajh26 adding ‘by my advisor, Dr. Hunt’ in the author's voice.",
                    "Ben separately confirmed ajh26 as his advisor in the November review.",
                ],
            ),
            limits="This is a lower bound on the surviving, identifiable evidence. Dates are visible event dates in the retained 2026 reviews. No claims of absence or a complete review on other dates.",
        ),
        model=dict(
            identity="OpenAI Codex, GPT-6 family, independent reviewer in the active Codex desktop conversation. An exact model snapshot and API sampling parameters are unavailable.",
            task="Original-text severity, presentation family, and same-issue interpretation of actual historical advisor actions.",
            prompt="Read each ajh26 message and its originating file, date, exact selection or explicitly labeled nearby passage, and each before/after edit. Identify the reader problem evidenced by the author's own words. Code each distinct edit −1. Code a clarification or wording request −1, a substantive request −2, a material factual or interpretive challenge −3, and unqualified praise 0. A mixed message receives the negative score. Prefer the narrowest claim supported by the wording. Mark low-information messages and partial edit content as bounded. Associate independently countable signals with one issue only when the source thread or common passage supports it. Do not infer criticism or a mark from an unannotated passage, repair, resolution, campaign weight, or later ruling.",
            inputs="Each signals[] record retains the original utterance or visible edit sides, its source locators, available anchors, and an immutable source handle. Git bytes can be retrieved by immutable_ref. EB captures are embedded with SHA-256 verification.",
            output="signals[].score, tags, rationale, issue_id, interpretation_confidence. The indicated score is the model's application of Ben's rule to the original draft. It is not a numeric rating supplied by the advisor.",
            reproducibility="The independent language-model classifications are frozen per record in history_extract.py, hash-bound to the ordered canonical event IDs. Running the generator reproduces those judgments without another model call.",
            ordered_input_id_sha256=digest("\n".join(r["id"] for r in rows).encode()),
        ),
        rubric=dict(
            tracked_edit=-1,
            comment={
                "0": "positive assessment with no change requested",
                "-1": "local wording, clarification, presentation, or minimal editorial mark",
                "-2": "substantive explanation, method, validation detail, narrative restructuring, or missing contribution",
                "-3": "material challenge to factual accuracy, scope, attribution, numerical meaning, or interpretation",
            },
            mixed="Retain the negative weight and record the praise component when present. A factual challenge is evidence of the advisor's concern, even when Ben disagrees or the hypothesis remains unverified.",
            counting="Count each original message and each distinct original tracked action once. Do not split one comment into extra points. A deletion and insertion that relocate a sentence remain two distinct actions when the panel lists both. A single replacement containing both sides is one action.",
            missing_text="Seven October 29 edit cards visibly establish author, file, minute, operation, and distinct occurrence; a collapsed side is null, and displayed_fragments gives only the visible prefix. Comments that survive only as an ambiguous fragment remain unscored.",
            anchors="An exact selection, historical local revision line, raw UI offset, and nearby technical paragraph have different status labels. Earlier content is preserved literally. Intentionally missing or detached selections are null. A merged source or later scope note must not silently supply an original selection.",
            status="A resolution, rewrite, or later score is not a new advisor signal and cannot earn back any of the original points.",
            comparison="July/September/October event rescoring and the contemporaneous plus/minus teaching ledger, 0–4 readiness table, and weighted publication-audit scores have different units.",
            issue_links="Distinct actions in one discussion or on the same demonstrable issue can share an issue_id. Tags indicate transferable editorial patterns across dates. A common tag alone does not merge event identities.",
        ),
        tag_catalog={name: desc for name, desc in TAGS.values()},
        sources=SOURCES,
        signals=rows,
        unscored_evidence=list(grouped.values()),
        coverage=dict(
            early="The dated ajh26 inventory contains 66 entries from July 21 and September 10–11: 49 comment messages and 17 tracked actions. Two further September comments on dated hardware changes are documented in the October 15 report. Two September question marks are real distinct messages in both inventories, but the individual selections cannot be paired. Their selections remain in unscored_evidence.",
            october_15="The original CHATGPT_REPORT.md capture preserves 46 open threads, 15 pending tracked entries and 32 resolved threads across authors. Its C/T/R identities and the separate 61-entry action ledger are cross-references; older campaign weights and non-advisor material are excluded.",
            october_22_23="Multiple same-project source snapshots were reconciled. The final October 23 capture has 137 cards across authors and dates, including 63 advisor edit cards and 11 Ben edit cards. Each advisor message within a thread is a separate event; current anchors and historical passages can differ. Duplicate offsets at EOF are never used alone to identify a selection.",
            october_26="The October 26 SNS replies file corroborates an older advisor message; its new reply and tracked edits are by Ben. The torque-review exports and the committee report document activity by multiple people, with separate project scope.",
            october_27_28=dict(
                final_resolution_log=controls,
                earlier_log={"dissertation": 9, "torque_paper": 17},
                detail="The archived operations resolve older comments. Fifty-seven advisor thread/opening references have unique matches. An October 23 question mark maps to two possible advisor message IDs. A second authorless excerpt repeats a verified Ben message. The independent Review panel identifies the authors of earlier resolution descriptions where a unique match survives. The final record marks the motor schematic issue as still open. No automatic extra event or score follows a resolution.",
            ),
            october_29_30="The two chapter snapshots cover the then-visible reviews in Chapters 4–5. The committed report contains 42 advisor opening messages and two advisor follow-ups (10:13 and 10:56), corroborated in the raw chapter panel. The 28 advisor edit cards on October 29 include seven showing at least one collapsed side. The October 30 dissertation panel preserves three October 29 message IDs and zero then-pending tracked changes. The 19 user comments on a local PDF, the torque paper, and the defense slides are different evidence streams.",
            duplicates="Prefer stable message IDs. In their absence match verified author, original date, file, time when retained, wording, and within-capture occurrence. Whitespace and the specific broken display encoding of NBSP are normalized for comparison only; the raw sources preserve original bytes. Alternate October 14 and October 23 edit representations at the same file/minute are held out where version identity cannot be established.",
            side_products="Rulings, completed fixes, checklist actions, numerical readiness ratings, automated findings, committee-editor instructions, and Ben or Leo remarks are not advisor marks. No continuous platform action history or deleted/missing reviews have been reconstructed.",
            attributed_records_by_source=source_counts,
            excluded_per_capture=[
                dict(source_id=k[0], author=k[1], type=k[2], count=v) for k, v in sorted(EXCLUDED.items())
            ],
            source_commit_support="Git records have full commit:path and blob ID; the October 28/29 machine-local review snapshots have pinned independent SHA-256 values and embedded original bytes. The Git entry points and the October 29 transcript corroborate their scope.",
        ),
        counts=dict(
            all=tally(rows),
            by_month={
                m: tally([r for r in rows if r["date"].startswith(m)])
                for m in sorted(set(r["date"][:7] for r in rows))
            },
            by_date={
                d: tally([r for r in rows if r["date"] == d]) for d in sorted(set(r["date"] for r in rows))
            },
            tag_counts={k: sum(k in r["tags"] for r in rows) for k in [t[0] for t in TAGS.values()]},
            tag_counts_october={
                k: sum(k in r["tags"] and r["date"].startswith("2026-10") for r in rows)
                for k in [t[0] for t in TAGS.values()]
            },
            unscored_evidence_items=len(grouped),
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--inspect",
        action="store_true",
        help="print ordered classification inputs instead of rebuilding the JSON",
    )
    parser.add_argument(
        "--check", action="store_true", help="recompute the artifact and compare it with the existing file"
    )
    args = parser.parse_args()
    rows = extract()
    if args.inspect:
        inspection(rows)
        return
    result = make_output(rows)
    path = HERE / "HISTORICAL_SIGNALS.json"
    encoded = (json.dumps(result, ensure_ascii=False, indent=2) + "\n").encode()
    if args.check:
        if path.read_bytes() != encoded:
            raise SystemExit("Recomputed JSON differs")
    else:
        path.write_bytes(encoded)
    print(f"{len(rows)} unique events; {result['counts']['all']['score']} points; SHA-256 {digest(encoded)}")


if __name__ == "__main__":
    import sys

    sys.stdout.reconfigure(encoding="utf-8")
    main()
