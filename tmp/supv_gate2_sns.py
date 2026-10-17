# Supervisor Gate-2 scan: AI-slop patterns over the ADDED lines of this batch's diffs.
# Reads the two diff logs captured from git diff (60-sns, 93-AppendixB), extracts added
# lines, strips LaTeX comments, runs rubric Tier A + Tier B regexes + structural checks.
import re, sys

DIFF_FILES = [
    r"D:\Github\Bipedal_Robot\tmp\supv_sns_diff2.log",
]

TIER_A = [
    (r"\bdelv(?:e|es|ed|ing)\b", "delve"),
    (r"\bit is important to note\b", "it is important to note"),
    (r"\bit is worth noting that\b", "it is worth noting that"),
    (r"\bit should be noted that\b", "it should be noted that"),
    (r"\bplays? (?:a|an) (?:crucial|pivotal|vital|key) role\b", "plays a crucial/pivotal/vital/key role"),
    (r"\blandscape\b", "landscape"),
    (r"\btapestry\b", "tapestry"),
    (r"\b(?:a|the) testament to\b", "a/the testament to"),
    (r"\bcutting-edge\b", "cutting-edge"),
    (r"\bseamless(?:ly)?\b", "seamless(ly)"),
    (r"\bmyriad\b", "myriad"),
    (r"\bplethora\b", "plethora"),
]

TIER_B = [
    (r"(?m)(?:^|\.\s+|\\par\b)(?:Moreover|Furthermore|Additionally),", "connective sentence-opener"),
    (r"\bholistic\b", "holistic"),
    (r"\bleverage(?:s|d)?\b", "leverage as verb"),
    (r"\bunderscore(?:s|d)?\b", "underscore(s/d)"),
    (r"\bshowcase(?:s|d)?\b", "showcase(s/d)"),
    (r"\bparadigm\b", "paradigm"),
    (r"\bcomprehensive\b", "comprehensive"),
    (r"\bmultifaceted\b", "multifaceted"),
    (r"\bin the realm of\b", "in the realm of"),
    (r"\bnavigate the complexities\b", "navigate the complexities"),
    (r"\bnot only\b[^.]*\bbut also\b", "not only ... but also"),
    (r"\bas requested\b|\bthe committee suggested\b|\bcommittee requested\b", "reviewer-text leak"),
]

added = []
for df in DIFF_FILES:
    with open(df, encoding="utf-8", errors="replace") as f:
        for ln in f:
            if ln.startswith("+") and not ln.startswith("+++"):
                added.append(ln[1:].rstrip("\n"))

def strip_comment(s):
    # remove content after an unescaped % (LaTeX comment)
    out, i = [], 0
    while i < len(s):
        c = s[i]
        if c == "\\" and i + 1 < len(s):
            out.append(s[i:i+2]); i += 2; continue
        if c == "%":
            break
        out.append(c); i += 1
    return "".join(out)

prose = [strip_comment(a) for a in added]
text = "\n".join(prose)
words = len(re.findall(r"[A-Za-z][A-Za-z'-]*", text))
print(f"added lines scanned: {len(added)}; ~{words} words of prose (comments stripped)")

fails, flags = [], []
for pat, name in TIER_A:
    for m in re.finditer(pat, text, re.I):
        ctx = text[max(0, m.start()-60):m.end()+60].replace("\n", " ")
        fails.append(f"TIER A [{name}] :: ...{ctx}...")
for pat, name in TIER_B:
    for m in re.finditer(pat, text, re.I):
        ctx = text[max(0, m.start()-60):m.end()+60].replace("\n", " ")
        flags.append(f"TIER B [{name}] :: ...{ctx}...")

# Em-dash count in changed prose
emd = text.count("---")
print(f"em-dash '---' occurrences in changed prose: {emd} (~{words} words)")

# Connective openers per paragraph (added text only, paragraph = blank-line split)
paras = [p for p in re.split(r"\n\s*\n", text) if p.strip()]
openers = sum(1 for p in paras if re.match(r"\s*(Moreover|Furthermore|Additionally),", p))
print(f"sentence-opening connectives found: {openers}")

print("\n=== TIER A (hard fail) ===")
print("\n".join(fails) if fails else "NONE")
print("\n=== TIER B (judged) ===")
print("\n".join(flags) if flags else "NONE")
