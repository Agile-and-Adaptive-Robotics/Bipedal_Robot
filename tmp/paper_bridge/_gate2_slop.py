# Gate 2 supervisor check: AI-slop patterns over Main.tex prose (comments stripped),
# per supervisor-rubric.md Tier A / Tier B + structural checks + BRIEF-specific bans.
import re, sys, io

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

MAIN = r"D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Knee_Torque_Test\ProofFinal\Main.tex"

with open(MAIN, encoding="utf-8", errors="replace") as f:
    raw = f.readlines()

prose = []  # (text, lineno) with comments stripped
for lineno, line in enumerate(raw, 1):
    t = re.sub(r"(?<!\\)%.*", "", line)
    prose.append((t, lineno))

TIER_A = [
    (r"\bdelv(e|es|ed|ing)\b", "delve"),
    (r"\bit is important to note\b", "it-is-important-to-note"),
    (r"\bit is worth noting that\b", "it-is-worth-noting"),
    (r"\bit should be noted that\b", "it-should-be-noted"),
    (r"\bplays? (a|an) (crucial|pivotal|vital|key) role\b", "plays-a-X-role"),
    (r"\blandscape\b", "landscape"),
    (r"\btapestry\b", "tapestry"),
    (r"\b(a|the) testament to\b", "testament-to"),
    (r"\bcutting-edge\b", "cutting-edge"),
    (r"\bseamless(ly)?\b", "seamless(ly)"),
    (r"\bmyriad\b", "myriad"),
    (r"\bplethora\b", "plethora"),
]
TIER_B = [
    (r"(?<![A-Za-z])(Moreover|Furthermore|Additionally),", "connective-opener"),
    (r"\bholistic\b", "holistic"),
    (r"\bleverag(e|es|ed|ing)\b", "leverage"),
    (r"\bunderscore(s|d)?\b", "underscore"),
    (r"\bshowcase(s|d)?\b", "showcase"),
    (r"\bparadigm\b", "paradigm"),
    (r"\bcomprehensive\b", "comprehensive"),
    (r"\bmultifaceted\b", "multifaceted"),
    (r"\bin the realm of\b", "in-the-realm-of"),
    (r"\bnavigate the complexities\b", "navigate-complexities"),
    (r"\bnot only\b", "not-only-x-but-also"),
    (r"\bcrucial\b", "crucial"),
    (r"\bpivotal\b", "pivotal"),
    (r"\brobust\b", "robust"),
    (r"\bnotably\b", "notably"),
    (r"\bremarkabl(e|y)\b", "remarkable"),
    (r"\bgroundbreaking\b", "groundbreaking"),
    (r"\bnovel\b", "novel"),
]

def scan(patterns, label):
    print(f"\n===== {label} =====")
    hits = 0
    for pat, name in patterns:
        for text, lineno in prose:
            for m in re.finditer(pat, text, re.I):
                s = max(0, m.start() - 45); e = min(len(text), m.end() + 45)
                ctx = text[s:e].replace("\n", " ").strip()
                print(f"  L{lineno} [{name}]: ...{ctx}...")
                hits += 1
    if hits == 0:
        print("  (no hits)")

scan(TIER_A, "TIER A (hard fail)")
scan(TIER_B, "TIER B (judged)")

# BRIEF-specific bans
print("\n===== BRIEF-SPECIFIC =====")
for pat, name in [(r"\bplants?\b", "word 'plant'"),
                  (r"\bcampaign\b", "campaign"),
                  (r"\bOverleaf\b", "overleaf"),
                  (r"\bZcode\b|\bagent\b|\bsession\b", "tooling-leak"),
                  (r"\d(mm|cm|kPa|kN|N\W)", "bare-number-unit (non-siunitx)"),
                  (r"20(2[0-9])[-/](0?[0-9]|1[0-2])[-/]", "ISO-style date"),
                  (r"\b2025\b|\b2026\b|\b2027\b", "year token"),
                  (r"\.(mat|m|csv|log|xlsx|txt|fig)\b", "file extension"),
                  (r"\bDig_out\b|\bOpt_run\b|minimize\w*|tmp/", "artifact/script name"),
                  (r"\brow 77\b|\bpick\b|\btrial \d+\b|\bwinner\b", "pick/trial language"),
                  (r"as requested|the committee|previously unclear|we have revised", "reviewer-leak")]:
    n = 0
    for text, lineno in prose:
        for m in re.finditer(pat, text, re.I):
            s = max(0, m.start() - 50); e = min(len(text), m.end() + 50)
            ctx = text[s:e].replace("\n", " ").strip()
            print(f"  L{lineno} [{name}]: ...{ctx}...")
            n += 1
    if n == 0:
        print(f"  [{name}] (no hits)")

# em-dash padding: count '---' in prose
em = sum(t.count("---") for t, _ in prose)
words = sum(len(re.findall(r"[A-Za-z]+", t)) for t, _ in prose)
print(f"\n'em-dashes (---)': {em} in ~{words} prose words (limit ~1 per 500)")
en = sum(len(re.findall(r"(?<![-\d])--(?![-\d])", t)) for t, _ in prose)
print(f"'en-dash --' occurrences: {en} (ranges, exempt)")

# paragraph rhythm survey (rough): paragraphs = blank-line separated blocks in body
body = "".join(t for t, _ in prose)
paras = [p.strip() for p in body.split("\n\n") if len(p.strip()) > 150]
print(f"\nparagraphs >150 chars: {len(paras)}")
