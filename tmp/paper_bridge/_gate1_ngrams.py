# Gate 1 supervisor check (v2, per-line state machine = accurate line numbers):
# n-gram overlap of the new manuscript against
# (a) the published Actuators paper (bolen_2026) - hard gate, >=10-word verbatim
# (b) dissertation 30-torque.tex (bolen_isometric_2025 subject) - adaptation check
# (c) dissertation 20-muscles.tex - force-model prose context
import re, sys, io

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

ROOT = r"D:\GitHub\Bipedal_Robot"
MAIN = ROOT + r"\Documentation\Reports and Papers\Knee_Torque_Test\ProofFinal\Main.tex"
ACT  = ROOT + r"\Documentation\Reports and Papers\Isometric_Force_Test\Final\Final_proofread\manuscript.v5\actuators-4207821-edited.tex"
CH30 = ROOT + r"\Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\30-torque.tex"
CH20 = ROOT + r"\Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\20-muscles.tex"

DROP_ENVS = {"equation", "align", "tabular", "tabularx", "eqnarray", "adjustwidth"}

def clean_lines(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        raw = f.readlines()
    out = []           # list of (cleaned_text, lineno); None = skipped env line
    env_stack = []
    for lineno, line in enumerate(raw, 1):
        text = re.sub(r"(?<!\\)%.*", "", line)
        begins = re.findall(r"\\begin\{([a-zA-Z*]+)\}", text)
        # decide: if inside a drop-env, skip (but keep tracking nesting)
        if env_stack and env_stack[-1] in DROP_ENVS:
            for b in begins:
                env_stack.append(b)
            ends = re.findall(r"\\end\{([a-zA-Z*]+)\}", text)
            for e in ends:
                if env_stack and env_stack[-1] == e:
                    env_stack.pop()
            out.append((None, lineno))
            continue
        # not inside drop-env: strip begin/end markers of drop-envs from text
        for b in begins:
            env_stack.append(b)
        ends = re.findall(r"\\end\{([a-zA-Z*]+)\}", text)
        for e in ends:
            if env_stack and env_stack[-1] == e:
                env_stack.pop()
        if any(b in DROP_ENVS for b in begins):
            out.append((None, lineno))
            continue
        # per-line substitutions (line-local => line numbers stay exact)
        text = re.sub(r"\$[^$]*\$", " MATH ", text)
        text = re.sub(r"\\(includegraphics|bibliography|graphicspath|input|bibliographystyle|externaldocument)\s*(\[[^\]]*\])?\{[^{}]*\}", " ", text)
        text = re.sub(r"\\(citep|citealp|citet|cite|ref|eqref|autoref|nameref|url|href| hyperlink)\s*(\[[^\]]*\])*\{[^{}]*\}", " CITE ", text)
        text = re.sub(r"\\[a-zA-Z@]+\*?", " ", text)
        out.append((text, lineno))
    return out

def tokenize(cleaned):
    toks = []
    for text, lineno in cleaned:
        if text is None:
            continue
        for m in re.finditer(r"[A-Za-z0-9]+", text):
            toks.append((m.group(0).lower(), lineno))
    return toks

def build_index(toks, n):
    idx = {}
    words = [w for w, _ in toks]
    for i in range(len(words) - n + 1):
        idx.setdefault(tuple(words[i:i+n]), []).append(toks[i][1])
    return idx

def matched_runs(main_toks, src_idx, n):
    mw = [w for w, _ in main_toks]
    matched = [False] * len(mw)
    for i in range(len(mw) - n + 1):
        if tuple(mw[i:i+n]) in src_idx:
            for k in range(i, i + n):
                matched[k] = True
    runs, i = [], 0
    while i < len(mw):
        if matched[i]:
            j = i
            while j < len(mw) and matched[j]:
                j += 1
            key = tuple(mw[i:i+n])
            lines = sorted(set(ln for _, ln in main_toks[i:j]))
            runs.append((mw[i:j], lines, src_idx[key][0]))
            i = j
        else:
            i += 1
    return runs

def report(srcname, main_toks, src_toks, nlist, minwords):
    for n in nlist:
        idx = build_index(src_toks, n)
        rs = [r for r in matched_runs(main_toks, idx, n) if len(r[0]) >= max(n, minwords)]
        print(f"\n=== vs {srcname} :: window {n} : {len(rs)} matched run(s) of >= {max(n, minwords)} words ===")
        for w, lines, srcline in rs:
            txt = " ".join(w)
            if len(txt) > 260:
                txt = txt[:257] + "..."
            print(f"  Main L{lines[0]}-{lines[-1]} | src L{srcline} | {len(w)}w: {txt}")

main_toks = tokenize(clean_lines(MAIN))
print(f"Main.tex prose tokens: {len(main_toks)}")
for label, path, wins in [("PUBLISHED-bolen_2026", ACT, [10]),
                          ("DISS-30-torque", CH30, [20, 12]),
                          ("DISS-20-muscles", CH20, [12])]:
    src_toks = tokenize(clean_lines(path))
    print(f"\n##### SOURCE {label}: {len(src_toks)} tokens #####")
    report(label, main_toks, src_toks, wins, wins[0])
