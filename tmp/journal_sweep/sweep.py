# -*- coding: utf-8 -*-
"""Journal-register mechanical sweep over the ten revised dissertation chapters.
Greps current text whole-file; 'NEW' checks are diff-aware vs tmp/diss_baseline_20261017
using whitespace-stripped whole-text membership (robust to rewrapping)."""
import re, os, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

ROOT = r"D:\Github\Bipedal_Robot"
CH = os.path.join(ROOT, r"Documentation\Reports and Papers\Dissertation\ProofFinal\chapters")
BASE = os.path.join(ROOT, r"tmp\diss_baseline_20261017")
CHAPTERS = ["10-introduction","15-background","20-muscles","30-torque","50-database",
            "60-sns","70-overall","80-conclusion","93-AppendixB","94-AppendixC"]

def load(p):
    for enc in ("utf-8","latin-1"):
        try:
            with open(p, encoding=enc) as f: return f.read(), enc
        except UnicodeDecodeError: continue
    return open(p, encoding='utf-8', errors='replace').read(), 'utf-8?'

def squash(s): return re.sub(r'\s+', '', s)

def lines_with(pat, text, flags=0):
    out = []
    for i, ln in enumerate(text.splitlines(), 1):
        if re.search(pat, ln, flags): out.append((i, ln.rstrip()))
    return out

UNIT_RE = re.compile(r'(?<![\w\\])(\d+(?:\.\d+)?)\s*(~|\\,|\s)?\s*(mm|kPa|MPa|Hz|kW|kN|N\\cdot|N\\!\\cdot|N\\,\\cdot|m/s|mm\^2|mm\^\{2\}|Nm\b|N\b|s\b|kg|degrees|deg)\b')
SI_CMD = re.compile(r'\\(qty|qtylist|SI|SIrange|si|num|ang|qtyrange)\s*(\[[^\]]*\])?\{')

def mask_si(text):
    """Replace contents of siunitx commands so unit grep doesn't hit them."""
    out = list(text)
    for m in SI_CMD.finditer(text):
        # crude brace matcher
        i = m.end()-1; depth=0; j=i
        while j < len(text):
            if text[j]=='{': depth+=1
            elif text[j]=='}':
                depth-=1
                if depth==0: break
            j+=1
        for k in range(i+1, j):
            if out[k].isalnum() or out[k] in '.;,/\\{}^~ ' or out[k].isdigit(): out[k] = ' '
    # also mask a second brace group (the unit group)
    return ''.join(out), ''.join(out)

def find_second_group(text):
    """Mask the 2nd brace group of qty/SI commands (units)."""
    res = list(text)
    for m in SI_CMD.finditer(text):
        # find end of first group
        i = m.end()-1; depth=0; j=i
        while j < len(text):
            if text[j]=='{': depth+=1
            elif text[j]=='}':
                depth-=1
                if depth==0: break
            j+=1
        # second group starts at j+1 if '{'
        if j+1 < len(text) and text[j+1]=='{':
            k=j+2; d=1
            while k < len(text):
                if text[k]=='{': d+=1
                elif text[k]=='}':
                    d-=1
                    if d==0: break
                k+=1
            for x in range(j+2,k):
                if res[x].isalnum() or res[x] in '.;,/\\{}^~ ' : res[x]=' '
    return ''.join(res)

for name in CHAPTERS:
    cur, enc = load(os.path.join(CH, name + ".tex"))
    try:
        base, _ = load(os.path.join(BASE, name + ".tex"))
        have_base = True
    except OSError:
        base = ""; have_base = False
    sq_cur, sq_base = squash(cur), squash(base)
    print("="*100)
    print(f"CHAPTER {name}  ({len(cur.splitlines())} lines, enc={enc}, baseline={'yes' if have_base else 'NO'})")

    # 1. campaign
    hits = lines_with(r'campaign', cur, re.I)
    print(f"-- 'campaign': {len(hits)} hits")
    for i,ln in hits: print(f"   L{i}: {ln.strip()[:160]}")

    # 2. of record
    hits = lines_with(r'\bof record\b', cur, re.I)
    print(f"-- 'of record': {len(hits)} hits")
    for i,ln in hits: print(f"   L{i}: {ln.strip()[:160]}")

    # 3. plant(s)
    hits = [(i,ln) for i,ln in lines_with(r'plant', cur, re.I)
            if not re.search(r'transplant|plantar|plantigrade', ln, re.I)]
    print(f"-- 'plant(s)' excluding transplant/plantar*: {len(hits)} hits")
    for i,ln in hits: print(f"   L{i}: {ln.strip()[:160]}")

    # 4. \phi
    hits = lines_with(r'\\phi', cur)
    print(f"-- '\\phi': {len(hits)} hits")
    for i,ln in hits: print(f"   L{i}: {ln.strip()[:160]}")

    # 5. bare units (siunitx masked)
    masked = find_second_group(mask_si(cur)[0])
    unit_hits = []
    for i, ln in enumerate(masked.splitlines(), 1):
        for m in UNIT_RE.finditer(ln):
            unit_hits.append((i, m.group(0), cur.splitlines()[i-1].strip()[:150]))
    print(f"-- bare-unit candidates outside siunitx: {len(unit_hits)}")
    for i,g,src in unit_hits[:40]: print(f"   L{i} [{g}]: {src}")

    # 6. SNN
    hits = lines_with(r'\bSNN', cur)
    print(f"-- 'SNN': {len(hits)} hits")
    for i,ln in hits: print(f"   L{i}: {ln.strip()[:160]}")

    # 7. pool
    hits = lines_with(r'\bpool', cur, re.I)
    print(f"-- 'pool*': {len(hits)} hits")
    for i,ln in hits: print(f"   L{i}: {ln.strip()[:160]}")

    # 8. NEW em-dash: occurrences of '---' in current whose 60-char squashed window is not in baseline
    new_em = []
    for m in re.finditer(r'---', cur):
        a, b = max(0, m.start()-50), min(len(sq_cur), 0)  # placeholder
    # do on squashed text
    for m in re.finditer(r'---', sq_cur):
        w = sq_cur[max(0,m.start()-45):m.start()+45]
        if w not in sq_base:
            new_em.append(w)
    print(f"-- '---' total in current: {len(re.findall(r'---', cur))}; NEW (not in baseline): {len(new_em)}")
    for w in new_em[:12]: print(f"   window: ...{w}...")

    # 9. slop phrases on current text, marked [NEW] if not in baseline
    slops = [r'\bFirstly\b', r'\bca\.', r'[Ii]t can be seen', r'[Ii]t should be noted',
             r'\bwe call\b', r'\bthe point is\b', r'\bthe distinction matters\b',
             r'\bdelve\b', r'important to note', r'crucial role', r'\blandscape\b',
             r'\bmoreover\b', r'\bfurthermore\b', r'\badditionally\b', r'\bnotably\b']
    print("-- slop phrases:")
    any_slop = False
    for pat in slops:
        for i, ln in lines_with(pat, cur, re.I):
            newsq = ' [NEW]' if squash(ln.strip(30*'%')) not in sq_base else ''
            # better new-check: window around match
            m2 = re.search(pat, ln, re.I)
            frag = squash(ln[max(0,m2.start()-40):m2.end()+40])
            newsq = ' [NEW]' if frag not in sq_base else ' [old]'
            print(f"   {pat} L{i}{newsq}: {ln.strip()[:150]}")
            any_slop = True
    if not any_slop: print("   none")

    # 10. doubled words (whole file, current)
    dbl = []
    for i, ln in enumerate(cur.splitlines(), 1):
        if ln.strip().startswith('%'): continue
        for m in re.finditer(r'\b([A-Za-z]{2,})\s+\1\b', ln, re.I):
            dbl.append((i, m.group(0), ln.strip()[:140]))
    print(f"-- doubled words: {len(dbl)}")
    for i,g,src in dbl: print(f"   L{i} [{g}]: {src}")

    # 11. \qtylist with single value
    qsingle = []
    for i, ln in enumerate(cur.splitlines(), 1):
        for m in re.finditer(r'\\qtylist\s*(?:\[[^\]]*\])?\{([^{}]*)\}', ln):
            if ';' not in m.group(1): qsingle.append((i, m.group(0), ln.strip()[:140]))
    print(f"-- \\qtylist single-value: {len(qsingle)}")
    for i,g,src in qsingle: print(f"   L{i} [{g}]: {src}")

    # 12. rhetorical questions (prose lines ending in ?) - whole file + NEW flag
    qs = [(i,ln) for i,ln in lines_with(r'\?', cur) if not ln.strip().startswith('%')]
    print(f"-- '?' lines (judge context): {len(qs)}")
    for i,ln in qs: print(f"   L{i}: {ln.strip()[:160]}")
print("="*100)
print("SWEEP DONE")
