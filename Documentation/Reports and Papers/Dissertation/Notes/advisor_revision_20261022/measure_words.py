#!/usr/bin/env python3
"""Measure narrative and caption word counts for dissertation chapters.

Counts words in LaTeX source, splitting:
  - narrative: all text outside \caption{...} blocks (comments stripped)
  - captions: text inside \caption[short]{long} long blocks only
Tables' body text is excluded from narrative (counted separately as table words)
so that trimming prose is measured against prose.

Usage: python measure_words.py file1.tex [file2.tex ...]
Prints: file, narrative_words, caption_words, table_words, total_lines
"""
import re
import sys


def strip_comments(text: str) -> str:
    out = []
    for line in text.split("\n"):
        # remove whole-line and trailing comments (naive but safe: no % in \url here)
        i = line.find("%")
        if i >= 0:
            line = line[:i]
        out.append(line)
    return "\n".join(out)


def find_matching_brace(text: str, open_idx: int) -> int:
    depth = 0
    i = open_idx
    while i < len(text):
        c = text[i]
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                return i
        i += 1
    return -1


def word_count(s: str) -> int:
    # drop latex commands (backslash words), keep their args' words in the stream
    s = re.sub(r"\\[a-zA-Z]+\*?", " ", s)
    s = re.sub(r"[{}~]", " ", s)
    words = [w for w in re.split(r"\s+", s) if re.search(r"[A-Za-z0-9]", w)]
    return len(words)


def measure(path: str):
    with open(path, "r", encoding="utf-8-sig", errors="replace") as f:
        raw = f.read()
    text = strip_comments(raw)
    n_lines = len(raw.split("\n"))

    # extract caption long blocks: \caption[...]{...} or \caption{...}
    caption_words = 0
    caption_spans = []
    for m in re.finditer(r"\\caption(?:\[[^\]]*\])?\s*\{", text):
        open_idx = m.end() - 1
        close_idx = find_matching_brace(text, open_idx)
        if close_idx < 0:
            continue
        caption_words += word_count(text[open_idx + 1 : close_idx])
        caption_spans.append((m.start(), close_idx))

    # mask captions out of narrative
    masked = list(text)
    for a, b in caption_spans:
        for i in range(a, b + 1):
            masked[i] = " "
    masked_s = "".join(masked)

    # extract tabular/tabularx bodies as table words
    table_words = 0
    env_spans = []
    for m in re.finditer(r"\\begin\{(tabularx?|sidewaystable)\*?\}", masked_s):
        # for sidewaystable, take through \end{sidewaystable}; for tabular through \end{tabular}
        env = m.group(1)
        end_marker = "\\end{%s}" % env
        e = masked_s.find(end_marker, m.end())
        if e < 0:
            continue
        # sidewaystable wraps a caption (already masked) + tabularx
        chunk = masked_s[m.end() : e]
        table_words += word_count(chunk)
        env_spans.append((m.start(), e + len(end_marker)))
    # second pass: mask table bodies out of narrative
    masked2 = list(masked_s)
    for a, b in env_spans:
        for i in range(a, b + 1):
            masked2[i] = " "
    narrative_s = "".join(masked2)

    # also mask figure includegraphics blocks? includegraphics has no prose. fine.
    return word_count(narrative_s), caption_words, table_words, n_lines


if __name__ == "__main__":
    print(f"{'file':<24}{'narrative':>10}{'captions':>10}{'tables':>9}{'lines':>7}")
    for p in sys.argv[1:]:
        n, c, t, l = measure(p)
        name = p.split("\\")[-1].split("/")[-1]
        print(f"{name:<24}{n:>10}{c:>10}{t:>9}{l:>7}")
