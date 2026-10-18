"""Minimal Word-RTF to text extractor for PlotPrompt.rtf (no deps).

Handles: \par -> newline, \'hh hex escapes, \{ \} \\ escapes, control words
skipped, {\fonttbl...}/{\colortbl...}/{\stylesheet...}/{\info...} groups
dropped wholesale, \*\... destination groups skipped.
"""
import re
import sys

SRC = r"D:\GitHub\Bipedal_Robot\Testing_Data\2022_02_Festo\PlotPrompt.rtf"
DST = r"D:\GitHub\Bipedal_Robot\tmp\paper_bridge\PlotPrompt_rtf.txt"

raw = open(SRC, "rb").read().decode("cp1252", errors="replace")

SKIP_GROUPS = ("fonttbl", "colortbl", "stylesheet", "info", "generator",
               "pict", "object", "header", "footer", "themedata", "datastore")

out = []
i = 0
skip_stack = []  # group-start indices where skipping began
depth = 0

def ctrl_word(s, j):
    m = re.match(r"\\([a-zA-Z]+)(-?\d+)? ?", s[j:])
    return m

while i < len(raw):
    ch = raw[i]
    if ch == "{":
        depth += 1
        # peek for destination marker \fonttbl etc. or \*
        m = ctrl_word(raw, i + 1)
        star = raw.startswith(r"\*", i + 1)
        if (m and m.group(1) in SKIP_GROUPS) or star:
            skip_stack.append(depth)
        i += 1
    elif ch == "}":
        if skip_stack and skip_stack[-1] == depth:
            skip_stack.pop()
        depth -= 1
        i += 1
    elif skip_stack:
        i += 1
    elif ch == "\\":
        if i + 1 < len(raw) and raw[i + 1] == "'":
            hh = raw[i + 2:i + 4]
            try:
                out.append(chr(int(hh, 16)))
            except ValueError:
                pass
            i += 4
            continue
        m = ctrl_word(raw, i)
        if m:
            word = m.group(1)
            if word == "par":
                out.append("\n")
            elif word in ("tab",):
                out.append("\t")
            i += m.end()
            continue
        if i + 1 < len(raw) and raw[i + 1] in "{}\\":
            out.append(raw[i + 1])
            i += 2
            continue
        i += 1
    elif ch in "\r\n":
        i += 1
    else:
        out.append(ch)
        i += 1

text = "".join(out)
text = re.sub(r"\n{3,}", "\n\n", text)
open(DST, "w", encoding="utf-8", newline="\n").write(text)
print("WROTE", DST, "chars", len(text))
