import re
import sys

from measure_words import strip_comments, word_count, find_matching_brace


def section_profile(path):
    raw = open(path, encoding="utf-8-sig").read()
    text = strip_comments(raw)
    # mask captions and table bodies so narrative is what's counted
    caption_spans = []
    for m in re.finditer(r"\\caption(?:\[[^\]]*\])?\s*\{", text):
        c = find_matching_brace(text, m.end() - 1)
        if c > 0:
            caption_spans.append((m.start(), c))
    masked = list(text)
    for a, b in caption_spans:
        for i in range(a, b + 1):
            masked[i] = " "
    s = "".join(masked)
    env_spans = []
    for m in re.finditer(r"\\begin\{(tabularx?|sidewaystable)\*?\}", s):
        env = m.group(1)
        e = s.find("\\end{%s}" % env, m.end())
        if e > 0:
            env_spans.append((m.start(), e + len("\\end{%s}" % env)))
    masked2 = list(s)
    for a, b in env_spans:
        for i in range(a, b + 1):
            masked2[i] = " "
    s2 = "".join(masked2)

    # split at section/subsection/chapter boundaries
    marks = [(m.start(), m.group(0)) for m in re.finditer(r"\\(sub)?section\*?\{[^}]*\}|\\chapter\{[^}]*\}", s2)]
    marks.append((len(s2), "END"))
    for i in range(len(marks) - 1):
        head = re.sub(r"\\(sub)?section\{?", "", marks[i][1]).rstrip("}").strip()
        chunk = s2[marks[i][0] : marks[i + 1][0]]
        n = word_count(chunk)
        print(f"{n:>6}  {head[:80]}")


if __name__ == "__main__":
    section_profile(sys.argv[1])
