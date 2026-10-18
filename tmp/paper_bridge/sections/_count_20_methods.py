import re
frag = open(r'D:/GitHub/Bipedal_Robot/tmp/paper_bridge/sections/20_methods.tex', encoding='utf-8').read()
lines = [l for l in frag.splitlines() if not l.lstrip().startswith('%')]
text = '\n'.join(lines)
text = re.sub(r'\\begin\{equation\}.*?\\end\{equation\}', ' ', text, flags=re.S)
text = re.sub(r'\\begin\{figure\}.*?\\end\{figure\}', ' ', text, flags=re.S)
text = re.sub(r'\\cite[tp]?\{[^}]*\}', ' CITATION ', text)
text = re.sub(r'\\(eq)?ref\{[^}]*\}', ' REF ', text)
text = re.sub(r'\\label\{[^}]*\}', ' ', text)
text = re.sub(r'\\qty\{([^}]*)\}\{[^}]*\}', r' \1 ', text)
text = re.sub(r'\\qtylist\{[^}]*\}\{[^}]*\}', ' NUM ', text)
text = re.sub(r'\\qtyrange\{[^}]*\}\{[^}]*\}\{[^}]*\}', ' RANGE ', text)
text = re.sub(r'\\si\{[^}]*\}', ' UNIT ', text)
text = re.sub(r'\\ang\{[^}]*\}', ' ANG ', text)
text = re.sub(r'\$[^$]*\$', ' MATH ', text)
text = re.sub(r'\\[a-zA-Z]+\*?(\[[^\]]*\])?', ' ', text)
text = re.sub(r'[{}~$%&]', ' ', text)
skip = ('CITATION', 'REF', 'NUM', 'RANGE', 'UNIT', 'ANG', 'MATH')
words = [w for w in re.split(r'\s+', text) if re.search(r'[A-Za-z]', w) and w not in skip]
print('accurate prose words (excl. equations, figure floats, math, citation keys):', len(words))

# per-subsection
parts = re.split(r'(\\subsection\{[^}]*\})', '\n'.join(lines))
import sys
buf = []
for p in parts:
    if p.startswith('\\subsection'):
        if buf:
            print('  ', buf[0][:55], '->', len(buf[1]))
        buf = [p, '']
    elif buf:
        buf[1] += p
if buf: print('  ', buf[0][:55], '->', len(buf[1]))
