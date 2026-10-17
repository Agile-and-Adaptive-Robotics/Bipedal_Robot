import re

p = r'Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\50-database.tex'
src = open(p, encoding='utf-8').read()
lines = src.splitlines()

print('lines', len(lines), 'words', len(src.split()))

# comment check: any line containing % that is not \% escape
for i, l in enumerate(lines, 1):
    if re.search(r'(?<!\\)%', l):
        print('COMMENT at', i, repr(l))

# brace balance (ignoring escaped \{ \})
body = re.sub(r'\\[{}]', '', src)
print('brace balance', body.count('{') - body.count('}'))

# begin/end environment balance
print('begin/end', len(re.findall(r'\\begin\{', src)), len(re.findall(r'\\end\{', src)))

# em-dash count
print('em-dash count', src.count('---'))

tier_a = [r'delv(e|es|ed|ing)', r'it is important to note', r'it is worth noting',
          r'it should be noted', r'plays? (a|an) (crucial|pivotal|vital|key) role',
          r'landscape', r'tapestry', r'(a|the) testament to', r'cutting-edge',
          r'seamless(ly)?', r'myriad', r'plethora']
tier_b = [r'\bmoreover\b', r'\bfurthermore\b', r'\badditionally\b', r'\bholistic\b',
          r'\bleverage[ds]?\b', r'\bunderscore[ds]?\b', r'\bshowcase[ds]?\b',
          r'\bparadigm\b', r'\bcomprehensive\b', r'\bmultifaceted\b',
          r'in the realm of', r'navigate the complexities',
          r'not only .{0,80}but also', r'\bplant\b']
for name, pats in [('A', tier_a), ('B', tier_b)]:
    for pat in pats:
        for m in re.finditer(pat, src, re.I):
            print('TIER', name, 'HIT', repr(m.group(0)))

brit = [r'behaviour', r'fibre', r'modelled', r'colour', r'\bcentre\b']
for pat in brit:
    for m in re.finditer(pat, src, re.I):
        print('BRITISH HIT', repr(m.group(0)))

# citation keys used
print('cite keys', sorted(set(re.findall(r'\\cite[tp]\{([^}]*)\}', src))))
