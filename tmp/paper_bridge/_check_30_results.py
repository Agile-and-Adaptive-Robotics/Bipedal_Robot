import re
src = open('30_results.tex', encoding='utf-8').read()
body = re.sub(r'(?m)%.*$', '', src)
body_nofloat = re.sub(r'\\begin\{figure\}.*?\\end\{figure\}', '', body, flags=re.S)
body_words = len(re.findall(r"[A-Za-z0-9']+", re.sub(r'\\[a-zA-Z]+|[${}~]', ' ', body_nofloat)))
print("body prose words (excl. floats/comments):", body_words)
print("em-dashes:", src.count('---'))
print("includegraphics:", re.findall(r'includegraphics\[[^]]*\]\{([^}]*)\}', src))
print("labels:", re.findall(r'\\label\{([^}]*)\}', src))
print("refs:", sorted(set(re.findall(r'\\ref\{([^}]*)\}', src))))
print("qty:", len(re.findall(r'\\qty\{', src)), "| ang:", len(re.findall(r'\\ang\{', src)),
      "| qtylist:", len(re.findall(r'\\qtylist\{', src)))
print("diameter uses:", src.count(r'\diameter'))
