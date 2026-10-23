# Corpus-vocabulary scan for watch_config.json keyword seeding (one-off, 2026-10-22).
import json, collections, re

d = json.load(open(r"D:\Github\Bipedal_Robot\SADb_audit\export\sadb_export.json", encoding="utf-8"))
papers = [r for r in d if r.get("recordType") == "paper"]
text = " ".join((r.get("title") or "") + " " + (r.get("notes") or "") for r in papers).lower()
text = re.sub(r"[^a-z\s-]", " ", text)
toks = text.split()
uni = collections.Counter(toks)
bi = collections.Counter(" ".join(p) for p in zip(toks, toks[1:]))
tri = collections.Counter(" ".join(p) for p in zip(toks, toks[1:], toks[2:]))

STOP = set("the a an of and or in on to for with from by at as is are was were be been its it this that these those we their his her during after before between via both all not no than then thus also can may might could would should using used use show shows showed suggest suggests role roles study studies new novel human humans rat rats cat cats mouse mice".split())

print("TOP UNIGRAMS (content):")
print(", ".join("%s(%d)" % (w, n) for w, n in uni.most_common(60) if w not in STOP and len(w) > 3))
print("\nTOP BIGRAMS:")
print("; ".join("%s (%d)" % (w, n) for w, n in bi.most_common(40) if not (w.split()[0] in STOP or w.split()[1] in STOP)))
print("\nTOP TRIGRAMS:")
print("; ".join("%s (%d)" % (w, n) for w, n in tri.most_common(25)))
