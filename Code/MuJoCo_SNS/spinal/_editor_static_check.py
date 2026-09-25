"""Static integrity check for the block editor: every getElementById
target must exist as id= in the HTML; every onclick-wired id must be
unique. Catches the duplicate-button / missing-palette class of bug."""
import re

src = open("connectome_block_editor.html", encoding="utf-8").read()
# NB: the embedded-template <script> (2026-09-24) precedes the main
# one — slice the LAST script block (index() would grab the embed's
# closer and hand back an empty js, making every check vacuously pass).
js = src[src.rindex("<script>") + 8:src.rindex("</script>")]
ids_in_html = set(re.findall(r'id="([^"]+)"', src))
get_ids = set(re.findall(r"getElementById\(\"([^\"]+)\"\)", js))
missing = get_ids - ids_in_html
print("missing ids:", missing if missing else "NONE")

# duplicates
import collections
all_ids = re.findall(r'id="([^"]+)"', src)
dupes = [i for i, n in collections.Counter(all_ids).items() if n > 1]
print("duplicate ids:", dupes if dupes else "NONE")

# functions called vs defined
defined = set(re.findall(r"function (\w+)", js)) | \
    set(re.findall(r"const (\w+) = ", js))
called = set(re.findall(r"(?:addEventListener|\.)?(\w+)\(\)", js))
# filter obvious non-functions
known_non = {"use", "strict"}
undef = {c for c in called if c not in defined and c not in known_non
         and not c[0].isupper()}
print("possibly-undefined calls:", undef if undef else "NONE")
print("palette build present:", "palette" in js and "pitem" in js)
print("templates present:", "const TEMPLATES" in js)
