"""Extract the MAIN <script> block (the last one; the embedded-template
script precedes it) to a .js file for `node --check`."""
import sys

src = open(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
           r"\connectome_block_editor.html", encoding="utf-8").read()
js = src[src.rindex("<script>") + 8:src.rindex("</script>")]
out = r"D:\Github\Bipedal_Robot\tmp\editor_main_20261003.js"
open(out, "w", encoding="utf-8").write(js)
print("wrote", out, len(js), "chars")
