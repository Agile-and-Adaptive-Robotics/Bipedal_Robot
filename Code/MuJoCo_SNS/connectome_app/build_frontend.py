"""Package the complete original editor and apply desktop extensions."""
from pathlib import Path
import hashlib,json
HERE=Path(__file__).parent
ORIGINAL=HERE.parent/'spinal/connectome_block_editor.html'
def build(destination):
    source=ORIGINAL.read_text(encoding='utf-8')
    # Disable browser-session restoration only in the generated desktop page.
    a=source.index('try {\n  const saved = JSON.parse(localStorage.getItem("cme_models_v1"));')
    b=source.index('\nfunction cur()',a)
    source=source[:a]+source[b:]
    a=source.index('function persist() {');b=source.index('\nfunction renderTabs()',a)
    source=source[:a]+"function persist() { if (window.desktopChanged) window.desktopChanged(); }\n"+source[b:]
    source=source.replace('Connectome Block Editor v3.3','Connectome Editor Desktop')
    source=source.replace('const TYPES = {','const TYPES = {\n  ADAPTER:{shape:"circle",color:"#FA8E80",w:2,grp:"Physical maps"},\n  SPINDLE:{shape:"rect",color:"#C7AA78",w:2,grp:"Physical maps"},\n  TENDON:{shape:"rect",color:"#C7AA78",w:2,grp:"Physical maps"},')
    source=source.replace('Math.hypot(hp.x - n.x, hp.y - n.y) < 24','desktopNodeHit(n,hp)')
    # The desktop drag wrapper snapshots before movement; do not add a second
    # post-movement snapshot, which would consume the first Undo without moving.
    source=source.replace('if (drag && drag.moved) {\n    snapshot();','if (drag && drag.moved) {')
    source=source.replace('<button id="bSave">Export JSON</button>','<button id="bSave">Save</button><button id="dSaveAs">Save As…</button><button id="dSubSave">Save Subnetwork As…</button><button id="dBuild">Build Neural Model</button><button id="dSvg">Export SVG</button>')
    source=source.replace('Load JSON from file…','Open…')
    source=source.replace('<div id="tabs"></div>','<div id="desktopPath" style="font-size:12px;padding:6px 0;color:#42566b">Unsaved project</div><div id="tabs"></div>')
    source=source.replace('tabs auto-saved','explicit project Save / Save As')
    source=source.replace('different models, auto-saved in this browser; templates load into a NEW','different models; save changes to project files. Templates load into a NEW')
    source=source.replace('// ---------------- boot ----------------','// ---------------- desktop extensions ----------------\n'+(HERE/'desktop_editor.js').read_text(encoding='utf-8')+'\n// ---------------- boot ----------------')
    source=source.replace('</head>','<script src="qrc:///qtwebchannel/qwebchannel.js"></script></head>')
    # Run bridge initialization after the original boot, retaining every command.
    source=source.replace('</body>','<script>desktopInitialize();</script></body>')
    path=Path(destination);path.parent.mkdir(parents=True,exist_ok=True);path.write_text(source,encoding='utf-8')
    return path
