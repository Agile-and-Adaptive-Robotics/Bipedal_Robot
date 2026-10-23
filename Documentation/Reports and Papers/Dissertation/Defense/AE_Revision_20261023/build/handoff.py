from pathlib import Path
import json,hashlib,zipfile
from lxml import etree
root=Path(__file__).resolve().parents[1]
specs=json.loads((root/'slide_plan.json').read_text(encoding='utf-8'))
assets=json.loads((root/'build/asset_manifest.json').read_text(encoding='utf-8'))
lines=['# Slide review and figure sources','', 'Main talk: slides 1–31. Hidden backup: slides 32–37. Full explanatory notes are embedded in the deck.','']
for item in specs:
    lines += [f"## {item['slide']}. {item['title']}",'',f"{item['section']}. {'Hidden backup.' if item['backup'] else 'Main talk.'}",f"Claim source: {item['cite']}."]
    for k in item['assets']:
        m=assets[k];lines.append(f"- Figure `{k}`: `{m['source']}`.")
        if 'crop' in m:lines.append(f"  Selected view: {m['crop']}.")
    if not item['assets']:lines.append('Native editable text, table, chart, or conceptual diagram. Numeric chart values come from the cited source prose, not image digitization.')
    lines += ['']
lines += ['# Original-slide dispositions','', 'Original slide numbers refer to the preserved 54-slide deck. Removal here means exclusion from the revised copy, not deletion of source files.','', '| Original | Revised decision |','|---|---|']
decisions={1:'Rebuilt title (1)',2:'Replace dense objectives with gap mapping (2)',3:'Refocus scope and attribution (3)',4:'Integrate motivation into opening notes',5:'Explain BPA in characterization notes; remove separate descriptive slide',6:'Merge with revised gap map (2)',7:'Enlarge characterization figure (4)',8:'Retain measured length evidence (5)',9:'Move equations and coefficients into notes/backup (32)',10:'Editable reported-value comparison (6)',11:'Retain normalized-curve evidence (7)',12:'Move fit/split details into notes',13:'Select representative existing comparison panels (8)',14:'Move force-model implementation detail to backup (32)',15:'Refocus apparatus purpose (9)',16:'Retain scoped hybrid diagnosis (10)',17:'Refocus correction mechanisms (11)',18:'Separate identification from redesign (12)',19:'Retain representative corrected traces (13)',20:'Combine flexor/extensor validation means (14)',21:'Move wrapping interpretation to model notes',22:'Move parameter nonuniqueness and provenance into notes',23:'Frame earlier shortfall as redesign motivation (15)',24:'Move earlier configuration to backup (33)',25:'Replace ladder with specific route/hardware evidence (16)',26:'Replace moment-arm emphasis with predicted torque envelope (17)',27:'Move modeled pulley extensions to backup (34)',28:'Replace stale corpus text with current resource scope (21)',29:'Move curation/AI workflow to notes',30:'Replace implied parameter mapping with source-interpretation limit (22)',31:'Separate initial model (23) from later bilateral extension (24)',32:'Scope ownership of existing packages and added interfaces (25)',33:'Retain matched rhythm evidence (27)',34:'Supported-motion backup (35); add explicit bridge boundary (28)',35:'Move generalization warning to scope notes; remove score catalogue',36:'Remove broad physiological-equivalence claim; representative reflex scope (26)',37:'Move selected spiking comparison to backup (36)',38:'Merge coding caveats into backup notes (36)',39:'Replace synthesis ladder with evidence/boundary comparison (29)',40:'Integrate limitations beside the relevant evidence',41:'Focus next experiments on measured gaps (30)',42:'Merge acknowledgments with questions (31) and notes',43:'Questions ends the main sequence (31)',44:'Force-model backup (32)',45:'Keep exact coefficients in dissertation; do not repeat table',46:'Keep complete model-comparison statistics in dissertation/notes',47:'Move mechanics explanation to notes and dissertation',48:'Move implementation detail to notes and dissertation',49:'Identification design visible in main (12), details in notes',50:'Move stiffness invariance to dissertation; no main claim needed',51:'Move bracket coordinates/terms to provenance notes',52:'Keep sign conventions in torque notes; target backup (37)',53:'Representative library/reflex main (25–26); coding backup (36)',54:'Remove stale walker score table; use explicit verification boundaries'}
for i in range(1,55):lines.append(f'| {i} | {decisions[i]} |')
(root/'SLIDE_REVIEW.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
# Preserve reference identity, enumerate slide visibility and scientific asset hashes.
src=root.parent/'Bolen_Dissertation_Defense.pptx'
ref=root/'reference/Original_Defense.pptx'
checks={'original_deck_unchanged':hashlib.sha256(src.read_bytes()).hexdigest()==hashlib.sha256(ref.read_bytes()).hexdigest(),'original_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'slide_count':len(specs),'main_slides':sum(not x['backup'] for x in specs),'hidden_backups':sum(x['backup'] for x in specs)}
deck=root/'output/Bolen_Defense_Assertion_Evidence_20261023_v6.pptx'
with zipfile.ZipFile(deck) as z:
 ns={'p':'http://schemas.openxmlformats.org/presentationml/2006/main'}
 hidden=[]
 for i in range(1,38):
  x=etree.fromstring(z.read(f'ppt/slides/slide{i}.xml'))
  if x.get('show')=='0':hidden.append(i)
 checks['actual_hidden_slides']=hidden
 checks['notes_count']=len([n for n in z.namelist() if n.startswith('ppt/notesSlides/notesSlide') and n.endswith('.xml')])
checks['final_deck_sha256']=hashlib.sha256(deck.read_bytes()).hexdigest()
(root/'build/conservation_checks.json').write_text(json.dumps(checks,indent=2))
print(checks)



