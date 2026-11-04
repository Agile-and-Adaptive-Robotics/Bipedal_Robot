from pathlib import Path
import hashlib,json,zipfile,posixpath
from lxml import etree as E
from PIL import Image
from verify_delivery import slide_parts,N
W=Path(__file__).resolve().parent;B=W/'build/plot_rotation_180'
deck=W/'output/Bolen_PhD_Defense_RoseGreen_20261102_v21.pptx'
old=(B/'before.png').read_bytes()
new=(W/'assets/mechanics/historical_bflh_common_scales.png').read_bytes()
sha=lambda b:hashlib.sha256(b).hexdigest()
assert Image.open(B/'before.png').size==Image.open(W/'assets/mechanics/historical_bflh_common_scales.png').size
def owners(z,targets):
    result={}
    for n,part in enumerate(slide_parts(z),1):
        relfile=posixpath.join(posixpath.dirname(part),'_rels',posixpath.basename(part)+'.rels')
        rels={r.get('Id'):posixpath.normpath(posixpath.join(posixpath.dirname(part),r.get('Target'))).lstrip('/') for r in E.fromstring(z.read(relfile)) if r.get('TargetMode')!='External'}
        pics=E.fromstring(z.read(part)).findall('.//a:blip',N)
        hits=[rels.get(p.get('{'+N['r']+'}embed')) for p in pics]
        hits=[x for x in hits if x in targets]
        if hits:result[str(n)]=hits
    return result
with zipfile.ZipFile(B/'artifact_images.pptx') as a:
    sdkparts=[n for n in a.namelist() if n.startswith('ppt/media/') and sha(a.read(n))==sha(new)]
    assert sorted(owners(a,sdkparts))==['33','8'],(sdkparts,owners(a,sdkparts))
    replacement=a.read(sdkparts[0])
with zipfile.ZipFile(deck) as z:
    targets=[n for n in z.namelist() if n.startswith('ppt/media/') and sha(z.read(n))==sha(old)]
    mapping=owners(z,targets);assert sorted(mapping)==['33','8']
    original={n:sha(z.read(n)) for n in z.namelist()}
    with zipfile.ZipFile(B/'candidate.pptx','w') as out:
        for entry in z.infolist():out.writestr(entry,replacement if entry.filename in targets else z.read(entry.filename))
with zipfile.ZipFile(B/'candidate.pptx') as z:
    changed=[n for n in z.namelist() if sha(z.read(n))!=original[n]]
assert set(changed)==set(targets)
report={'original_deck_sha256':sha(deck.read_bytes()),'candidate_sha256':sha((B/'candidate.pptx').read_bytes()),'changed_parts':changed,'affected_slides':[8,33],'owners':mapping,'all_other_package_parts_identical':True,'old_image_sha256':sha(old),'new_image_sha256':sha(new),'original_part_hashes':original}
(B/'package_check.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k!='original_part_hashes'}))
