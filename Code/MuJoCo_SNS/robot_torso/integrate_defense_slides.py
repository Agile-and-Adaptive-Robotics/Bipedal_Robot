"""Insert the verified supplement into a new copy of Ben's defense deck."""
from pathlib import Path
import json,hashlib,shutil,zipfile,xml.etree.ElementTree as ET
import win32com.client
ROOT=Path(__file__).resolve().parents[3]
DEFENSE=ROOT/'Documentation/Reports and Papers/Dissertation/Defense'
ASSETS=DEFENSE/'slideshow/AH001_20261028'
source=DEFENSE/'AE_Revision_20261023/output/Bolen_Defense_Assertion_Evidence_20261023_v6.pptx'
supplement=ASSETS/'presentation/AH001_supported_control_20261028.pptx'
candidate=ASSETS/'.build/combined_candidate.pptx'
final=source.with_name('Bolen_Defense_Assertion_Evidence_20261028_AH001_v8.pptx')
assert not final.exists(),'Do not overwrite an existing final deck.'
original_hash=hashlib.sha256(source.read_bytes()).hexdigest()
shutil.copy2(source,candidate)
app=win32com.client.Dispatch('PowerPoint.Application')
deck=app.Presentations.Open(str(candidate),ReadOnly=False,Untitled=False,WithWindow=False)
def signature(slide):
    texts=[]
    for shape in slide.Shapes:
        if shape.HasTextFrame and shape.TextFrame.HasText:
            text=shape.TextFrame.TextRange.Text
            if not (text.strip().isdigit() and shape.Top>500 and shape.Left>850):texts.append(text)
    notes=[]
    for shape in slide.NotesPage.Shapes:
        if shape.HasTextFrame and shape.TextFrame.HasText:notes.append(shape.TextFrame.TextRange.Text)
    return dict(id=slide.SlideID,texts=texts,notes=notes)
try:
    before=[signature(s) for s in deck.Slides]
    assert len(before)==37
    assert deck.Slides.InsertFromFile(str(supplement),28,1,1)==1
    assert deck.Slides.InsertFromFile(str(supplement),deck.Slides.Count,2,2)==1
    deck.Slides(29).SlideShowTransition.Hidden=0
    deck.Slides(39).SlideShowTransition.Hidden=-1
    # Native PowerPoint renumbers fields; explicitly update plain page labels.
    for number,slide in enumerate(deck.Slides,1):
        for shape in slide.Shapes:
            if shape.HasTextFrame and shape.TextFrame.HasText:
                if shape.TextFrame.TextRange.Text.strip().isdigit() and shape.Top>500 and shape.Left>850:
                    shape.TextFrame.TextRange.Text=str(number)
                    if number in (29,39):
                        shape.Left=918;shape.Width=30;shape.TextFrame.WordWrap=0
    originals=[signature(s) for i,s in enumerate(deck.Slides,1) if i not in (29,39)]
    assert originals==before,'Original slide text, notes, or identity changed.'
    deck.SaveAs(str(final),24)
    for number in (29,39):deck.Slides(number).Export(str(ASSETS/'presentation'/('defense_slide_%d.png'%number)),'PNG',1600,900)
    count=deck.Slides.Count
finally:deck.Close()
assert hashlib.sha256(source.read_bytes()).hexdigest()==original_hash
with zipfile.ZipFile(final) as z:
    assert z.testzip() is None
    gif_hashes={p:hashlib.sha256(z.read(p)).hexdigest() for p in z.namelist() if p.startswith('ppt/media/') and p.lower().endswith('.gif')}
    asset_hash=hashlib.sha256((ASSETS/'ah001_abductor_control.gif').read_bytes()).hexdigest()
    assert asset_hash in gif_hashes.values(),'Animated GIF not preserved byte for byte.'
result=dict(source=str(source),source_unchanged=True,source_sha256=original_hash,output=str(final),
    final_sha256=hashlib.sha256(final.read_bytes()).hexdigest(),slides=count,inserted_animation_slide=29,
    inserted_hidden_still_comparison_slide=39,all_original_slide_text_notes_and_identity_preserved=True,
    animated_gif_embedded_byte_exact=True,package_crc_pass=True)
(ASSETS/'presentation/integration_verification.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
