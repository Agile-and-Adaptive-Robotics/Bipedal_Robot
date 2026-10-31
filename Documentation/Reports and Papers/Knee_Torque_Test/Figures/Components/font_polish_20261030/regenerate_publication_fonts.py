"""Reformat frozen vector plots without re-evaluating any scientific series.

Run with Python containing PyMuPDF. Sources are the pre-formatting PDFs in the
paper's Deprecated/font_polish_20261030 folder. Original trace operators and
clipping are retained; only text and the original legend graphics are replaced.
"""
from pathlib import Path
import hashlib
import json
import math
import re
import fitz

FIGURES = Path(__file__).resolve().parents[2]
SOURCE = FIGURES / 'Deprecated/font_polish_20261030'
ARIAL = Path('C:/Windows/Fonts/arial.ttf')
ARIAL_BOLD = Path('C:/Windows/Fonts/arialbd.ttf')
FONT = fitz.Font(fontfile=str(ARIAL_BOLD))
SIZE = 17.4

def strip_legend(stream, rectangle):
    token = re.escape(rectangle) + r' re '
    pattern = token + r'f.*?' + token + r'S'
    out, count = re.subn(pattern, '[]0 d .129883 G 0 j 20.025 w 1 g ', stream, count=1, flags=re.S)
    assert count == 1, rectangle
    return out

def draw_text(page, text, center, size=SIZE, bold=True, angle=0, align='center'):
    font = FONT if bold else fitz.Font(fontfile=str(ARIAL))
    width = font.text_length(text, fontsize=size)
    x,y = center
    x -= width/2 if align == 'center' else width if align == 'right' else 0
    baseline = y + size*(font.ascender+font.descender)/2
    morph = (fitz.Point(*center), fitz.Matrix(angle)) if angle else None
    page.insert_text((x,baseline),text,fontname='ArialBold' if bold else 'Arial',fontfile=str(ARIAL_BOLD if bold else ARIAL),fontsize=size,color=(.13,.13,.13),morph=morph)

def legend(page, x, y, label, kind, color, dash=None):
    if kind=='point':page.draw_circle((x+16,y),3.9,color=None,fill=color)
    else:page.draw_line((x,y),(x+32,y),color=color,width=2.2,dashes=dash)
    draw_text(page,label,(x+41,y),align='left')

def angle_label(page, x, y):
    # Arial has theta and degree glyphs but no subscript-k glyph.
    draw_text(page,'θ',(x-19,y))
    draw_text(page,'k',(x-5,y+5))
    draw_text(page,'(°)',(x+19,y))

def regenerate(name):
    source=SOURCE/name
    raw=source.read_bytes()
    original=fitz.open(stream=raw,filetype='pdf')
    op=original[0]
    spans=[(line,span) for block in op.get_text('dict')['blocks'] for line in block.get('lines',[]) for span in line['spans']]
    cleaned=fitz.open(stream=raw,filetype='pdf')
    cp=cleaned[0]
    stream=cp.read_contents().decode('latin1')
    if name.startswith('ExtPin'):
        stream=strip_legend(stream,'6465 4762.5 870 405')
        width,height,dx,dy=800,655,40,27
    else:
        stream=strip_legend(stream,'967.5 495 1620 397.5')
        stream=strip_legend(stream,'5715 2730 1695 405')
        width,height,dx,dy=787,505,15,35
    # Delete text objects only. All original trace/path coordinates, colors,
    # linewidths, dash patterns, and clipping operators remain in the stream.
    stream=re.sub(r'\bBT\b.*?\bET\b','',stream,flags=re.S)
    cp.clean_contents()
    xref=cp.get_contents()[0]
    cleaned.update_stream(xref,stream.encode('latin1'))
    output=fitz.open();page=output.new_page(width=width,height=height)
    page.show_pdf_page(fitz.Rect(dx,dy,dx+op.rect.width,dy+op.rect.height),cleaned,0)
    # The original degree sign was vectorized separately from the text.
    page.draw_rect(fitz.Rect(370+dx,op.rect.height-22+dy,450+dx,op.rect.height+2+dy),color=None,fill=(1,1,1))
    for line,span in spans:
        text=span['text'].strip();box=fitz.Rect(span['bbox'])
        if not text or span['size']<10 or text in ['θ','k',',','(A)','(B)','(C)','(D','(D)','(E)','(F)'] or text=='Torque, N m':continue
        center=((box.x0+box.x1)/2+dx,(box.y0+box.y1)/2+dy)
        if text in ['Training','Validation'] or 'cm' in text or 'route variants' in text:
            draw_text(page,text,center)
        elif re.fullmatch(r'-?\d+',text):
            angle=30 if line['dir'][1]<-.1 else 0
            # Right alignment keeps enlarged y ticks clear of their axes.
            if not angle and box.width<45 and (name.startswith('ExtPin') or box.x1<89 or 400<box.x1<469):
                draw_text(page,text,(box.x1+dx,center[1]),align='right')
            else:draw_text(page,text,center,angle=angle)
    if name.startswith('ExtPin'):
        for label,x,y in [('A',.75,20),('B',259.5,20),('C',506.25,20),('D',.75,297),('E',259.5,297),('F',506.25,297)]:
            draw_text(page,f'({label})',(x+dx,y+dy),align='left')
        draw_text(page,'Torque (N m)',(18,282),angle=90)
        angle_label(page,400,591)
        entries=[(40,'Hybrid','point',(1,.88,.25),None),(210,'Measured','point',(.25,.25,1),None),(400,'Unoptimized','line',(.98,.529,.46),'[6 6]'),(610,'Optimized','line',(.801,.2,.711),None)]
        for x,label,kind,color,dash in entries:legend(page,x,624,label,kind,color,dash)
    else:
        draw_text(page,'(A)',(16,20),align='left');draw_text(page,'(B)',(416,20),align='left')
        draw_text(page,'Torque (N m)',(31,191),angle=90)
        angle_label(page,405,379)
        groups=[(65,[('Direct extensor route',(.98,.529,.46),None),('Cam extensor (design of record)',(.801,.2,.711),'[9 9]'),('Required envelope (+5%)',(.4,.4,.4),'[2.25 2.25]'),('Human target (vasti)',(0,0,0),'[6 1.5 3 1.5]')]),(430,[('Direct flexor route',(.98,.529,.46),None),('Mono-pulley flexor (insertion side)',(1,.84,0),'[9 9]'),('Required envelope (+5%)',(.4,.4,.4),'[2.25 2.25]'),('Human target (bifemsh)',(0,0,0),'[6 1.5 3 1.5]')])]
        for x,group in groups:
            for i,(label,color,dash) in enumerate(group):legend(page,x,411+23*i,label,'line',color,dash)
    target=FIGURES/name
    output.save(target,garbage=4,deflate=True)
    check=fitz.open(target)
    text_sizes=[s['size'] for b in check[0].get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans'] if s['text'].strip()]
    assert min(text_sizes)>17.3
    return {'file':str(target),'source_sha256':hashlib.sha256(raw).hexdigest(),'native_width_pt':width,'native_height_pt':height,'minimum_native_font_pt':min(text_sizes),'minimum_at_503pt_width':min(text_sizes)*503/width,'source_series':'Frozen original vector paths, no evaluator or simulation run','legend_layout':'Outside plot axes'}

if __name__=='__main__':
    result=[regenerate(name) for name in ['ExtPin_complete_publication.pdf','RouteVariants_20261020_publication.pdf']]
    (Path(__file__).parent/'FONT_VERIFICATION.json').write_text(json.dumps(result,indent=2),encoding='utf8')
    print(json.dumps(result,indent=2))
