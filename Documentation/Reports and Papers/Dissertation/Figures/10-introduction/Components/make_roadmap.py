"""Export the dissertation roadmap as editable SVG and vector PDF."""
from pathlib import Path
from math import hypot
from html import escape
import argparse
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.colors import HexColor

WIDTH, HEIGHT = 468, 242.04
INK, LIGHT = '#243343', '#F1F5F8'
SIZE = 10.8
ARROW_LENGTH, ARROW_HALF_WIDTH = 8, 3.5

def build(pdf_path, svg_path):
    pdfmetrics.registerFont(TTFont('Arial',r'C:\Windows\Fonts\arial.ttf'))
    pdfmetrics.registerFont(TTFont('Arial-Bold',r'C:\Windows\Fonts\arialbd.ttf'))
    pdf=canvas.Canvas(str(pdf_path),pagesize=(WIDTH,HEIGHT),initialFontName='Arial',invariant=1)
    pdf.setTitle('Dissertation roadmap')
    pdf.setAuthor('Ben Bolen')
    svg=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}pt" height="{HEIGHT}pt" viewBox="0 0 {WIDTH} {HEIGHT}">',
         '<rect width="100%" height="100%" fill="white"/>']

    def line(points,dashed=False,arrow=False):
        pdf.setStrokeColor(HexColor(INK));pdf.setFillColor(HexColor(INK));pdf.setLineWidth(1)
        pdf.setDash([4,3] if dashed else [])
        path=pdf.beginPath();path.moveTo(points[0][0],HEIGHT-points[0][1])
        for x,y in points[1:]:path.lineTo(x,HEIGHT-y)
        pdf.drawPath(path)
        style=' stroke-dasharray="4 3"' if dashed else ''
        svg.append(f'<polyline points="'+ ' '.join(f'{x},{y}' for x,y in points)+f'" fill="none" stroke="{INK}" stroke-width="1"{style}/>')
        pdf.setDash([])
        if arrow:
            ax,ay=points[-2];bx,by=points[-1]
            length=hypot(bx-ax,by-ay);ux,uy=(bx-ax)/length,(by-ay)/length
            tri=[(bx,by),
                 (bx-ARROW_LENGTH*ux+ARROW_HALF_WIDTH*uy,by-ARROW_LENGTH*uy-ARROW_HALF_WIDTH*ux),
                 (bx-ARROW_LENGTH*ux-ARROW_HALF_WIDTH*uy,by-ARROW_LENGTH*uy+ARROW_HALF_WIDTH*ux)]
            p=pdf.beginPath();p.moveTo(tri[0][0],HEIGHT-tri[0][1])
            for x,y in tri[1:]:p.lineTo(x,HEIGHT-y)
            p.close();pdf.drawPath(p,stroke=0,fill=1)
            svg.append(f'<polygon points="'+ ' '.join(f'{x},{y}' for x,y in tri)+f'" fill="{INK}"/>')

    def text(x,top,value,bold=False,anchor='middle',size=SIZE):
        font='Arial-Bold' if bold else 'Arial'
        baseline=top+pdfmetrics.getAscent(font)*size/1000
        pdf.setFont(font,size);pdf.setFillColor(HexColor(INK))
        if anchor=='middle':pdf.drawCentredString(x,HEIGHT-baseline,value)
        elif anchor=='end':pdf.drawRightString(x,HEIGHT-baseline,value)
        else:pdf.drawString(x,HEIGHT-baseline,value)
        weight='bold' if bold else 'normal'
        svg.append(f'<text x="{x}" y="{baseline}" text-anchor="{anchor}" font-family="Arial, sans-serif" font-size="{size}" font-weight="{weight}" fill="{INK}">{escape(value)}</text>')

    def box(x,y,title,subtitle,dashed=False):
        w,h=112,51
        pdf.setStrokeColor(HexColor(INK));pdf.setFillColor(HexColor(LIGHT));pdf.setLineWidth(1)
        pdf.setDash([4,3] if dashed else [])
        pdf.rect(x,HEIGHT-y-h,w,h,stroke=1,fill=1);pdf.setDash([])
        style=' stroke-dasharray="4 3"' if dashed else ''
        svg.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="{LIGHT}" stroke="{INK}" stroke-width="1"{style}/>')
        a,b=title.split('\n'),subtitle.split('\n')
        block_height=(len(a)+len(b))*SIZE*1.10+5
        top=y+(h-block_height)/2+.5
        for i,value in enumerate(a):text(x+w/2,top+i*SIZE*1.10,value,bold=True)
        top+=len(a)*SIZE*1.10+5
        for i,value in enumerate(b):text(x+w/2,top+i*SIZE*1.10,value)

    line([(15,29),(36,29)]);text(43,22,'Chapter relationships',anchor='start')
    line([(295,29),(316,29)],dashed=True);text(323,22,'Future integration',anchor='start')
    nodes=[(28,53,'Background','(Ch. 2)',False),
           (178,53,'BPA force models','Characterization\n(Ch. 3)',False),
           (328,53,'Torque model','Design factors\nand routes (Ch. 4)',False),
           (28,153,'Sensory Afferent\nDatabase','Pathways (Ch. 5)',False),
           (178,153,'SNS simulation','Software tools\n(Ch. 6)',False),
           (328,153,'BPA robot','SNS control\n(future work)',True)]
    for node in nodes:box(*node)
    line([(140,78.5),(178,78.5)],arrow=True)
    line([(84,104),(84,153)],arrow=True)
    line([(290,78.5),(328,78.5)],arrow=True)
    line([(140,178.5),(178,178.5)],arrow=True)
    line([(290,178.5),(328,178.5)],dashed=True,arrow=True)
    line([(384,104),(384,126),(234,126),(234,153)],arrow=True)
    line([(384,204),(384,225),(12,225),(12,78.5),(28,78.5)],dashed=True,arrow=True)
    pdf.showPage();pdf.save()
    svg.append('</svg>')
    svg_path.write_text('\n'.join(svg),encoding='utf-8')

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--pdf',type=Path,default=Path(__file__).resolve().parent.parent/'roadmap.pdf')
    parser.add_argument('--svg',type=Path,default=Path(__file__).resolve().with_name('roadmap.svg'))
    args=parser.parse_args()
    build(args.pdf,args.svg)
