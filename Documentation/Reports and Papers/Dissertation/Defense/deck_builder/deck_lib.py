"""Helpers for the Bolen dissertation defense deck (C3NS assertion-evidence style)."""
import os
from PIL import Image
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE, MSO_CONNECTOR
from pptx.chart.data import CategoryChartData
from pptx.enum.chart import XL_CHART_TYPE, XL_LABEL_POSITION, XL_TICK_MARK, XL_LEGEND_POSITION

# ---- palette (C3NS talk) ----
BROWN  = RGBColor(0xAF, 0x7B, 0x51)   # assertion titles
TEXT   = RGBColor(0x22, 0x22, 0x22)
MUTED  = RGBColor(0x59, 0x59, 0x59)
FAINT  = RGBColor(0x7F, 0x7F, 0x7F)
TEAL   = RGBColor(0x00, 0x79, 0x6B)   # labels / annotations
GREEN  = RGBColor(0x38, 0x76, 0x1D)   # ours / good
GREEN2 = RGBColor(0x6A, 0xA8, 0x4F)
RED    = RGBColor(0xC0, 0x00, 0x00)   # problem / baseline
BLUE   = RGBColor(0x2E, 0x5A, 0x88)
CARD   = RGBColor(0xF5, 0xF2, 0xEE)   # light warm card
CARD2  = RGBColor(0xEF, 0xF3, 0xF1)   # light teal-gray card
WHITE  = RGBColor(0xFF, 0xFF, 0xFF)
LINEGR = RGBColor(0xD9, 0xD4, 0xCE)

FONT = "Calibri"
W, H = 13.333, 7.5
MX = 0.55

ASSETS = r"D:\GitHub\Bipedal_Robot\.zcode_tmp\deck_assets"
FIG = os.path.join(ASSETS, "figs")
EQ  = os.path.join(ASSETS, "eqs")
MEDIA = os.path.join(ASSETS, "talk_media")

_AR = {}
def ar(path):
    if path not in _AR:
        with Image.open(path) as im:
            _AR[path] = im.size[0] / im.size[1]
    return _AR[path]

def new_deck():
    prs = Presentation()
    prs.slide_width = Inches(W)
    prs.slide_height = Inches(H)
    return prs

def _bg(slide):
    slide.background.fill.solid()
    slide.background.fill.fore_color.rgb = WHITE

bg = _bg  # public alias (star-import skips underscore names)

def tx(slide, x, y, w, h, text, size=16, color=TEXT, bold=False, italic=False,
       align=PP_ALIGN.LEFT, font=FONT, anchor=MSO_ANCHOR.TOP, spacing=1.0, wrap=True):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = wrap
    tf.vertical_anchor = anchor
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    lines = text.split("\n") if isinstance(text, str) else text
    for i, ln in enumerate(lines):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        p.line_spacing = spacing
        r = p.add_run(); r.text = ln
        r.font.size = Pt(size); r.font.bold = bold; r.font.italic = italic
        r.font.color.rgb = color; r.font.name = font
    return tb

def rich(slide, x, y, w, h, paras, align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.TOP, spacing=1.0):
    """paras: list of (text, size, color, bold, italic, space_after) tuples or dicts."""
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame; tf.word_wrap = True; tf.vertical_anchor = anchor
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    for i, spec in enumerate(paras):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = spec.get("align", align); p.line_spacing = spec.get("spacing", spacing)
        if spec.get("space_after") is not None:
            p.space_after = Pt(spec["space_after"])
        if spec.get("space_before") is not None:
            p.space_before = Pt(spec["space_before"])
        r = p.add_run(); r.text = spec["t"]
        r.font.size = Pt(spec.get("size", 16)); r.font.bold = spec.get("bold", False)
        r.font.italic = spec.get("italic", False)
        r.font.color.rgb = spec.get("color", TEXT); r.font.name = FONT
    return tb

def card(slide, x, y, w, h, fill=CARD, line=None, radius=0.08, shadow=False):
    sp = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(x), Inches(y), Inches(w), Inches(h))
    try:
        sp.adjustments[0] = radius
    except Exception:
        pass
    sp.fill.solid(); sp.fill.fore_color.rgb = fill
    if line is None:
        sp.line.fill.background()
    else:
        sp.line.color.rgb = line; sp.line.width = Pt(1)
    sp.shadow.inherit = False
    sp.text_frame.paragraphs[0].text = ""
    return sp

def pic(slide, path, x=0, y=0, w=None, h=None, box=None, align="center", valign="middle"):
    """Fit image by aspect ratio. Give w or h, or box=(bx,by,bw,bh)."""
    if box:
        bx, by, bw, bh = box
        a = ar(path)
        pw, ph = bw, bw / a
        if ph > bh:
            ph, pw = bh, bh * a
        px = bx + (bw - pw) / 2 if align == "center" else (bx if align == "left" else bx + bw - pw)
        py = by + (bh - ph) / 2 if valign == "middle" else (by if valign == "top" else by + bh - ph)
        return slide.shapes.add_picture(path, Inches(px), Inches(py), Inches(pw), Inches(ph))
    a = ar(path)
    if w and not h:
        h = w / a
    if h and not w:
        w = h * a
    return slide.shapes.add_picture(path, Inches(x), Inches(y), Inches(w), Inches(h))

def eq(slide, name, x=0, y=0, h=0.5, box=None):
    p = os.path.join(EQ, f"{name}.png")
    return pic(slide, p, x, y, box=box, h=h if not box else None)

def content_slide(prs, kicker, assertion, cite=None, notes=None, assert_size=26):
    s = prs.slides.add_slide(prs.slide_layouts[6])
    _bg(s)
    if kicker:
        tx(s, MX, 0.30, 12.2, 0.3, kicker.upper(), size=12, color=TEAL, bold=True)
    tx(s, MX, 0.60, 12.25, 1.25, assertion, size=assert_size, color=BROWN, bold=True, spacing=0.98)
    if cite:
        tx(s, MX, 7.08, 11.4, 0.32, cite, size=10, color=FAINT, italic=True)
    if notes:
        s.notes_slide.notes_text_frame.text = notes
    return s

def pageno(slide, n):
    tx(slide, 12.55, 7.08, 0.55, 0.32, str(n), size=11, color=FAINT, align=PP_ALIGN.RIGHT)

def stat(slide, x, y, w, value, label, color=BROWN, vsize=40, lsize=13, align=PP_ALIGN.CENTER):
    rich(slide, x, y, w, 1.6, [
        {"t": value, "size": vsize, "color": color, "bold": True, "align": align, "space_after": 2},
        {"t": label, "size": lsize, "color": MUTED, "align": align, "spacing": 0.95},
    ])

def node(slide, x, y, w, h, text, fill=CARD2, line=None, size=13, color=TEXT, bold=False, radius=0.12):
    c = card(slide, x, y, w, h, fill=fill, line=line, radius=radius)
    tf = c.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = Inches(0.06)
    tf.margin_top = tf.margin_bottom = Inches(0.03)
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    lines = text.split("\n")
    for i, ln in enumerate(lines):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = PP_ALIGN.CENTER
        p.line_spacing = 0.95
        r = p.add_run(); r.text = ln
        r.font.size = Pt(size); r.font.bold = bold; r.font.color.rgb = color; r.font.name = FONT
    return c

def arrow(slide, x1, y1, x2, y2, color=MUTED, width=2.0):
    cn = slide.shapes.add_connector(MSO_CONNECTOR.STRAIGHT, Inches(x1), Inches(y1), Inches(x2), Inches(y2))
    cn.line.color.rgb = color
    cn.line.width = Pt(width)
    ln = cn.line._get_or_add_ln()
    from pptx.oxml.ns import qn
    tail = ln.find(qn("a:tailEnd"))
    if tail is None:
        tail = ln.makeelement(qn("a:tailEnd"), {})
        ln.append(tail)
    tail.set("type", "triangle"); tail.set("w", "med"); tail.set("len", "med")
    return cn

def bar_chart(slide, x, y, w, h, cats, vals, colors, title=None, num_fmt="0.0",
              label_size=13, cat_size=13, max_scale=None, gap=60):
    cd = CategoryChartData()
    cd.categories = cats
    cd.add_series("s", vals)
    gf = slide.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(x), Inches(y), Inches(w), Inches(h), cd)
    ch = gf.chart
    ch.has_legend = False
    ch.has_title = bool(title)
    if title:
        ch.chart_title.text_frame.text = title
        for r in ch.chart_title.text_frame.paragraphs[0].runs:
            r.font.size = Pt(13); r.font.bold = False; r.font.color.rgb = MUTED; r.font.name = FONT
    plot = ch.plots[0]
    plot.gap_width = gap
    plot.vary_by_categories = False
    ser = plot.series[0]
    for i, pt in enumerate(ser.points):
        pt.format.fill.solid()
        pt.format.fill.fore_color.rgb = colors[i % len(colors)]
    plot.has_data_labels = True
    dl = plot.data_labels
    dl.number_format = num_fmt; dl.number_format_is_linked = False
    dl.position = XL_LABEL_POSITION.OUTSIDE_END
    dl.font.size = Pt(label_size); dl.font.bold = True; dl.font.color.rgb = TEXT; dl.font.name = FONT
    ca = ch.category_axis
    ca.tick_labels.font.size = Pt(cat_size); ca.tick_labels.font.name = FONT
    ca.tick_labels.font.color.rgb = TEXT
    ca.major_tick_mark = XL_TICK_MARK.NONE
    ca.format.line.color.rgb = LINEGR
    va = ch.value_axis
    va.has_major_gridlines = False
    va.visible = False
    va.major_tick_mark = XL_TICK_MARK.NONE
    if max_scale:
        va.maximum_scale = max_scale
    return ch
