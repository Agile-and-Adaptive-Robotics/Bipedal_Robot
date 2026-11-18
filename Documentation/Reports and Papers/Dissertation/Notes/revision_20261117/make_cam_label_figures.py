"""Place Appendix C4/C5 callouts on the original MATLAB vector artwork.

Run with the bundled document Python. Write candidates to the named directory.
The CAM_FIGURE_ORIGINALS archive is always the input, including on a rerun.
"""
import argparse
import hashlib
import io
import json
import math
from pathlib import Path
import re
import zipfile

import pdfplumber
from pypdf import PdfReader, PdfWriter
from pypdf.generic import DecodedStreamObject, NameObject
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas


HERE = Path(__file__).resolve().parent
ASSETS = 'Figures/94-AppendixC'
SIZES = {'AF': (490, 480), 'GH': (490, 360)}
PANELS = {'AF': range(6), 'GH': range(6, 8)}
TEXT = re.compile(rb'\bBT\b.*?\bET\b', re.DOTALL)
GLYPHS = re.compile(rb'<([0-9a-fA-F]{4})>\s*Tj')
GLYPH_DIGITS = {5: 1, 4: 2, 9: 3, 8: 4, 7: 5, 6: 6, 12: 7, 13: 8, 14: 9}
CALL_FONT = 'CamArialBold'
CALL_SIZE = 11.3
PLACED_WIDTH_PT = 468

# Centres of editorial labels in the original 490-point figures, with y measured
# down from the top. The actual leader endpoints come from the original vector.
LABEL_CENTRES = {
    'A': {3:(167,139), 4:(171,155), 5:(164,171), 6:(138,182), 7:(119,190), 8:(99,179), 9:(99,117)},
    'B': {3:(315,140), 4:(315,157), 6:(297,176), 7:(279,191), 8:(255,187), 9:(243,132)},
    'C': {3:(460,145), 6:(453,172), 7:(441,192), 8:(417,191), 9:(393,146)},
    'D': {3:(164,371), 6:(163,389), 7:(155,409), 8:(132,421), 9:(101,385)},
    'E': {6:(313,382), 7:(316,399), 8:(305,416), 9:(250,398)},
    'F': {6:(465,378), 7:(469,396), 8:(464,414), 9:(402,408)},
    'G': {7:(155,176), 8:(156,199), 9:(99,212)},
    'H': {7:(305,166), 8:(308,195), 9:(264,213)},
}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def drawing(page):
    return page['/Resources']['/XObject']['/fzFrm0']['/Resources']['/XObject']['/fullpage']


def filter_original_labels(content, expected):
    removed = []

    def remove(match):
        original = match.group()
        codes = [int(c, 16) for c in GLYPHS.findall(original)]
        if 11 not in codes or b'/F7 13 Tf' not in original:
            return original
        # Registered MATLAB label is [space, p, digit] and occupies its own BT.
        assert len(codes) == 3 and codes[:2] == [10, 11], codes
        removed.append(GLYPH_DIGITS[codes[2]])
        return b''

    preserved = TEXT.sub(remove, content)
    assert removed == expected, (removed, expected)
    return preserved


def existing_page_stream(data, key):
    if key != 'GH':
        return data
    return data.replace(b'/helv 10 Tf',b'/helv 11 Tf').replace(b'[4 3] d',b'[4 3] 0 d')


def blue_routes(data):
    with pdfplumber.open(io.BytesIO(data)) as doc:
        page = doc.pages[0]
        strokes = [curve['pts'] for curve in page.curves
                   if curve['stroking_color'] == (0., 0., 1.)
                   and curve['width'] > 15 and curve['height'] > 20]
    assert len(strokes) == 8
    return strokes


def validate_model_geometry(panel, path):
    """Fit the equal x/y plotted scale and translations to all active points."""
    points = panel['femur_xy_m']
    assert len(path) == len(panel['point_ids']) == len(points)
    n = len(points)
    mx, my = (sum(v[j] for v in points) / n for j in range(2))
    px, py = (sum(v[j] for v in path) / n for j in range(2))
    numerator = sum((x-mx)*(u-px) - (y-my)*(v-py)
                    for (x, y), (u, v) in zip(points, path))
    denominator = sum((x-mx)**2 + (y-my)**2 for x, y in points)
    scale = numerator / denominator
    tx, ty = px - scale*mx, py + scale*my
    err = max(math.hypot(tx + scale*x - u, ty - scale*y - v)
              for (x, y), (u, v) in zip(points, path))
    assert scale > 0 and err < 0.0001, (panel['panel'], err)
    return scale, err


def label_figure(key, routes, panels):
    width, height = SIZES[key]
    buffer = io.BytesIO()
    fig = canvas.Canvas(buffer, pagesize=(width, height), invariant=True,
                        pageCompression=1)
    fig.setFont(CALL_FONT, CALL_SIZE)
    asc, desc = pdfmetrics.getAscentDescent(CALL_FONT, CALL_SIZE)
    boxes = []
    records = []
    for i in PANELS[key]:
        panel = panels[i]
        label = panel['panel']
        scale, err = validate_model_geometry(panel, routes[i])
        records.append({'panel':label, 'angle_deg':panel['angle_deg'],
                        'scale_figure_pt_per_m':scale,
                        'maximum_model_to_vector_error_pt':err, 'points':[]})
        for pid, point in zip(panel['point_ids'], routes[i]):
            x, top = point
            if pid in (1, 2):
                cx, cy = x + 16.5, top - 2
            else:
                cx, cy = LABEL_CENTRES[label][pid]
            name = f'p{pid}'
            w = pdfmetrics.stringWidth(name,CALL_FONT,CALL_SIZE)
            hw, hh = w/2, (asc-desc)/2
            box = (cx-hw, cy-hh, cx+hw, cy+hh)
            assert 0 <= box[0] < box[2] <= width and 0 <= box[1] < box[3] <= height
            for code, (l,t,r,b) in boxes:
                assert box[2]+2 <= l or box[0]-2 >= r or box[3]+2 <= t or box[1]-2 >= b, (label,name,code)
            boxes.append((f'{label} {name}',box))

            # Leader is aimed at the mathematical centre. Stop on the original
            # marker's perimeter, so the printed fill identifies its material.
            dx, dy = cx-x, cy-top
            distance = math.hypot(dx,dy)
            unit_x, unit_y = dx/distance, dy/distance
            radius = 5.1 if pid in (1,9) else 4.5 if pid in (6,7,8) else 2.6
            label_pad = min((hw+1.3)/abs(unit_x) if unit_x else math.inf,
                            (hh+1.3)/abs(unit_y) if unit_y else math.inf)
            sx, sy = x + radius*unit_x, top + radius*unit_y
            ex, ey = cx - label_pad*unit_x, cy - label_pad*unit_y
            if distance > radius + label_pad + 0.5:
                fig.setStrokeColorRGB(.24,.24,.24)
                fig.setLineWidth(.68)
                fig.line(sx,height-sy,ex,height-ey)
            fig.setFillColorRGB(.09,.09,.09)
            fig.drawCentredString(cx,height-cy-(asc+desc)/2,name)
            records[-1]['points'].append({'id':name,'vector_xy_pt_top':point,
                                           'label_center_pt_top':[cx,cy]})
    fig.showPage()
    fig.save()
    return buffer.getvalue(), records


def rebuild(old, panels, key):
    original = PdfReader(io.BytesIO(old))
    writer = PdfWriter(clone_from=original)
    page = writer.pages[0]
    inner = drawing(page)
    old_vector = inner.get_data()
    expected = [p for panel in panels for p in panel['point_ids']]
    new_vector = filter_original_labels(old_vector, expected)
    inner.set_data(new_vector)

    # Older PDF did not supply the required dash phase for a C5 legend swatch.
    # Repair the PDF operator and enlarge existing small legend and note text.
    changes = {'key_or_note_font_changes':0, 'dash_phase_repairs':0}
    for ref in page['/Contents']:
        stream = ref.get_object()
        data = stream.get_data()
        source = data
        if key == 'GH':
            changes['key_or_note_font_changes'] += data.count(b'/helv 10 Tf')
            changes['dash_phase_repairs'] += data.count(b'[4 3] d')
            data = existing_page_stream(data,key)
        if source != data:
            stream.set_data(data)
    assert changes['key_or_note_font_changes'] == (13 if key == 'GH' else 0)
    assert changes['dash_phase_repairs'] == (1 if key == 'GH' else 0)

    routes = blue_routes(old)
    annotation, panels_checked = label_figure(key,routes,panels)
    anno = PdfReader(io.BytesIO(annotation)).pages[0]
    # Put the callouts in a separate form. Original plot placement and clipping
    # are preserved, as are all original page streams outside the edits above.
    form = DecodedStreamObject()
    form.set_data(anno.get_contents().get_data())
    form.update({NameObject('/Type'):NameObject('/XObject'),
                 NameObject('/Subtype'):NameObject('/Form'),
                 NameObject('/BBox'):anno.mediabox.clone(writer),
                 NameObject('/Resources'):anno['/Resources'].clone(writer)})
    page['/Resources']['/XObject'][NameObject('/CamPointLabels')] = writer._add_object(form.flate_encode())
    call = DecodedStreamObject()
    call.set_data(b'q\n/CamPointLabels Do\nQ\n')
    page['/Contents'].append(writer._add_object(call))
    out = io.BytesIO()
    writer.write(out)
    result = out.getvalue()

    reread = PdfReader(io.BytesIO(result)).pages[0]
    old_page = original.pages[0]
    assert list(reread.mediabox) == list(old_page.mediabox)
    old_streams = [s.get_object().get_data() for s in old_page['/Contents']]
    new_streams = [s.get_object().get_data() for s in reread['/Contents']]
    assert new_streams[:-1] == [existing_page_stream(s,key) for s in old_streams]
    assert new_streams[-1] == b'q\n/CamPointLabels Do\nQ\n'
    assert drawing(reread).get_data() == new_vector
    assert blue_routes(result) == routes
    wrapper = old_page['/Resources']['/XObject']['/fzFrm0']
    revised_wrapper = reread['/Resources']['/XObject']['/fzFrm0']
    assert revised_wrapper.get_data() == wrapper.get_data()
    for field in ('/BBox','/Matrix'):
        assert revised_wrapper[field] == wrapper[field]
        assert drawing(reread).get(field) == drawing(old_page).get(field)
    return result, dict(original_sha256=digest(old),revised_sha256=digest(result),
                        original_nested_plot_sha256=digest(old_vector),
                        nested_plot_sha256_after_deleting_old_labels=digest(new_vector),
                        point_label_size_pt_at_468pt=CALL_SIZE*PLACED_WIDTH_PT/490,
                        key_label_size_pt_at_468pt=11*PLACED_WIDTH_PT/490,
                        media_box=list(SIZES[key]), plot_clip=list(wrapper['/BBox']),
                        plot_placement_matrix=list(wrapper['/Matrix']),
                        retained_page_stream_count=len(old_streams),
                        original_plot_drawing_sha256_after_filter=digest(new_vector),
                        revised_plot_drawing_sha256=digest(drawing(reread).get_data()),
                        changes=changes,panels=panels_checked)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('output_dir',type=Path)
    ap.add_argument('--font',type=Path,default=Path('C:/Windows/Fonts/arialbd.ttf'))
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=True)
    pdfmetrics.registerFont(TTFont(CALL_FONT,str(args.font)))
    panels = json.loads((HERE/'CAM_LABEL_GEOMETRY.json').read_text())['panels']
    audit = {}
    with zipfile.ZipFile(HERE/'CAM_FIGURE_ORIGINALS.zip') as archive:
        manifest = {i['path']:i for i in json.loads(archive.read('manifest.json'))}
        for key in SIZES:
            name = f'extensor_cam_{key}_bones.pdf'
            member = f'{ASSETS}/{name}'
            source = archive.read(member)
            assert digest(source) == manifest[member]['sha256']
            data, audit[key] = rebuild(source,panels,key)
            (args.output_dir/name).write_bytes(data)
    (args.output_dir/'CAM_LABEL_VERIFICATION.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')


if __name__ == '__main__':
    main()
