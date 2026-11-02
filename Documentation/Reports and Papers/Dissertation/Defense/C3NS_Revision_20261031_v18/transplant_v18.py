from __future__ import annotations

from io import BytesIO
from pathlib import Path, PurePosixPath
import copy
import json
import posixpath
import zipfile

from lxml import etree as E
from openpyxl import load_workbook

W = Path(__file__).resolve().parent
BASE = W.parent / 'C3NS_Revision_20261031_v17/output/Bolen_PhD_Defense_RoseGreen_20261031_v17.pptx'
SOURCE = W / 'build/Bolen_PhD_Defense_RoseGreen_20261031_v18_artifact_draft.pptx'
OUTPUT = W / 'build/candidate.pptx'

NS = {
    'p': 'http://schemas.openxmlformats.org/presentationml/2006/main',
    'a': 'http://schemas.openxmlformats.org/drawingml/2006/main',
    'r': 'http://schemas.openxmlformats.org/officeDocument/2006/relationships',
    'pr': 'http://schemas.openxmlformats.org/package/2006/relationships',
    'c': 'http://schemas.openxmlformats.org/drawingml/2006/chart',
}
REL_NS = NS['pr']
R_NS = NS['r']
EMU_PER_POINT = 12700


def rels_part(part: str) -> str:
    p = PurePosixPath(part)
    return str(p.parent / '_rels' / (p.name + '.rels'))


def resolve_target(owner_part: str, target: str) -> str:
    if target.startswith('/'):
        return target.lstrip('/')
    return posixpath.normpath(posixpath.join(posixpath.dirname(owner_part), target))


def relative_target(owner_part: str, target_part: str) -> str:
    return posixpath.relpath(target_part, posixpath.dirname(owner_part))


def presentation_slide_parts(package: dict[str, bytes]) -> list[str]:
    pres = E.fromstring(package['ppt/presentation.xml'])
    rels = E.fromstring(package['ppt/_rels/presentation.xml.rels'])
    by_id = {node.get('Id'): node.get('Target') for node in rels}
    parts = []
    for node in pres.findall('p:sldIdLst/p:sldId', NS):
        rid = node.get(f'{{{R_NS}}}id')
        target = by_id[rid]
        parts.append(target.lstrip('/') if target.startswith('/') else posixpath.normpath(posixpath.join('ppt', target)))
    return parts


def read_package(path: Path) -> dict[str, bytes]:
    with zipfile.ZipFile(path) as z:
        return {name: z.read(name) for name in z.namelist()}


def relation_type(node) -> str:
    return (node.get('Type') or '').rsplit('/', 1)[-1]


def replace_notes_text(package: dict[str, bytes], notes_part: str, text: str) -> None:
    root = E.fromstring(package[notes_part])
    body = None
    for shape in root.findall('.//p:sp', NS):
        ph = shape.find('p:nvSpPr/p:nvPr/p:ph', NS)
        if ph is not None and ph.get('type') in (None, 'body'):
            if shape.findall('.//a:t', NS):
                body = shape
                if ph.get('type') == 'body':
                    break
    if body is None:
        return
    nodes = body.findall('.//a:t', NS)
    nodes[0].text = text
    for node in nodes[1:]:
        node.text = ''
    package[notes_part] = E.tostring(root, xml_declaration=True, encoding='UTF-8', standalone=True)


def transplant_slide(number: int, src: dict[str, bytes], dst: dict[str, bytes], src_parts: list[str], dst_parts: list[str]) -> dict:
    src_part = src_parts[number - 1]
    dst_part = dst_parts[number - 1]
    src_rels_part = rels_part(src_part)
    dst_rels_part = rels_part(dst_part)
    src_rels = E.fromstring(src[src_rels_part])
    old_dst_rels = E.fromstring(dst[dst_rels_part])
    old_by_type = {}
    for rel in old_dst_rels:
        old_by_type.setdefault(relation_type(rel), []).append(rel)

    media_copied = []
    notes_copied = False
    for rel in src_rels:
        kind = relation_type(rel)
        mode = rel.get('TargetMode')
        if mode == 'External':
            continue
        if kind == 'slideLayout' and old_by_type.get(kind):
            rel.set('Target', old_by_type[kind][0].get('Target'))
        elif kind == 'notesSlide' and old_by_type.get(kind):
            src_notes = resolve_target(src_part, rel.get('Target'))
            dst_target = old_by_type[kind][0].get('Target')
            dst_notes = resolve_target(dst_part, dst_target)
            rel.set('Target', dst_target)
            if src_notes in src and dst_notes in dst:
                # The notes XML contains no media relationships in this deck. Keep
                # the destination notes relationships so its slide backlink remains valid.
                dst[dst_notes] = src[src_notes]
                notes_copied = True
        elif kind in {'image', 'audio', 'video', 'media'}:
            src_media = resolve_target(src_part, rel.get('Target'))
            suffix = PurePosixPath(src_media).suffix
            stem = PurePosixPath(src_media).stem
            unique_name = f'v18_s{number:03d}_{rel.get("Id")}_{stem}{suffix}'
            dst_media = f'ppt/media/{unique_name}'
            dst[dst_media] = src[src_media]
            rel.set('Target', relative_target(dst_part, dst_media))
            media_copied.append(dst_media)
        else:
            target = resolve_target(src_part, rel.get('Target'))
            if target not in dst:
                raise RuntimeError(f'Unhandled internal relationship on slide {number}: {kind} -> {target}')
            rel.set('Target', relative_target(dst_part, target))

    dst[dst_part] = src[src_part]
    dst[dst_rels_part] = E.tostring(src_rels, xml_declaration=True, encoding='UTF-8', standalone=True)
    return {'slide': number, 'media': media_copied, 'notes': notes_copied, 'part': dst_part}


def set_slide_title(package: dict[str, bytes], slide_part: str, title: str) -> None:
    root = E.fromstring(package[slide_part])
    selected = None
    for shape in root.findall('p:cSld/p:spTree/p:sp', NS):
        xfrm = shape.find('p:spPr/a:xfrm', NS)
        if xfrm is None:
            continue
        off = xfrm.find('a:off', NS)
        ext = xfrm.find('a:ext', NS)
        if off is None or ext is None:
            continue
        if int(off.get('y', '99999999')) < 40 * EMU_PER_POINT and int(ext.get('cx', '0')) > 500 * EMU_PER_POINT:
            selected = shape
            break
    if selected is None:
        raise RuntimeError('Could not locate slide 44 title')
    nodes = selected.findall('.//a:t', NS)
    nodes[0].text = title
    for node in nodes[1:]:
        node.text = ''
    for props in selected.findall('.//a:rPr', NS) + selected.findall('.//a:defRPr', NS):
        props.set('sz', '3000')
    package[slide_part] = E.tostring(root, xml_declaration=True, encoding='UTF-8', standalone=True)


def update_slide44_chart(package: dict[str, bytes], slide_part: str) -> dict:
    slide_rels = E.fromstring(package[rels_part(slide_part)])
    chart_rel = next(rel for rel in slide_rels if relation_type(rel) == 'chart')
    chart_part = resolve_target(slide_part, chart_rel.get('Target'))
    chart_root = E.fromstring(package[chart_part])

    new_values = {
        'Redesign': {'Torque coverage': 65.81, 'Active-force coverage': 99.31},
        'Baseline': {'Torque coverage': 37.35, 'Active-force coverage': 93.00},
    }
    cache_updates = {}
    for ser in chart_root.findall('.//c:ser', NS):
        name_node = ser.find('.//c:tx//c:v', NS)
        if name_node is None or name_node.text not in new_values:
            continue
        name = name_node.text
        categories = {}
        for point in ser.findall('./c:cat//c:pt', NS):
            value = point.find('c:v', NS)
            if value is not None:
                categories[point.get('idx')] = value.text
        for point in ser.findall('./c:val//c:pt', NS):
            label = categories.get(point.get('idx'))
            if label in new_values[name]:
                value = point.find('c:v', NS)
                value.text = f'{new_values[name][label]:.4f}'
                cache_updates[(name, label)] = float(value.text)
        format_code = ser.find('./c:val//c:formatCode', NS)
        if format_code is not None:
            format_code.text = '0.00'

    package[chart_part] = E.tostring(chart_root, xml_declaration=True, encoding='UTF-8', standalone=True)

    chart_rels_part = rels_part(chart_part)
    chart_rels = E.fromstring(package[chart_rels_part])
    workbook_rel = next(rel for rel in chart_rels if relation_type(rel) == 'package')
    workbook_part = resolve_target(chart_part, workbook_rel.get('Target'))
    workbook = load_workbook(BytesIO(package[workbook_part]))
    updated_cells = []
    for ws in workbook.worksheets:
        headers = {}
        categories = {}
        for row in ws.iter_rows():
            for cell in row:
                if cell.value in new_values:
                    headers[cell.value] = cell.column
                if cell.value in ('Torque coverage', 'Active-force coverage'):
                    categories[cell.value] = cell.row
        if headers and categories:
            for series, col in headers.items():
                for category, row in categories.items():
                    ws.cell(row=row, column=col, value=new_values[series][category])
                    ws.cell(row=row, column=col).number_format = '0.00"%"'
                    updated_cells.append(f'{ws.title}!{ws.cell(row=row, column=col).coordinate}')
    if len(updated_cells) != 4:
        raise RuntimeError(f'Expected four workbook cells on slide 44, updated {updated_cells}')
    out = BytesIO()
    workbook.save(out)
    package[workbook_part] = out.getvalue()
    return {
        'chart_part': chart_part,
        'workbook_part': workbook_part,
        'cells': updated_cells,
        'cache': [{'series': key[0], 'category': key[1], 'value': value} for key, value in cache_updates.items()],
    }


src = read_package(SOURCE)
dst = read_package(BASE)
src_parts = presentation_slide_parts(src)
dst_parts = presentation_slide_parts(dst)
if len(src_parts) != 123 or len(dst_parts) != 123:
    raise RuntimeError(f'Unexpected slide counts: source={len(src_parts)} destination={len(dst_parts)}')

changed = [3, 4, 6, 8, 11, 49, 53, 54, 55, 56, 57, 58, 59, 60, 64, 65, 69]
report = [transplant_slide(number, src, dst, src_parts, dst_parts) for number in changed]

slide44 = dst_parts[43]
set_slide_title(dst, slide44, 'Route redesign raises mean torque coverage from 37.35% to 65.81%')
chart_report = update_slide44_chart(dst, slide44)

slide44_rels = E.fromstring(dst[rels_part(slide44)])
notes_rel = next((rel for rel in slide44_rels if relation_type(rel) == 'notesSlide'), None)
if notes_rel is not None:
    notes_part = resolve_target(slide44, notes_rel.get('Target'))
    replace_notes_text(
        dst,
        notes_part,
        'The submitted dissertation corrects the whole-leg comparison: mean torque coverage rises from 37.35% to 65.81%, while mean active-force coverage rises from 93.00% to 99.31%. Emphasize that route redesign expands the joint-angle conditions in which the robot model meets or exceeds the human benchmark.\n\nSource: Submitted dissertation, Chapter 4, Table 4.12, and Appendix C.',
    )

tmp = OUTPUT.with_suffix('.tmp')
with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as z:
    for name, data in dst.items():
        z.writestr(name, data)
tmp.replace(OUTPUT)

summary = {
    'baseline': str(BASE),
    'source': str(SOURCE),
    'output': str(OUTPUT),
    'transplanted': report,
    'slide44': chart_report,
    'preserved_from_v17': ['native charts and workbooks', 'hidden slide state', 'GIF animations', 'two AnimatLab MP4 recordings'],
}
(W / 'build/transplant_report.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
print(json.dumps({'slides_transplanted': len(report), 'output': str(OUTPUT), 'slide44_workbook_cells': chart_report['cells']}))
