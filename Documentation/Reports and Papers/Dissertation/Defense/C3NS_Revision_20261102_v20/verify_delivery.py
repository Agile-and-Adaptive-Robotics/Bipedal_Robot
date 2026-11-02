"""Read-only package checks; static verification does not certify playback."""
from pathlib import Path
import argparse, hashlib, io, json, posixpath, re, zipfile
from lxml import etree as E
from PIL import Image
from pypdf import PdfReader

W = Path(__file__).resolve().parent
TOTAL, MAIN = 127, 74
CHART_OWNERS = [15, 25, 30, 52, 101, 102]
SUBMITTED_SHA = '447c2cc5a9dd219fa0857e298727fb68d4797af8f406251865eb22f4b126361f'
N = {'p':'http://schemas.openxmlformats.org/presentationml/2006/main',
     'a':'http://schemas.openxmlformats.org/drawingml/2006/main',
     'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships',
     'c':'http://schemas.openxmlformats.org/drawingml/2006/chart'}

def sha(data):
    return hashlib.sha256(data).hexdigest()

def resolve(owner, target):
    return (posixpath.normpath(target.lstrip('/')) if target.startswith('/')
            else posixpath.normpath(posixpath.join(posixpath.dirname(owner), target)))

def relationships(z, part):
    rp = posixpath.join(posixpath.dirname(part), '_rels', posixpath.basename(part)+'.rels')
    if rp not in z.namelist():
        return {}
    return {r.get('Id'):{'kind':r.get('Type','').rsplit('/',1)[-1],
                        'part':resolve(part,r.get('Target')),
                        'external':r.get('TargetMode')=='External'}
            for r in E.fromstring(z.read(rp))}

def slide_parts(z):
    p = E.fromstring(z.read('ppt/presentation.xml'))
    rm = relationships(z, 'ppt/presentation.xml')
    return [rm[s.get('{'+N['r']+'}id')]['part'] for s in p.findall('p:sldIdLst/p:sldId',N)]

def text_of(z, part):
    return ' '.join(E.fromstring(z.read(part)).xpath('//a:t/text()',namespaces=N))

def linked_video_hashes(z, part):
    root, rm = E.fromstring(z.read(part)), relationships(z,part)
    used = {value for node in root.iter() for name,value in node.attrib.items()
            if name in ('{'+N['r']+'}embed','{'+N['r']+'}link')}
    return sorted({sha(z.read(rm[rid]['part'])) for rid in used if rid in rm
                   and not rm[rid]['external'] and rm[rid]['part'].lower().endswith('.mp4')
                   and rm[rid]['part'] in z.namelist()})

def branding_images(z, part):
    """Hash existing top-right logo/footer pictures and their exact placement."""
    root, rm, found = E.fromstring(z.read(part)), relationships(z,part), set()
    for pic in root.findall('./p:cSld/p:spTree/p:pic',N):
        xf, blip = pic.find('p:spPr/a:xfrm',N), pic.find('p:blipFill/a:blip',N)
        if xf is None or blip is None:
            continue
        off, ext = xf.find('a:off',N), xf.find('a:ext',N)
        if off is None or ext is None:
            continue
        box = tuple(int(v) for v in (off.get('x'),off.get('y'),ext.get('cx'),ext.get('cy')))
        x,y,width,height = box
        if not ((x >= .90*9144000 and y < .12*5143500) or y >= .95*5143500):
            continue
        rid = blip.get('{'+N['r']+'}embed')
        if rid in rm and not rm[rid]['external']:
            found.add((sha(z.read(rm[rid]['part'])),box))
    return found

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--target',type=Path,default=W/'output/Bolen_PhD_Defense_RoseGreen_20261102_v20.pptx')
    parser.add_argument('--skip-pdf',action='store_true',help='Package-only preflight; no claim of a current PDF render.')
    args = parser.parse_args()
    target = args.target.resolve()
    baseline = W.parent/'C3NS_Revision_20261102_v19/output/Bolen_PhD_Defense_RoseGreen_20261102_v19.pptx'
    failures = []
    def require(condition,label):
        if not condition:
            failures.append(label)
    with zipfile.ZipFile(target) as z, zipfile.ZipFile(baseline) as b:
        parts, oldparts = slide_parts(z), slide_parts(b)
        hidden, rejected, sadb_main, titles, missing = [],[],[],[],[]
        chart_owners, workbook_owners, video_owners = {},{},{}
        for i,part in enumerate(parts,1):
            root, rm = E.fromstring(z.read(part)), relationships(z,part)
            words = text_of(z,part)
            notes = ' '.join(text_of(z,r['part']) for r in rm.values() if r['kind']=='notesSlide' and not r['external'])
            if root.get('show') in ('0','false'):
                hidden.append(i)
            titles.append({'slide':i,'text':words,'notes':notes})
            for location,text in [('slide',words),('notes',notes)]:
                if re.search(r'\btransmission\s+(?:models?|paths?)\b',text,re.I):
                    rejected.append({'slide':i,'location':location,'text':text})
                if i<=MAIN and re.search(r'\bSADb\b',text,re.I):
                    sadb_main.append({'slide':i,'location':location,'text':text})
            for rel in rm.values():
                if not rel['external'] and rel['part'] not in z.namelist():
                    missing.append({'slide':i,'part':rel['part']})
            cps = []
            for c in root.findall('.//c:chart',N):
                rel = rm.get(c.get('{'+N['r']+'}id'))
                require(rel is not None and rel['kind']=='chart' and not rel['external'],f'Slide {i}: native chart relationship')
                if rel is not None:
                    cps.append(rel['part'])
            if cps:
                chart_owners[i] = cps
                attached_books = []
                for cp in cps:
                    cr, crm = E.fromstring(z.read(cp)), relationships(z,cp)
                    for ed in cr.findall('c:externalData',N):
                        rel = crm.get(ed.get('{'+N['r']+'}id'))
                        if rel and not rel['external'] and rel['part'].endswith('.xlsx'):
                            attached_books.append(rel['part'])
                if attached_books:
                    workbook_owners[i] = attached_books
            clips = linked_video_hashes(z,part)
            if clips:
                video_owners[i] = sorted(clips)
        charts = sorted(n for n in z.namelist() if re.fullmatch(r'ppt/charts/chart\d+\.xml',n))
        oldcharts = sorted(n for n in b.namelist() if re.fullmatch(r'ppt/charts/chart\d+\.xml',n))
        books = sorted(n for n in z.namelist() if n.startswith('ppt/embeddings/') and n.endswith('.xlsx'))
        oldbooks = sorted(n for n in b.namelist() if n.startswith('ppt/embeddings/') and n.endswith('.xlsx'))
        chart_equal = {n:n in b.namelist() and z.read(n)==b.read(n) for n in charts}
        book_equal = {n:n in b.namelist() and z.read(n)==b.read(n) for n in books}
        protected = [n for n in b.namelist() if n.startswith(('ppt/slideMasters/','ppt/slideLayouts/','ppt/theme/'))]
        templates_equal = {n:n in z.namelist() and z.read(n)==b.read(n) for n in protected}
        template_media = {r['part'] for p in protected if p.endswith('.xml')
                          for r in relationships(b,p).values()
                          if not r['external'] and r['part'].startswith('ppt/media/')}
        template_media_equal = {n:n in z.namelist() and z.read(n)==b.read(n) for n in template_media}
        layout_failures = [i for i,p in enumerate(parts,1)
                           if any(r['kind']=='slideLayout' and r['part'] not in b.namelist()
                                  for r in relationships(z,p).values())]
        reference_branding = branding_images(b,oldparts[1])
        branding_failures = [i for i,part in enumerate(parts[1:],2) if not reference_branding.issubset(branding_images(z,part))]
        title_unchanged = z.read(parts[0])==b.read(oldparts[0])
        videos = {sha(z.read(n)) for n in z.namelist() if n.lower().endswith('.mp4')}
        oldvideos = {sha(b.read(n)) for n in b.namelist() if n.lower().endswith('.mp4')}
        active_videos = {v for values in video_owners.values() for v in values}
        mapping = json.loads((W/'build/package_report.json').read_text(encoding='utf-8'))['source_to_final']
        expected_video_owners = {}
        for source,part in enumerate(oldparts,1):
            clips = linked_video_hashes(b,part)
            if clips and str(source) in mapping:
                expected_video_owners[int(mapping[str(source)])] = clips
        expected_active_videos = {v for values in expected_video_owners.values() for v in values}
        gifs = []
        for n in z.namelist():
            if n.lower().endswith('.gif'):
                with Image.open(io.BytesIO(z.read(n))) as im:
                    gifs.append({'part':n,'frames':im.n_frames})
    submitted = Path(r'C:\Users\Ben Bolen\Documents\ME603 Dissertation\Defense\Bolen_Dissertation.pdf')
    submitted_unchanged = sha(submitted.read_bytes())==SUBMITTED_SHA
    pdf = W/'output/Bolen_PhD_Defense_RoseGreen_20261102_v20.pdf'
    pdf_pages = None if args.skip_pdf else len(PdfReader(pdf).pages)
    require(len(parts)==TOTAL,'127 slides')
    require(hidden==list(range(MAIN+1,TOTAL+1)),'Main 1–74; hidden 75–127')
    require(not rejected,'No rejected transmission-model/path wording')
    require(not sadb_main,'SADb absent from main slide text and notes')
    require(not missing,'All slide relationships resolve')
    require(sorted(chart_owners)==CHART_OWNERS,'Exact native chart owners')
    require(sorted(workbook_owners)==CHART_OWNERS,'Exact embedded-workbook chart owners')
    require(charts==oldcharts and all(chart_equal.values()),'Native chart inventory and bytes preserved')
    require(books==oldbooks and all(book_equal.values()),'Embedded workbook inventory and bytes preserved')
    require(videos==oldvideos and active_videos==expected_active_videos,'Original video assets and actively used clips preserved')
    require(video_owners==expected_video_owners,'Native video shapes remain on their mapped original slides')
    require(all(templates_equal.values()),'Original masters, layouts, and themes preserved')
    require(not branding_failures and not layout_failures and bool(template_media_equal)
            and all(template_media_equal.values()),'Original logo/footer template media and layout preserved')
    require(title_unchanged,'Original title/honeycomb slide preserved')
    require(submitted_unchanged,'Submitted dissertation PDF unchanged')
    if not args.skip_pdf:
        require(pdf_pages==TOTAL,'Review PDF includes all 127 slides')
    report = {'deck':str(target),'sha256':sha(target.read_bytes()),'bytes':target.stat().st_size,
              'slides':len(parts),'main_slides':MAIN,'hidden_slides':hidden,'rejected_wording':rejected,
              'sadb_main':sadb_main,'missing_slide_relationships':missing,
              'native_chart_owners':chart_owners,'embedded_workbook_owners':workbook_owners,
              'native_charts_unchanged':chart_equal,'embedded_workbooks_unchanged':book_equal,
              'native_video_hashes_preserved':videos==oldvideos,'native_videos_slide_linked':active_videos==expected_active_videos,
              'video_owner_slides':video_owners,'videos':len(videos),'animated_gifs':gifs,
              'active_video_count':len(active_videos),'preserved_inactive_video_assets':len(videos-active_videos),
              'native_video_owner_mapping_preserved':video_owners==expected_video_owners,
              'original_template_parts_unchanged':templates_equal,'original_title_unchanged':title_unchanged,
              'branding_reference_image_count':len(reference_branding),'branding_mismatch_slides':branding_failures,
              'branding_basis':'Original masters/layouts, their image bytes, and unchanged title; direct edge images additionally checked when present.',
              'template_media_unchanged':template_media_equal,'nonreference_layout_slides':layout_failures,
              'submitted_pdf_unchanged':submitted_unchanged,'submitted_pdf_sha256':SUBMITTED_SHA,
              'review_pdf_pages':pdf_pages,'pdf_check_skipped':args.skip_pdf,'media_playback_verified':False,
              'scope':'Package preservation and static PDF count only; visual review, media playback, and timed rehearsal are separate.',
              'failures':failures,'passed':not failures}
    prefix = 'candidate' if target.name=='candidate.pptx' else 'delivery'
    (W/f'build/{prefix}_checks.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    if prefix=='delivery':
        (W/'build/final_slide_text.json').write_text(json.dumps(titles,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps({k:report[k] for k in ['sha256','bytes','slides','videos','review_pdf_pages','passed','failures']}))
    if failures:
        raise SystemExit(1)

if __name__=='__main__':
    main()
