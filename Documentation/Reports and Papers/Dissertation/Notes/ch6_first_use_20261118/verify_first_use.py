"""Check first uses in Chapter 6, verify the assembled proof, and extract the reading pages."""
import argparse
from collections import Counter
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import re
import shutil
import sys
import unicodedata
from zipfile import ZipFile

import pymupdf
from PIL import Image, ImageDraw


HERE = Path(__file__).resolve().parent
NOTES = HERE.parent
DISS = NOTES.parent
PRIOR = NOTES / 'ch6_sequence_20261118'
CHAPTER = DISS / 'Final/chapters/60-sns.tex'
PACKAGE = 'Chapter_6_and_Appendix_D_Overleaf.zip'
sys.path.insert(0, str(NOTES/'advisor_feedback_history'))
from review import passages


# The words after each name identify its job or distinguish it from a related model.
# The search includes captions and table cells. Section titles can name the topic.
TERMS = {
    'AnimatLab': ('graphical program', 'couples a neural circuit to a moving body'),
    'CPG': ('central pattern generator', 'produces rhythmic motor commands'),
    'Li--Deng': ('sagittal-plane biped', 'two-level central pattern generator'),
    'RG': ('rhythm generator', 'set the period'),
    'PF': ('Pattern-formation', 'determined when motor neurons'),
    'MNs': ('motor neurons', 'recruited each muscle group'),
    'Hill-type': ('contractile and elastic elements',),
    'MuJoCo': ('Multi-Joint dynamics with Contact', 'physics engine', 'joint motion and ground contact forces'),
    'SNS-Toolbox': ('Python library', 'neurons and synapses'),
    'W2L': ('Walker Two-Layer', 'each leg its own rhythm generator', 'heel and toe contact'),
    'OpenSim': ('software for modeling how muscles move the skeleton',),
    'Gait2392': ('human model', '23 degrees of freedom', '92 musculotendon actuators'),
    'BPAs': ('braided pneumatic actuators', 'shorten under pressure'),
    'SOLIDWORKS': ('computer-aided design',),
    'SNS-Simscape': ('I built', 'library of neural and actuator blocks'),
    'Simulink': ('block-diagram simulator',),
    'MathWorks bridge': ('sends actuator commands', 'returns the resulting motion and forces'),
    'BPA20': ('diameter BPAs', '20'),
    'MyoConverter': ('translates OpenSim muscle paths and joint geometry',),
    'static optimization': ('estimates how muscles share the load',),
    'Simscape Multibody': ('library that represents rigid bodies and joints',),
    'W2Lvar': ('Gait2392 human model', 'coordinates the hip muscles', 'coordinates the knee and ankle muscles'),
    'motor pool': ('group of motor neurons',),
    'SNSGui': ('circuit editor',),
    'LIF': ('leaky integrate-and-fire', 'emits a spike at threshold'),
    'Ia': ('muscle spindles', 'respond to stretch'),
    'Ib': ('tendon organs', 'respond to tension'),
    'SADb': ('Sensory Afferent Database', 'which sensory nerve fibers an experiment stimulated'),
    'NumPy': ('numerical array library for Python',),
    'deafferented': ('sensory return to the controller is disabled',),
    'afferented': ('means it is enabled',),
    'BiW2Lvar+': ('muscles spanning both joints', "Hunt's method", 'suspended Gait2392 human model'),
    's3k': ('reference controller', 'four pattern-formation cells per leg'),
    'syn6': ('six-synergy controller', 'six channels per leg', 'group of muscles with fitted weights'),
    'W2L108': ('W2L controller','twelve-muscle walker', '108 neurons and 216 connections', 'moving treadmill'),
    'Hunt-PT+RC': ('balance pathway', 'center of mass','Renshaw cells', 'inhibition back to the same motor pool'),
}


def token(word):
    return re.compile(r'(?<![a-zA-Z0-9_])' + re.escape(word) + r'(?![a-zA-Z0-9_])', re.I)


def checked_source():
    protected = re.compile(r'\\(?:label|ref|cite[tp]?|path)\{[^}]*\}')
    entries = [(start, protected.sub('', content)) for start, _, _, content in passages(CHAPTER)]
    found = {}
    for name, context in TERMS.items():
        match = token(name)
        start, content = next((n, s) for n, s in entries if match.search(s))
        missing = [phrase for phrase in context if phrase.casefold() not in content.casefold()]
        assert not missing, (name, start, missing)
        found[name] = {'definition_tex_line': start}
    return found


def plain(text):
    s = unicodedata.normalize('NFKC', text).replace('\u00ad', '').replace('–', '-')
    s = re.sub(r'(?<=\w)-\s*\n\s*(?=\w)', '', s)
    return ' '.join(s.split())


def read_helpers():
    spec = importlib.util.spec_from_file_location('previous_proof', PRIOR/'verify_proof.py')
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    return helper


def compare_sources(assembled):
    def hash_bytes(data):
        return sha256(data).hexdigest()
    with ZipFile(PRIOR/PACKAGE) as earlier, ZipFile(HERE/PACKAGE) as updated, ZipFile(DISS/'Bolen_Dissertation.zip') as incoming:
        chapter = 'chapters/60-sns.tex'
        assert updated.testzip() is None and set(earlier.namelist()) == set(updated.namelist())
        member_count = len(updated.namelist())
        changed = [n for n in updated.namelist() if earlier.read(n) != updated.read(n)]
        assert changed == [chapter] and updated.read(chapter) == CHAPTER.read_bytes()
        original = earlier.read(chapter).decode()
        original_paths = {item.filename for item in incoming.infolist() if not item.is_dir()}
        intact = original_paths-set(updated.namelist())
        assert all((assembled/n).read_bytes() == incoming.read(n) for n in intact)
        assert all((assembled/n).read_bytes() == updated.read(n) for n in updated.namelist())
        protected = {name:hash_bytes(incoming.read(name)) for name in sorted(intact)
                     if name.startswith('chapters/') and name.endswith('.tex')}
    current = CHAPTER.read_text(encoding='utf-8')
    for pattern in (r'\\label\{[^}]+\}', r'\\includegraphics\[[^]]+\]\{[^}]+\}'):
        assert Counter(re.findall(pattern, original)) == Counter(re.findall(pattern, current))
    def cited(text):
        return {item for key in re.findall(r'\\cite[pt]?\{([^}]+)\}', text) for item in key.split(',')}
    assert cited(original) <= cited(current)
    units = re.compile(r'\\(?:qtyrange|qty|ang|num)\{[^}]*\}(?:\{(?:[^{}]|\{[^}]*\})*\}){0,2}')
    before = Counter(units.findall(original))
    after = Counter(units.findall(current))
    assert not (before-after), before-after
    assert after-before == {r'\qty{20}{\milli\meter}': 1}, after-before
    return {'members_in_transfer': member_count, 'only_changed_transfer_members': changed,
            'members_preserved_from_incoming_export': len(intact), 'preserved_chapter_sha256': protected,
            'added_measurement_notation': list((after-before).elements()),
            'additional_citation_keys': sorted(cited(current)-cited(original)),
            'source_sha256': hash_bytes(CHAPTER.read_bytes()),'transfer_sha256':hash_bytes((HERE/PACKAGE).read_bytes())}


def visual_pages(doc, numbers, destination):
    destination.mkdir(exist_ok=True, parents=True)
    for i, number in enumerate(numbers):
        if i % 6 == 0:
            sheet = Image.new('RGB', (1200, 1060), '#e5e7e8')
            draw = ImageDraw.Draw(sheet)
        path = destination/f'page_{number:03}.png'
        doc[number-1].get_pixmap(dpi=126, alpha=False).save(path)
        with Image.open(path) as im:
            im.thumbnail((385, 490))
            x, y = i % 3 * 400 + 7, (i % 6)//3*530
            draw.text((x+2,y+2),f'PDF page {number}',fill='black')
            sheet.paste(im,(x,y+20))
        if i % 6 == 5 or i == len(numbers)-1:
            sheet.save(destination/f'sheet_{i//6+1:02}.jpg',quality=90)


def check_proof(build, labels, found, renders):
    filename = HERE/'Bolen_Dissertation_revised.pdf'
    shutil.copyfile(build/'main.pdf', filename)
    with pymupdf.open(filename) as doc:
        toc = doc.get_toc()
        first = next(n for lev, title, n in toc if lev == 1 and 'Tools for Constructing' in title)
        last = next(n for lev, _, n in toc if lev == 1 and n > first)-1
        offset = first-labels['ch:sns'][1]
        texts = {n:plain(doc[n-1].get_text(sort=True)) for n in range(first,last+1)}
        for name, record in found.items():
            # TeX prints double hyphens as en dashes. Boundaries distinguish W2L
            # from W2Lvar, and LIF from ordinary words such as "simplified".
            term = name.replace('--','-')
            match = token(term)
            pages = [n for n, data in texts.items() if match.search(data)]
            assert pages, name
            record.update(first_pdf_page=pages[0], first_printed_page=pages[0]-offset)
        def figure(key):
            return labels[key][1]+offset
        assert max(found[name]['first_pdf_page'] for name in ('W2L','MuJoCo','SNS-Toolbox','Gait2392','SNS-Simscape')) <= figure('fig:simulation_workflow')
        assert found['W2Lvar']['first_pdf_page'] <= figure('tab:sns_capability')
        assert found['BiW2Lvar+']['first_pdf_page'] <= figure('fig:sns_biw2l_air')
        assert found['Hunt-PT+RC']['first_pdf_page'] <= figure('tab:sns_runtime')
        for name, later in [('s3k','None of the twelve'), ('syn6','None of the twelve'), ('W2L108','eight qualifying left')]:
            page = found[name]['first_pdf_page']
            tail = ' '.join(texts[n] for n in range(page, min(page+2,last)+1))
            assert tail.index(name) < tail.index(later), name
        conclusion = next(n for _, title, n in toc if first <= n <= last and re.search(r'\bConclusion$',title))
        assert conclusion == last and 'A future researcher' in texts[conclusion] and 'force and motion.' in texts[conclusion], (conclusion,last)
        wanted = [('The AnimatLab and Li–Deng walker','AnimatLab'),('MuJoCo, SNS-Toolbox, and Walker Two-Layer','W2L'),
                  ('OpenSim, Gait2392, and the CAD workflow','Gait2392'),('W2Lvar on the human model','W2Lvar'),
                  ('BiW2Lvar+ and its suspended experiment','BiW2Lvar+'),('The s3k and syn6 controller comparison','syn6'),
                  ('W2L108 on the moving treadmill','W2L108')]
        with pymupdf.open() as brief:
            bookmarks = []
            selected = []
            for title, name in wanted:
                n = found[name]['first_pdf_page']
                if n not in selected:
                    brief.insert_pdf(doc,from_page=n-1,to_page=n-1,links=False)
                    selected.append(n)
                bookmarks.append([1,title,selected.index(n)+1])
            brief.set_toc(bookmarks)
            brief.set_metadata({'title':'Chapter 6 • model and software introductions','author':'Benjamin Bolen'})
            brief.save(HERE/'Chapter_6_terms_reading_pages.pdf',deflate=True,garbage=4,no_new_id=True)
        with pymupdf.open(HERE/'Chapter_6_terms_reading_pages.pdf') as brief:
            assert len(brief) == len(selected)
            assert all(plain(brief[i].get_text()) == plain(doc[n-1].get_text()) for i,n in enumerate(selected))
        visual_pages(doc,list(range(first,last+2)),renders)
        return {'pages':len(doc),'chapter_pdf_pages':[first,last], 'chapter_printed_pages':[first-offset,last-offset],
                'key_reading_printed_pages':[n-offset for n in selected], 'conclusion_printed_page':conclusion-offset,
                'figure_6_4_pdf_page':figure('fig:sns_reflex_current'), 'terms':found,
                'pdf_sha256':sha256(filename.read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ('working_build','export_build','assembled_source','renders'):
        parser.add_argument(arg,type=Path)
    args = parser.parse_args()
    working, export, source, renders = [x.resolve() for x in vars(args).values()]
    found = checked_source()
    result = compare_sources(source)
    helper = read_helpers()
    result['builds'] = {'working':helper.build_log(working), 'export':helper.build_log(export)}
    assets = helper.used_assets(working,DISS/'Final',DISS/'Figures')
    assert assets == helper.used_assets(export,source,source/'figs')
    chapter_assets = {p:h for p,h in assets.items() if p.split('/')[0] in helper.AUDITED}
    older = json.loads((PRIOR/'PROOF_VERIFICATION.json').read_text(encoding='utf-8'))
    assert chapter_assets == older['chapter6_7_appendices_active_files']
    active_roots = {p.relative_to(DISS/'Figures').as_posix() for stem in helper.AUDITED
                    for p in (DISS/'Figures'/stem).iterdir() if p.is_file() and p.suffix.lower() in helper.GRAPHIC_EXTENSIONS}
    assert active_roots == chapter_assets.keys()
    labels = helper.aux_labels(export)
    assert labels['fig:sns_reflex_current'][0]=='6.4' and labels['tab:sns_model_configurations'][0]=='6.1'
    for label,(number,_) in older['numbered_chapter6_figures'].items():
        assert labels[label][0] == number
    result.update(all_compiled_graphics=len(assets), audited_chapter_graphics=len(chapter_assets),
                  proof=check_proof(export,labels,found,renders))
    with pymupdf.open(working/'main.pdf') as doc:
        result['working_proof_pages'] = len(doc)
    (HERE/'PROOF_VERIFICATION.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(len(found),'first uses,',result['proof']['pages'],'pages, reading pages',result['proof']['key_reading_printed_pages'])


if __name__ == '__main__':
    main()
