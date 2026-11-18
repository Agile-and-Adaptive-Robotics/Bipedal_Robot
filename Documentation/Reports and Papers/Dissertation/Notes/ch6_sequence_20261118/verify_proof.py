"""Inspect the complete local and export-based builds and make the bookmarked reading copy."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import shutil
from zipfile import ZipFile

import pymupdf as pdf
from PIL import Image, ImageDraw


HERE = Path(__file__).resolve().parent
DISSERTATION = HERE.parents[1]
FIGURES = DISSERTATION/'Figures'
CHAPTERS = DISSERTATION/'Final/chapters'
AUDITED = {'60-sns', '70-overall', '92-AppendixA', '93-AppendixB', '94-AppendixC', '95-AppendixD', '96-AppendixE'}
GRAPHIC_EXTENSIONS = {'.pdf', '.png', '.jpg', '.jpeg', '.eps'}


def digest(data):
    return sha256(data).hexdigest()


def build_log(build):
    log=(build/'main.log').read_text(encoding='utf-8',errors='replace')
    blg=(build/'main.blg').read_text(encoding='utf-8',errors='replace')
    bad=[msg for msg in ('undefined on input line','There were undefined','multiply defined',
                          'multiply-defined','Fatal error','Overfull ','Missing character:',
                          'Float too large','Rerun to get','Label(s) may have changed')
         if msg.lower() in log.lower()]
    assert not bad, (build,bad)
    assert not re.search(r'^(?:!|LaTeX Error)',log,re.M)
    assert not re.search(r"Warning--|error message|I couldn't open",blg,re.I)
    return {'latex_passes':len(list(build.glob('latex*.stdout'))),'bibtex':(build/'main.bbl').stat().st_size>0,
            'pdftex_notices':log.count('pdfTeX warning'), 'latex_warnings':log.count('LaTeX Warning'),
            'overfull_boxes':0,'underfull_hboxes':log.count(r'Underfull \hbox'),
            'underfull_vboxes':log.count(r'Underfull \vbox')}


def used_assets(build,source,root):
    used={}
    for line in (build/'main.fls').read_text(encoding='utf-8',errors='replace').splitlines():
        if line.startswith('INPUT '):
            path=(source/line[6:].strip().strip('"')).resolve()
            if path.is_relative_to(root) and path.suffix.lower() in GRAPHIC_EXTENSIONS:
                rel=path.relative_to(root)
                assert len(rel.parts)==2,rel
                if rel.as_posix() not in used:
                    used[rel.as_posix()]=digest(path.read_bytes())
    return used


def aux_labels(build):
    labels={}
    for name in ('60-sns','70-overall','95-AppendixD'):
        data=(build/'chapters'/(name+'.aux')).read_text(encoding='utf-8')
        for key,number,page in re.findall(r'\\newlabel\{([^}]+)\}\{\{([^}]*)\}\{(\d+)\}',data):
            labels[key]=(number,int(page))
    return labels


def proof_files(build,render_dir,labels):
    render_dir.mkdir(exist_ok=True)
    manuscript=HERE/'Bolen_Dissertation_revised.pdf'
    shutil.copyfile(build/'main.pdf',manuscript)
    mapping={}
    with pdf.open(manuscript) as doc:
        toc=doc.get_toc()
        first=next(p for level,title,p in toc if level==1 and 'Tools for Constructing' in title)
        last=next(p for level,title,p in toc if level==1 and p>first)-1
        offset=first-labels['ch:sns'][1]
        def page_for(label):
            return labels[label][1]+offset
        conclusion=next(p for level,title,p in toc if first<=p<=last and re.search(r'\bConclusion$',title))
        key=[('The AnimatLab problem',first),
             ('The two routes to MuJoCo',page_for('fig:simulation_workflow')),
             ('Mechanical models in order',page_for('tab:sns_model_configurations')),
             ('Designing a large neural circuit',page_for('sec:sns_new_tools')),
             ('Stimulus and response on a common clock',page_for('fig:sns_reflex_current')),
             ('Retesting and optimizing the controller',page_for('fig:sns_design_fit_workflow')),
             ('Conclusion',conclusion)]
        final_phrase=next(n for n in range(conclusion,last+1) if 'future researcher' in doc[n-1].get_text().lower())
        assert final_phrase==conclusion==last
        assert 'change in force and motion.' in ' '.join(doc[conclusion-1].get_text().split())
        brief=pdf.open()
        bookmarks=[]
        for title,num in key:
            brief.insert_pdf(doc,from_page=num-1,to_page=num-1,links=False)
            bookmarks.append([1,title,len(brief)])
        brief.set_toc(bookmarks)
        brief.set_metadata({'title':'Chapter 6 reading pages','author':'Benjamin Bolen',
                            'subject':'The revised chapter opening, workflows, model table, knee experiment, and conclusion'})
        brief.save(HERE/'Chapter_6_reading_pages.pdf',garbage=4,deflate=True,no_new_id=True)
        brief.close()
        appendix_page=next(i+1 for i,p in enumerate(doc) if i+1>last and
                           '0.88 and 0.98' in ' '.join(p.get_text().split()))
        to_render=list(range(first,last+2))+[appendix_page]
        files=[]
        for n in to_render:
            page=doc[n-1]
            dest=render_dir/f'page_{n:03}.png'
            page.get_pixmap(dpi=126,alpha=False).save(dest)
            files.append(dest)
        for k in range(0,len(files),6):
            sheet=Image.new('RGB',(1200,1060),'#e5e7e8')
            draw=ImageDraw.Draw(sheet)
            for i,filename in enumerate(files[k:k+6]):
                with Image.open(filename) as im:
                    im.thumbnail((385,490))
                    x=(i%3)*400+7
                    y=(i//3)*530
                    draw.text((x+2,y+2),filename.stem,fill='black')
                    sheet.paste(im,(x,y+20))
            sheet.save(render_dir/f'sheet_{k//6+1:02}.jpg',quality=88)
        for name,n in key:
            assert 0<len(doc[n-1].get_text()) and first<=n<=last
            mapping[name]={'pdf_page':n,'printed_page':n-offset}
        return {'pages':len(doc),'chapter6_pdf_start':first,'chapter6_pdf_end':last,
                'physical_minus_arabic_number':offset,'appendix_correction_pdf_page':appendix_page,
                'reading_pages':mapping,'figure_6_4_pdf_page':page_for('fig:sns_reflex_current')}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('working_build','export_build','export_source','renders'):
        parser.add_argument(name,type=Path)
    args=parser.parse_args()
    local,assembled,source,renders=[p.resolve() for p in (args.working_build,args.export_build,args.export_source,args.renders)]
    builds={'working':build_log(local),'latest_export_with_replacements':build_log(assembled)}
    assets=used_assets(local,DISSERTATION/'Final',FIGURES)
    assert assets==used_assets(assembled,source,source/'figs')
    target={name:sha for name,sha in assets.items() if name.split('/')[0] in AUDITED}
    roots={path.relative_to(FIGURES).as_posix() for stem in AUDITED for path in (FIGURES/stem).glob('*')
           if path.is_file() and path.suffix.lower() in GRAPHIC_EXTENSIONS}
    assert roots==target.keys(), roots^target.keys()
    prior=json.loads((HERE.parent/'figure_cleanup_20261118/COMPILED_FIGURES.json').read_text(encoding='utf-8'))
    expected={r['chapter']+'/'+name for r in prior['rows'] for name in r['graphic']}
    assert target.keys()-expected=={'60-sns/controller_retest_and_optimization.pdf'} and not expected-target.keys()
    with ZipFile(DISSERTATION/'Bolen_Dissertation.zip') as z:
        original=z.read('chapters/60-sns.tex').decode('utf-8')
        original_d=z.read('chapters/95-AppendixD.tex')
        assert original_d.replace(b'between 0.87 and 0.98',b'between 0.88 and 0.98')==(CHAPTERS/'95-AppendixD.tex').read_bytes()
    current=(CHAPTERS/'60-sns.tex').read_text(encoding='utf-8')
    model_rows=['Li--Deng walker in AnimatLab','W2L walker in MuJoCo','Gait2392 human',
                'Gait2392 robot counterpart','CAD robot','Single-joint knee','Single-joint arm']
    start=current.index(r'\label{tab:sns_model_configurations}')
    first_table=current[start:current.index(r'\end{table}',start)]
    row_names=[line.split(' & ')[0] for line in first_table.splitlines() if ' & ' in line][1:]
    assert row_names==model_rows,row_names
    introductions=['the sagittal-plane biped','In a later walker model, W2L','gait2392_simbody.osim',
                   'gait2392_robotbody.osim','Building the Mechanical Model from Robot CAD',
                   'In the reduced knee','The loaded-cup model uses']
    positions=[current.index(phrase) for phrase in introductions]
    assert positions==sorted(positions),positions
    assert 'A future researcher' in current.split(r'\section{Conclusion}')[-1]
    pattern=r'\\label\{([^}]+)\}'
    old_labels=set(re.findall(pattern,original))
    new_labels=set(re.findall(pattern,current))
    assert not old_labels-new_labels,old_labels-new_labels
    def citations(s):
        return {key.strip() for group in re.findall(r'\\cite[pt]?\{([^}]+)\}',s) for key in group.split(',')}
    assert not citations(original)-citations(current)
    labels=aux_labels(assembled)
    assert labels['fig:sns_reflex_current'][0]=='6.4' and labels['tab:sns_model_configurations'][0]=='6.1'
    assert labels['fig:simulation_workflow'][0]=='6.2' and labels['fig:sns_design_fit_workflow'][0]=='6.6'
    layout=proof_files(assembled,renders,labels)
    with pdf.open(local/'main.pdf') as doc:
        working_pages=len(doc)
    with pdf.open(HERE/'Chapter_6_reading_pages.pdf') as doc:
        assert len(doc)==len(layout['reading_pages'])
    out={'builds':builds,'working_proof_pages':working_pages,'delivered_layout':layout,
         'table_6_1_and_narrative_order':model_rows,'conclusion_complete_on_one_page':True,
         'revised_chapter_labels':sorted(new_labels-old_labels),'old_chapter_label_count':len(old_labels),
         'old_citation_keys':sorted(citations(original)),'added_citation_keys':sorted(citations(current)-citations(original)),
         'all_dissertation_graphic_files':len(assets),'chapter6_7_appendices_active_files':target,
         'chapter6_7_appendices_numbered_figures':prior['numbered_figures']+1,
         'chapter6_7_appendices_placements':prior['placements']+1,
         'numbered_chapter6_figures':{k:list(v) for k,v in labels.items() if k.startswith('fig:') and v[0].startswith('6.')},
         'sha256_full_manuscript':digest((HERE/'Bolen_Dissertation_revised.pdf').read_bytes())}
    (HERE/'PROOF_VERIFICATION.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({**layout,'builds':builds,'graphic_files':len(target)},indent=2))


if __name__=='__main__':
    main()
