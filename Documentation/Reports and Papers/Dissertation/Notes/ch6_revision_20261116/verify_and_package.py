"""Build a Chapter 6 review from the authoritative downloaded source package.

Chapter 6, Appendix D, and their selected graphics are overlaid. One linked
Chapter 7 drive-description correction is carried forward from the existing
Final text. All other downloaded sources are preserved byte for byte.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import zipfile
from pypdf import PdfReader, PdfWriter
import pypdfium2 as pdfium

ROOT = Path(__file__).resolve().parents[5]
DISS = ROOT/'Documentation/Reports and Papers/Dissertation'
NOTES = DISS/'Notes/ch6_revision_20261116'
SCRATCH = ROOT/'tmp/ch6_20261116'
TEXBIN = Path(r'C:\Users\Ben\AppData\Local\Programs\MiKTeX\miktex\bin\x64')
OLD_DRIVE = ('Human-model stepping persists after removal of externally supplied neural drive, '
             'in place under the support rig; the robot model progresses under its harness with reduced drive')
NEW_DRIVE = ('In the human-model experiment, stepping persists after the supplementary left/right '
             'hand-drive channels are removed, in place under the support rig; the shared locomotor-drive '
             'input remains active. In the separate robot-model experiment, the robot model progresses '
             'under its harness with reduced drive')

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def run(cmd, cwd, log, env):
    with log.open('w', encoding='utf-8') as stream:
        proc = subprocess.run([str(x) for x in cmd], cwd=cwd, env=env,
                              stdout=stream, stderr=subprocess.STDOUT)
    if proc.returncode:
        tail = log.read_text(encoding='utf-8', errors='replace').splitlines()[-45:]
        print('\n'.join(tail))
        raise RuntimeError(f'Build command failed ({proc.returncode}): {cmd[0]}')

def prepare(stage):
    with zipfile.ZipFile(DISS/'Bolen_Dissertation.zip') as archive:
        archive.extractall(stage)
    for name in ('60-sns.tex', '95-AppendixD.tex'):
        shutil.copy2(DISS/'Final/chapters'/name, stage/'chapters'/name)
    for folder in ('60-sns', '95-AppendixD'):
        for src in (DISS/'Figures'/folder).glob('*20261116*'):
            if src.is_file():
                shutil.copy2(src, stage/'figs'/folder/src.name)
    overall = stage/'chapters/70-overall.tex'
    data = overall.read_bytes()
    if data.count(OLD_DRIVE.encode('utf-8'))!=1:
        raise RuntimeError('Chapter 7 drive correction needs reconciliation')
    overall.write_bytes(data.replace(OLD_DRIVE.encode('utf-8'), NEW_DRIVE.encode('utf-8')))
    return stage

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--revision', default='r1')
    parser.add_argument('--prepare-only', action='store_true')
    parser.add_argument('--tex-bin', type=Path, default=TEXBIN)
    parser.add_argument('--isolated-miktex', action='store_true')
    args = parser.parse_args()
    stage = SCRATCH/f'review_source_{args.revision}'
    build = SCRATCH/f'build_{args.revision}'
    stage.mkdir(parents=True, exist_ok=True)
    build.mkdir(parents=True, exist_ok=True)
    (build/'chapters').mkdir(exist_ok=True)
    prepare(stage)
    if args.prepare_only:
        print(stage)
        return
    env = os.environ.copy()
    env['BIBINPUTS'] = str(stage)+';'
    env['BSTINPUTS'] = str(stage)+';'
    if args.isolated_miktex:
        probe = SCRATCH/'tex_probe'
        env.update({'MIKTEX_USERCONFIG':str(probe/'config'),
                    'MIKTEX_USERDATA':str(probe/'data'),
                    'MIKTEX_USERINSTALL':str(args.tex_bin.parents[2]),
                    'MIKTEX_USERSTARTUPFILE':str(probe/'startup.ini'),
                    'MIKTEX_COMMONSTARTUPFILE':str(probe/'startup.ini')})
        shutil.copy2(probe/'work/pdflatex.fmt', stage/'pdflatex.fmt')
    flags = ['--miktex-disable-maintenance','--miktex-disable-installer']
    latex = [args.tex_bin/'pdflatex.exe', *flags, '-interaction=nonstopmode', '-halt-on-error',
             '-file-line-error', f'-output-directory={build}', 'main.tex']
    run(latex, stage, build/'pass1.txt', env)
    run([args.tex_bin/'bibtex.exe', *flags, 'main'], build, build/'bibtex_pass.txt', env)
    for n in (2,3,4):
        run(latex, stage, build/f'pass{n}.txt', env)
    log = (build/'main.log').read_text(encoding='utf-8', errors='replace')
    gates = {
        'undefined_references': bool(re.search(r'There were undefined references|Reference .* undefined|Citation .* undefined', log)),
        'duplicate_labels': 'multiply defined' in log,
        'missing_files': bool(re.search(r'File .* not found|LaTeX Error: File', log)),
        'oversized_floats': 'Float too large' in log,
        'pending_rerun': bool(re.search(r'Label\(s\) may have changed|Rerun to get', log)),
    }
    if any(gates.values()):
        print(json.dumps(gates, indent=2))
        raise RuntimeError('A final LaTeX gate failed')
    reader = PdfReader(build/'main.pdf')
    chapter_pages = []
    appendix_pages = []
    for n, page in enumerate(reader.pages):
        text = page.extract_text() or ''
        if re.search(r'Chapter 6\b', text) or ('Tools for Constructing' in text and n>20):
            if 'Contents' not in text:
                chapter_pages.append(n)
        if 'Appendix D:' in text and 'Supplementary' in text:
            appendix_pages.append(n)
    # Include chapter/appendix pages from their .aux labels, which are stable even
    # when extraction of a multi-line chapter heading changes between TeX versions.
    mainaux = (build/'main.aux').read_text(errors='replace')
    aux = '\n'.join(p.read_text(errors='replace') for p in (build/'chapters').glob('*.aux'))
    labels = {}
    for label, printed in re.findall(r'\\newlabel\{([^}]+)\}\{\{[^}]*\}\{([^}]+)\}', aux):
        labels[label] = printed
    output = ROOT/'output/pdf'
    output.mkdir(parents=True, exist_ok=True)
    pdf_out = output/'Bolen_Dissertation_Chapter6_Review_20261116.pdf'
    shutil.copy2(build/'main.pdf', pdf_out)
    sources = ROOT/'output/packages'
    sources.mkdir(parents=True, exist_ok=True)
    zip_out = sources/'Bolen_Dissertation_Chapter6_Review_20261116.zip'
    # Keep only referenced graphics in the delivery, retaining matched EPS/PDF
    # conversion caches when present. The complete baseline remains untouched.
    active = set()
    missing = []
    for tex in (stage/'chapters').glob('*.tex'):
        text = tex.read_text(encoding='utf-8', errors='replace')
        paths = re.findall(r'\{(figs/[^{}]+/)\}', text)
        for graphic in re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^{}]+)\}', text):
            candidates = [stage/p/graphic for p in paths]+[stage/'figs'/graphic, stage/graphic]
            found = []
            for candidate in candidates:
                if candidate.suffix:
                    if candidate.exists(): found.append(candidate)
                else:
                    found.extend(candidate.with_suffix(s) for s in ('.pdf','.png','.jpg','.jpeg','.eps') if candidate.with_suffix(s).exists())
            if not found:
                missing.append((tex.name, graphic))
            for item in found:
                active.add(item)
                if item.suffix.lower()=='.eps':
                    cache = item.with_name(item.stem+'-eps-converted-to.pdf')
                    if cache.exists(): active.add(cache)
    if missing:
        raise RuntimeError(f'Cannot package unresolved graphics: {missing}')
    with zipfile.ZipFile(zip_out, 'w', zipfile.ZIP_DEFLATED) as archive:
        for path in stage.rglob('*'):
            if not path.is_file(): continue
            rel = path.relative_to(stage)
            if rel.parts[0]=='figs':
                if path not in active: continue
            elif path.suffix.lower() not in ('.tex','.cls','.bst','.bib'):
                continue
            archive.write(path, rel.as_posix())
    unchanged = []
    with zipfile.ZipFile(DISS/'Bolen_Dissertation.zip') as baseline:
        for member in baseline.namelist():
            if Path(member).suffix.lower() in ('.tex','.bib','.bst','.cls') and member not in ('chapters/60-sns.tex','chapters/95-AppendixD.tex'):
                expected = baseline.read(member)
                if member=='chapters/70-overall.tex':
                    expected = expected.decode('utf-8').replace(OLD_DRIVE, NEW_DRIVE).encode('utf-8')
                actual = (stage/member).read_bytes()
                if expected.replace(b'\r\n',b'\n')!=actual.replace(b'\r\n',b'\n'):
                    raise RuntimeError(f'Unrelated source changed: {member}')
                if member!='chapters/70-overall.tex':
                    unchanged.append(member)
    report = {'revision':args.revision, 'build':str(build.relative_to(ROOT)),
              'pages':len(reader.pages), 'gates':gates,
              'overfull_box_count':log.count('Overfull \\hbox'),
              'underfull_box_count':log.count('Underfull \\hbox'),
              'unchanged_source_files':unchanged,
              'bounded_consistency_correction':{'file':'chapters/70-overall.tex','old':OLD_DRIVE,'new':NEW_DRIVE},
              'input_zip_sha256':sha(DISS/'Bolen_Dissertation.zip'),
              'output_pdf':str(pdf_out.relative_to(ROOT)), 'output_pdf_sha256':sha(pdf_out),
              'output_zip':str(zip_out.relative_to(ROOT)), 'output_zip_sha256':sha(zip_out),
              'packaged_graphics':len(active), 'labels':labels}
    (NOTES/'BUILD_RESULT.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('labels','unchanged_source_files')}, indent=2))

if __name__=='__main__':
    main()
