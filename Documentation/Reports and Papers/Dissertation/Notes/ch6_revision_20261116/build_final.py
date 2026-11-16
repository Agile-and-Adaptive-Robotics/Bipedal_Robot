"""Compile the canonical Final tree without staging or packaging a ZIP.

Only generated build products are written below the repository's ignored tmp.
Use --isolated-miktex only for the scratch runtime described in TEX_RUNTIME_NOTE.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[5]
DISS = ROOT/'Documentation/Reports and Papers/Dissertation'
NOTES = DISS/'Notes/ch6_revision_20261116'
FINAL = DISS/'Final'
SCRATCH = ROOT/'tmp/ch6_20261116'
DEFAULT_TEXBIN = Path(r'C:\Users\Ben\AppData\Local\Programs\MiKTeX\miktex\bin\x64')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tex-bin', type=Path, default=DEFAULT_TEXBIN)
    parser.add_argument('--isolated-miktex', action='store_true')
    parser.add_argument('--revision', default='canonical_r1')
    parser.add_argument('--timeout', type=int, default=180)
    args = parser.parse_args()
    build = SCRATCH/f'build_final_{args.revision}'
    build.mkdir(parents=True, exist_ok=True)
    (build/'chapters').mkdir(exist_ok=True)
    env = os.environ.copy()
    env['BIBINPUTS'] = str(FINAL)+';'
    env['BSTINPUTS'] = str(FINAL)+';'
    flags = ['--miktex-disable-maintenance','--miktex-disable-installer']
    fmtargs = []
    owned_fmt = None
    if args.isolated_miktex:
        probe = SCRATCH/'tex_probe'
        env.update({'MIKTEX_USERCONFIG':str(probe/'config'),
                    'MIKTEX_USERDATA':str(probe/'data'),
                    'MIKTEX_USERINSTALL':str(args.tex_bin.parents[2]),
                    'MIKTEX_USERSTARTUPFILE':str(probe/'startup.ini'),
                    'MIKTEX_COMMONSTARTUPFILE':str(probe/'startup.ini')})
        cache = probe/'data/fonts'
        install = args.tex_bin.parents[2]/'fonts'
        env['TFMFONTS'] = ';'.join([str(cache/'tfm/public/mathabx'),str(cache/'tfm/public/knuth-lib'),str(install/'tfm')+'//',''])
        env['PKFONTS'] = ';'.join(str(p) for p in (cache/'pk').rglob('dpi*') if p.is_dir())+';'
        # Replace only the scratch-registered stale format. This avoids automatic
        # regeneration when an output directory suppresses current-dir lookup.
        scratch_fmt = probe/'data/miktex/data/le/pdftex/pdflatex.fmt'
        scratch_fmt.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(probe/'work/pdflatex.fmt',scratch_fmt)
    sources = list(FINAL.rglob('*.tex'))+list(FINAL.glob('*.bib'))
    before = {str(p.relative_to(FINAL)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    report = {'source_root':str(FINAL),'build':str(build),'status':'running',
              'isolated_runtime':args.isolated_miktex,'commands':[]}
    start = time.monotonic()
    def call(cmd,cwd,name):
        report['commands'].append([str(x) for x in cmd])
        with (build/name).open('w',encoding='utf-8') as out:
            try:
                proc = subprocess.run([str(x) for x in cmd],cwd=cwd,env=env,
                                      stdout=out,stderr=subprocess.STDOUT,timeout=args.timeout)
            except subprocess.TimeoutExpired:
                report.update(status='timeout',failed_step=name)
                raise
        if proc.returncode:
            report.update(status='failed',failed_step=name,exit_code=proc.returncode)
            raise RuntimeError(f'{name}: exit {proc.returncode}')
    latex = [args.tex_bin/'pdflatex.exe',*flags,*fmtargs,'-interaction=nonstopmode',
             '-halt-on-error','-file-line-error','-output-directory='+str(build),'main.tex']
    try:
        call(latex,FINAL,'pass1.txt')
        call([args.tex_bin/'bibtex.exe',*flags,'main'],build,'bibtex.txt')
        call(latex,FINAL,'pass2.txt')
        call(latex,FINAL,'pass3.txt')
        log = (build/'main.log').read_text(encoding='utf-8',errors='replace')
        gates = {'undefined_references': bool(re.search(r'There were undefined references|Reference .* undefined|Citation .* undefined',log)),
                 'duplicate_labels':'multiply defined' in log,
                 'missing_files':bool(re.search(r'File .* not found|LaTeX Error: File',log)),
                 'oversized_floats':'Float too large' in log,
                 'pending_rerun':bool(re.search(r'Label\(s\) may have changed|Rerun to get',log))}
        report['gates'] = gates
        report['status'] = 'gate_failed' if any(gates.values()) else 'compiled'
        if any(gates.values()):
            raise RuntimeError('Final LaTeX gate failed')
        from pypdf import PdfReader
        report['pages'] = len(PdfReader(build/'main.pdf').pages)
        delivery = ROOT/'output/pdf/Bolen_Dissertation_Final_Chapter6_20261116.pdf'
        delivery.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(build/'main.pdf',delivery)
        report['pdf'] = str(delivery)
        report['pdf_sha256'] = hashlib.sha256(delivery.read_bytes()).hexdigest()
    except KeyboardInterrupt:
        report.update(status='interrupted',error='Build interrupted before verification')
    except Exception as exc:
        report['error'] = str(exc)
        print(str(exc))
    finally:
        after = {str(p.relative_to(FINAL)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
        report['source_conserved'] = before==after
        report['elapsed_s'] = round(time.monotonic()-start,3)
        if owned_fmt is not None and owned_fmt.exists():
            if hashlib.sha256(owned_fmt.read_bytes()).hexdigest()==hashlib.sha256((SCRATCH/'tex_probe/work/pdflatex.fmt').read_bytes()).hexdigest():
                owned_fmt.unlink()
                report['temporary_format_removed'] = True
        (NOTES/'FINAL_BUILD_RESULT.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
        print(json.dumps(report,indent=2))
    if report['status']!='compiled':
        raise SystemExit(1)

if __name__=='__main__':
    main()
