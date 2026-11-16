"""Build canonical Final/main.tex without changing or packaging its sources.

Run with D:\\Anaconda\\python.exe. Generated products stay outside Final.
No package installation, registry changes, or system-format regeneration is used.
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
DISS = ROOT / 'Documentation/Reports and Papers/Dissertation'
FINAL = DISS / 'Final'
NOTES = Path(__file__).resolve().parent
TEXBIN = Path(r'D:\MiKTeX\miktex\bin\x64')


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_hashes():
    return {str(p.relative_to(FINAL)): sha256(p)
            for p in sorted(FINAL.rglob('*')) if p.is_file()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--revision', default='r1')
    parser.add_argument('--timeout', type=int, default=600)
    parser.add_argument('--portable-probe', action='store_true',
                        help='Use only the scratch-local MiKTeX startup configuration.')
    parser.add_argument('--priority-log',type=Path,
                        help='Prioritize only installed directories observed in this prior build log.')
    parser.add_argument('--no-recursive-inputs',action='store_true',
                        help='Use only observed installed package directories, failing explicitly if a needed file is absent.')
    parser.add_argument('--cache-installed-support',action='store_true',
                        help='Use byte-identical scratch copies of observed installed support files to avoid repeated directory searches.')
    parser.add_argument('--private-install',action='store_true',
                        help='Use fully scratch-local core roots with byte-identical installed MiKTeX configuration copies.')
    args = parser.parse_args()
    build = ROOT / 'tmp/ch6_authoritative_20261116' / ('build_final_' + args.revision)
    if build.exists():
        raise SystemExit('Use a fresh --revision; build directory already exists: ' + str(build))
    (build / 'chapters').mkdir(parents=True)
    env = os.environ.copy()
    runtime_info = None
    env['BIBINPUTS'] = str(FINAL) + ';'
    env['BSTINPUTS'] = str(FINAL) + ';'
    if args.portable_probe:
        probe = build / 'runtime'
        (probe/'config/miktex/config').mkdir(parents=True)
        (probe/'data/miktex/data/le/pdftex').mkdir(parents=True)
        (probe/'startup.ini').write_text('[Auto]\nConfig=Portable\n[Paths]\n'
            'UserConfig='+str(probe/'config').replace('\\','/')+'\n'
            'UserData='+str(probe/'data').replace('\\','/')+'\n'
            'UserInstall=D:/MiKTeX\n',encoding='ascii')
        (probe/'config/miktex/config/miktex.ini').write_text(
            '[TeXandFriends]\nRenewFormatsOnUpdate=false\n',encoding='ascii')
        existing_data=Path(r'C:\Users\Ben Bolen\AppData\Local\MiKTeX\miktex\data\le')
        for database in existing_data.glob('*.fndb-5*'):
            shutil.copy2(database,probe/'data/miktex/data/le'/database.name)
        original_format=Path(r'C:\Users\Ben Bolen\AppData\Local\MiKTeX\miktex\data\le\pdftex\pdflatex.fmt')
        scratch_format=probe/'data/miktex/data/le/pdftex/pdflatex.fmt'
        scratch_format.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(original_format,scratch_format)
        os.utime(scratch_format,None)
        env.update(MIKTEX_USERCONFIG=str(probe / 'config'),
                   MIKTEX_USERDATA=str(probe / 'data'),
                   MIKTEX_USERINSTALL=r'D:\MiKTeX',
                   MIKTEX_USERSTARTUPFILE=str(probe / 'startup.ini'),
                   MIKTEX_COMMONSTARTUPFILE=str(probe / 'startup.ini'),
                   TEXINPUTS=r'.;D:\MiKTeX\tex\latex\base;D:\MiKTeX\tex\generic\config;D:\MiKTeX\tex//;',
                   TFMFONTS=r'C:\Users\Ben Bolen\AppData\Local\MiKTeX\fonts\tfm//;D:\MiKTeX\fonts\tfm//;',
                   TEXFONTMAPS=r'C:\Users\Ben Bolen\AppData\Local\MiKTeX\fonts\map\pdftex;D:\MiKTeX\fonts\map//;',
                   T1FONTS=r'D:\MiKTeX\fonts\type1//;',
                   ENCFONTS=r'D:\MiKTeX\fonts\enc//;',
                   VFFONTS=r'D:\MiKTeX\fonts\vf//;',
                   TEXFORMATS=str(scratch_format.parent))
        runtime_info={'format_source':str(original_format),'format_sha256':sha256(original_format),
                      'scratch_format_sha256':sha256(scratch_format),
                      'scratch_timestamp_updated':True,
                      'installed_files_changed':False,
                      'explicit_environment':{key:env[key] for key in (
                          'MIKTEX_USERCONFIG','MIKTEX_USERDATA','MIKTEX_USERINSTALL',
                          'MIKTEX_USERSTARTUPFILE','MIKTEX_COMMONSTARTUPFILE','TEXINPUTS',
                          'TFMFONTS','TEXFONTMAPS','T1FONTS','ENCFONTS','VFFONTS','TEXFORMATS')}}
    before = source_hashes()
    priority_info=None
    if args.priority_log:
        from BUILD_texpaths import texinputs
        env['TEXINPUTS'],priority_info=texinputs(FINAL,args.priority_log,not args.no_recursive_inputs)
        if runtime_info: runtime_info['explicit_environment']['TEXINPUTS']=env['TEXINPUTS']
    support_cache=None
    if args.cache_installed_support:
        if not args.portable_probe or not args.priority_log:
            raise SystemExit('--cache-installed-support requires --portable-probe and --priority-log')
        from BUILD_support_cache import prepare
        texcache,fontcache,support_cache=prepare(FINAL,args.priority_log,probe)
        env['TEXINPUTS']=str(FINAL)+';.;'+str(texcache)
        for key in ('TFMFONTS','VFFONTS','T1FONTS','TEXFONTMAPS','ENCFONTS'):
            env[key]=str(fontcache)
        runtime_info['explicit_environment'].update({key:env[key] for key in (
            'TEXINPUTS','TFMFONTS','VFFONTS','T1FONTS','TEXFONTMAPS','ENCFONTS')})
    if args.private_install:
        if not args.cache_installed_support:
            raise SystemExit('--private-install requires --cache-installed-support')
        install=probe/'private_install'
        installed_config=Path(r'D:\MiKTeX\miktex\config')
        shutil.copytree(installed_config,install/'miktex/config')
        config_records=[]
        for original in installed_config.rglob('*'):
            if original.is_file():
                target=install/'miktex/config'/original.relative_to(installed_config)
                original_sha=sha256(original)
                if sha256(target)!=original_sha: raise RuntimeError('Core configuration cache differs: '+str(original))
                config_records.append({'source':str(original),'cache':str(target),'sha256':original_sha})
        config_manifest=probe/'private_install_config_manifest.json'
        config_manifest.write_text(json.dumps(config_records,indent=2),encoding='utf-8')
        roots={'UserInstall':install,'CommonInstall':install,'UserConfig':probe/'config',
               'CommonConfig':probe/'config','UserData':probe/'data','CommonData':probe/'data'}
        (probe/'startup.ini').write_text('[Auto]\nConfig=Portable\n[Paths]\n'+''.join(
            key+'='+str(value).replace('\\','/')+'\n' for key,value in roots.items()),encoding='ascii')
        env.update(MIKTEX_USERINSTALL=str(install),MIKTEX_COMMONINSTALL=str(install),
                   MIKTEX_COMMONCONFIG=str(probe/'config'),MIKTEX_COMMONDATA=str(probe/'data'))
        env['BIBINPUTS']=str(FINAL)
        env['BSTINPUTS']=str(FINAL)+';'+str(texcache)
        os.utime(scratch_format,None)
        runtime_info.update(private_install=True,private_config_manifest=str(config_manifest),
                            private_config_manifest_sha256=sha256(config_manifest))
        runtime_info['explicit_environment'].update({key:env[key] for key in (
            'MIKTEX_USERINSTALL','MIKTEX_COMMONINSTALL','MIKTEX_COMMONCONFIG','MIKTEX_COMMONDATA')})
    result = {'canonical_entry': str(FINAL / 'main.tex'),
              'build_directory': str(build), 'status': 'running',
              'portable_probe': args.portable_probe,
              'scratch_runtime':runtime_info,
              'package_search_priority':priority_info,
              'support_cache':support_cache,
              'bibliography_environment': {k: env[k] for k in ('BIBINPUTS', 'BSTINPUTS')},
              'commands': [], 'sources_before': before}
    started = time.monotonic()
    flags = ['--miktex-disable-maintenance', '--miktex-disable-installer']
    latex = [str(TEXBIN / 'pdflatex.exe'), *flags, '-interaction=nonstopmode',
             '-halt-on-error', '-file-line-error', '-output-directory=' + str(build), 'main.tex']
    def run(command, cwd, name):
        if args.portable_probe and Path(command[0]).name.lower()=='pdflatex.exe':
            # Writable FNDB/core caches may be updated by the prior TeX/BibTeX
            # process; keep only this byte-identical scratch format newer so
            # MiKTeX does not launch unneeded automatic format generation.
            os.utime(scratch_format,None)
        entry = {'command': command, 'cwd': str(cwd), 'log': str(build / name)}
        result['commands'].append(entry)
        tick = time.monotonic()
        with (build / name).open('w', encoding='utf-8') as stream:
            proc = subprocess.Popen(command, cwd=cwd, env=env, stdout=stream,
                                    stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)
            entry['pid'] = proc.pid
            try:
                code = proc.wait(timeout=args.timeout)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait(timeout=10)
                entry['timed_out'] = True
                raise RuntimeError(name + ': bounded timeout; owned process terminated')
        entry.update(exit_code=code, elapsed_s=round(time.monotonic() - tick, 3))
        if code:
            raise RuntimeError(name + ': exit ' + str(code))
    try:
        run(latex, FINAL, 'pass1.txt')
        run([str(TEXBIN / 'bibtex.exe'), *flags, 'main'], build, 'bibtex.txt')
        run(latex, FINAL, 'pass2.txt')
        run(latex, FINAL, 'pass3.txt')
        log = (build / 'main.log').read_text(encoding='utf-8', errors='replace')
        def gates():
            return {'errors': bool(re.search(r'^!|LaTeX Error:|Emergency stop|Fatal error', log, re.M)),
                    'undefined_references_or_citations': bool(re.search(
                        r'There were undefined references|Reference .* undefined|Citation .* undefined', log)),
                    'duplicate_labels': 'multiply defined' in log,
                    'missing_graphics_or_files': bool(re.search(r'File .* not found|LaTeX Error: File', log)),
                    'oversized_floats': 'Float too large' in log,
                    'rerun_required': bool(re.search(r'Label\(s\) may have changed|Rerun to get', log))}
        result['gates'] = gates()
        if result['gates']['rerun_required']:
            run(latex, FINAL, 'pass4.txt')
            log = (build / 'main.log').read_text(encoding='utf-8', errors='replace')
            result['gates'] = gates()
        result['box_notices'] = re.findall(r'(?:Over|Under)full [^\n]+', log)
        if any(result['gates'].values()):
            raise RuntimeError('Compiled output failed one or more final build gates')
        import fitz
        doc = fitz.open(build / 'main.pdf')
        result.update(status='compiled', pages=len(doc), pdf=str(build / 'main.pdf'),
                      pdf_sha256=sha256(build / 'main.pdf'))
        aux = (build / 'main.aux').read_text(encoding='utf-8', errors='replace')
        toc = (build / 'main.toc').read_text(encoding='utf-8', errors='replace')
        (NOTES / 'BUILD_TOC.txt').write_text(toc, encoding='utf-8')
        result['chapter_aux_files'] = sorted(str(p) for p in (build / 'chapters').glob('*.aux'))
        result['target_label_pages'] = {}
        for name in ('60-sns', '95-AppendixD'):
            targetaux = (build / 'chapters' / (name + '.aux')).read_text(encoding='utf-8', errors='replace')
            result['target_label_pages'][name] = re.findall(r'\\newlabel\{([^}]+)\}\{\{[^}]*\}\{([^}]+)\}', targetaux)
        doc.close()
    except KeyboardInterrupt:
        result.update(status='interrupted',error='Stopped owned attempt to correct runtime lookup paths; not an accepted build')
    except Exception as exc:
        result.update(status='blocked' if not (build / 'main.pdf').exists() else 'failed',
                      error=str(exc))
        logs = '\n'.join(p.read_text(encoding='utf-8', errors='replace')
                         for p in build.glob('*.txt'))
        result['exact_failure_lines'] = [line for line in logs.splitlines()
            if any(term in line for term in ('Access is denied', 'RegOpenKeyExW', 'Result: 5',
                'memory dump file', 'Permission denied:', 'Sorry,', 'FATAL', 'Error:'))][-35:]
    finally:
        after = source_hashes()
        result['sources_conserved_during_build'] = before == after
        result['elapsed_s'] = round(time.monotonic() - started, 3)
        (NOTES / ('BUILD_RESULT_' + args.revision + '.json')).write_text(
            json.dumps(result, indent=2), encoding='utf-8')
        print(json.dumps({k: v for k, v in result.items()
                          if k not in ('sources_before', 'chapter_aux_files', 'target_label_pages')}, indent=2))
    return 0 if result['status'] == 'compiled' else 1


if __name__ == '__main__':
    raise SystemExit(main())
