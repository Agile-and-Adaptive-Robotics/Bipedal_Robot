"""Bounded read-only-installed-tree MiKTeX startup probe, no source changes."""
from pathlib import Path
import argparse, json, os, subprocess, time

root=Path(r'D:\GitHub\Bipedal_Robot')
probe=root/'tmp/ch6_authoritative_build_probe_20261116'
parser=argparse.ArgumentParser()
parser.add_argument('--preamble',action='store_true')
parser.add_argument('--refresh-fndb',action='store_true')
parser.add_argument('--priority-log')
args=parser.parse_args()
final=root/'Documentation/Reports and Papers/Dissertation/Final'
for fmt in (probe/'pdflatex.fmt',probe/'data/miktex/data/le/pdftex/pdflatex.fmt'):
    if fmt.exists(): os.utime(fmt,None)
env=os.environ.copy()
env.update(MIKTEX_USERCONFIG=str(probe/'config'), MIKTEX_USERDATA=str(probe/'data'),
           MIKTEX_USERINSTALL=r'D:\MiKTeX', MIKTEX_USERSTARTUPFILE=str(probe/'startup.ini'),
           MIKTEX_COMMONSTARTUPFILE=str(probe/'startup.ini'),
           TEXINPUTS=str(final)+r';.;D:\MiKTeX\tex\latex\base;D:\MiKTeX\tex\generic\config;D:\MiKTeX\tex//;',
           TFMFONTS=r'C:\Users\Ben Bolen\AppData\Local\MiKTeX\fonts\tfm//;D:\MiKTeX\fonts\tfm//;',
           TEXFONTMAPS=r'C:\Users\Ben Bolen\AppData\Local\MiKTeX\fonts\map\pdftex;D:\MiKTeX\fonts\map//;',
           T1FONTS=r'D:\MiKTeX\fonts\type1//;',
           ENCFONTS=r'D:\MiKTeX\fonts\enc//;',
           VFFONTS=r'D:\MiKTeX\fonts\vf//;',TEXFORMATS=str(probe))
if args.priority_log:
    from BUILD_texpaths import texinputs
    env['TEXINPUTS'],search_info=texinputs(final,args.priority_log)
    print(json.dumps(search_info))
inputfile=probe/'probe.tex'
if args.preamble:
    inputfile=probe/'preamble_probe.tex'
    preamble=(final/'main.tex').read_text(encoding='utf-8').split(r'\begin{document}',1)[0]
    inputfile.write_text(preamble+'\n'+r'\begin{document}Preamble and font probe. $\diameter\, x=1.$\end{document}',encoding='utf-8')
command=[r'D:\MiKTeX\miktex\bin\x64\pdflatex.exe','--miktex-disable-maintenance',
         '--miktex-disable-installer','-interaction=nonstopmode','-halt-on-error','-file-line-error',
         str(inputfile)]
if args.refresh_fndb:
    command=[r'D:\MiKTeX\miktex\bin\x64\initexmf.exe','--update-fndb=D:\MiKTeX','--disable-installer']
result={'command':command,'cwd':str(probe),'timeout_s':35}
start=time.monotonic()
logfile=probe/('preamble_probe.txt' if args.preamble else 'portable_probe3.txt')
if args.refresh_fndb: logfile=probe/'refresh_fndb.txt'
with logfile.open('w',encoding='utf-8') as out:
    process=subprocess.Popen(command,cwd=probe,env=env,stdout=out,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL)
    result['pid']=process.pid
    try:
        result['exit_code']=process.wait(timeout=35)
    except subprocess.TimeoutExpired:
        process.kill();process.wait(timeout=10)
        result.update(timed_out=True,owned_process_terminated=True)
result['elapsed_s']=round(time.monotonic()-start,3)
result['last_output']=logfile.read_text(errors='replace')[-5000:]
print(json.dumps(result,indent=2))
resultname='BUILD_FNDB_RESULT.json' if args.refresh_fndb else ('BUILD_PREAMBLE_RESULT.json' if args.preamble else 'BUILD_PROBE_RESULT.json')
(Path(__file__).parent/resultname).write_text(json.dumps(result,indent=2),encoding='utf-8')
