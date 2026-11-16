"""Test writable MiKTeX core-cache roots using unchanged installed config bytes."""
from pathlib import Path
import hashlib, json, os, shutil, subprocess, time
root=Path(r'D:\GitHub\Bipedal_Robot')
notes=Path(__file__).resolve().parent
probe=root/'tmp/ch6_authoritative_private_install_probe_20261116'
cache=root/'tmp/ch6_authoritative_20261116/build_final_r4/runtime'
final=root/'Documentation/Reports and Papers/Dissertation/Final'
probe.mkdir(exist_ok=True)
install=probe/'install'
shutil.copytree(Path(r'D:\MiKTeX\miktex\config'),install/'miktex/config',dirs_exist_ok=True)
(probe/'config/miktex/config').mkdir(parents=True,exist_ok=True)
(probe/'data/miktex/data/le/pdftex').mkdir(parents=True,exist_ok=True)
startup=probe/'startup.ini'
extra_fonts=probe/'extra_fonts';extra_fonts.mkdir(exist_ok=True)
extra_tex=probe/'extra_tex';extra_tex.mkdir(exist_ok=True)
shutil.copy2(Path(r'D:\MiKTeX\tex\generic\mathabx\mathabx.dcl'),extra_tex/'mathabx.dcl')
font_records=[]
for fontroot in (Path(r'D:\MiKTeX\fonts\tfm'),Path(r'D:\MiKTeX\fonts\type1')):
    for prefix in ('lasy','ectt','sftt'):
        for fontfile in fontroot.rglob(prefix+'*'):
            if fontfile.is_file() and fontfile.suffix in ('.tfm','.pfb'):
                shutil.copy2(fontfile,extra_fonts/fontfile.name)
                font_records.append({'source':str(fontfile),'sha256':hashlib.sha256(fontfile.read_bytes()).hexdigest()})
startup.write_text('[Auto]\nConfig=Portable\n[Paths]\n'+''.join(key+'='+str(value).replace('\\','/')+'\n' for key,value in {
    'UserInstall':install,'CommonInstall':install,'UserConfig':probe/'config','CommonConfig':probe/'config',
    'UserData':probe/'data','CommonData':probe/'data'}.items()),encoding='ascii')
(probe/'config/miktex/config/miktex.ini').write_text('[TeXandFriends]\nRenewFormatsOnUpdate=false\n',encoding='ascii')
fmt=probe/'data/miktex/data/le/pdftex/pdflatex.fmt'
shutil.copy2(Path(r'C:\Users\Ben Bolen\AppData\Local\MiKTeX\miktex\data\le\pdftex\pdflatex.fmt'),fmt)
os.utime(fmt,None)
env=os.environ.copy()
env.update(MIKTEX_USERCONFIG=str(probe/'config'),MIKTEX_COMMONCONFIG=str(probe/'config'),
    MIKTEX_USERDATA=str(probe/'data'),MIKTEX_COMMONDATA=str(probe/'data'),
    MIKTEX_USERINSTALL=str(install),MIKTEX_COMMONINSTALL=str(install),
    MIKTEX_USERSTARTUPFILE=str(startup),MIKTEX_COMMONSTARTUPFILE=str(startup),
    TEXINPUTS=str(final)+';.;'+str(extra_tex)+';'+str(cache/'package_cache'),TEXFORMATS=str(fmt.parent),
    TFMFONTS=str(extra_fonts)+';'+str(cache/'font_cache'),VFFONTS=str(cache/'font_cache'),T1FONTS=str(extra_fonts)+';'+str(cache/'font_cache'),
    TEXFONTMAPS=str(cache/'font_cache'),ENCFONTS=str(cache/'font_cache'))
inputfile=probe/'probe.tex'
preamble=(final/'main.tex').read_text(encoding='utf-8').split(r'\begin{document}',1)[0]
inputfile.write_text(preamble+'\n'+r'\begin{document}Preamble and font probe. $\diameter\,x=1.$\end{document}',encoding='utf-8')
cmd=[r'D:\MiKTeX\miktex\bin\x64\pdflatex.exe','--miktex-disable-maintenance','--miktex-disable-installer',
    '-interaction=nonstopmode','-halt-on-error','-file-line-error',str(inputfile)]
result={'command':cmd,'cwd':str(probe),'timeout_s':45,'installed_files_changed':False,'extra_fonts':font_records}
tick=time.monotonic()
with (probe/'probe.txt').open('w',encoding='utf-8') as stream:
    process=subprocess.Popen(cmd,cwd=probe,env=env,stdout=stream,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL)
    result['pid']=process.pid
    try: result['exit_code']=process.wait(timeout=45)
    except subprocess.TimeoutExpired:
        process.kill();process.wait(timeout=10)
        result.update(timed_out=True,owned_process_terminated=True)
result['elapsed_s']=round(time.monotonic()-tick,3)
result['last_output']=(probe/'probe.txt').read_text(errors='replace')[-5000:]
(notes/'BUILD_PRIVATE_INSTALL_PROBE_RESULT.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
