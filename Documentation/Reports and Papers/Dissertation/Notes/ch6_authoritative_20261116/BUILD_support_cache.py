"""Copy byte-identical installed support files to a writable runtime cache.

Only directories observed in the reference TeX log are eligible. Conflicting
basenames are excluded, retaining normal directory resolution for those files.
No manuscript, installed support file, or scientific data is changed.
"""
from pathlib import Path
import hashlib, json, re, shutil
from BUILD_texpaths import texinputs

def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def prepare(final,reference_log,runtime):
    search,info=texinputs(final,reference_log,False)
    dirs=[Path(p) for p in search.split(';') if p and p!='.' and Path(p)!=final]
    texcache=runtime/'package_cache';texcache.mkdir()
    fontcache=runtime/'font_cache';fontcache.mkdir()
    seen={};conflicts=set();records=[]
    def copy(file,out):
        key=(out,file.name.lower())
        sha=digest(file)
        if key in conflicts: return
        if key in seen:
            if seen[key]['sha256']==sha: return
            (out/file.name).unlink(missing_ok=True)
            conflicts.add(key)
            return
        target=out/file.name
        shutil.copy2(file,target)
        if digest(target)!=sha:
            raise RuntimeError('Installed support cache differs from its source: '+str(file))
        record={'source':str(file),'cache':str(target),'sha256':sha}
        seen[key]=record
    for directory in dirs:
        for file in directory.iterdir():
            if file.is_file() and file.suffix.lower() in ('.tex','.ltx','.sty','.cls','.clo','.fd','.def','.cfg','.trsl','.mkii','.dcl'):
                copy(file,texcache)
    prefixes=re.compile(r'^(?:ppl|zpl|cm|lasy|ectt|sftt|math[abx]|ms[ab]|line|lcircle|upl|fpl|pazo|eur|eus)',re.I)
    for root in [Path(r'D:\MiKTeX\fonts\tfm'),Path(r'D:\MiKTeX\fonts\vf'),
                 Path(r'D:\MiKTeX\fonts\type1'),Path(r'C:\Users\Ben Bolen\AppData\Local\MiKTeX\fonts\tfm')]:
        for file in root.rglob('*'):
            if file.suffix.lower() in ('.tfm','.vf','.pfb') and prefixes.match(file.name) and file.is_file():
                copy(file,fontcache)
    for file in (Path(r'C:\Users\Ben Bolen\AppData\Local\MiKTeX\fonts\map\pdftex\pdftex.map'),
                 Path(r'D:\MiKTeX\fonts\enc\dvips\base\8r.enc')):
        copy(file,fontcache)
    copy(Path(r'D:\MiKTeX\bibtex\bst\natbib\plainnat.bst'),texcache)
    records=[record for key,record in seen.items() if key not in conflicts]
    manifest={'package_directory':str(texcache),'font_directory':str(fontcache),
        'copied_file_count':len(records),'files':records,
        'conflicting_names_excluded':sorted(name for _,name in conflicts),
        'scope':'Byte-identical copies of installed support/font files; no source or installed-file changes.'}
    manifestpath=runtime/'support_cache_manifest.json'
    manifestpath.write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    return texcache,fontcache,{'manifest':str(manifestpath),'manifest_sha256':digest(manifestpath),
        'copied_file_count':len(records),'conflicting_names_excluded':manifest['conflicting_names_excluded']}
