"""Synchronize other Final chapters with Ben's latest downloaded source.

Keep the revised Chapter 6 and Appendix D, and the verified drive clarification.
No bibliography, figure, model, or result writes. No ZIP deliverable is created.
"""
from pathlib import Path, PurePosixPath
import argparse
import difflib
import hashlib
import json
import re
import shutil
import zipfile

ROOT=Path(__file__).resolve().parents[5]
DISS=ROOT/'Documentation/Reports and Papers/Dissertation'
FINAL=DISS/'Final'
NOTES=DISS/'Notes/ch6_revision_20261116'
SOURCE=DISS/'Bolen_Dissertation.zip'
BACKUP=DISS/'Archive/final_source_sync_20261116'
KEEP={'chapters/60-sns.tex','chapters/95-AppendixD.tex'}
OLD=('Human-model stepping persists after removal of externally supplied neural drive, '
     'in place under the support rig; the robot model progresses under its harness with reduced drive')
NEW=('In the human-model experiment, stepping persists after the supplementary left/right '
     'hand-drive channels are removed, in place under the support rig; the shared locomotor-drive '
     'input remains active. In the separate robot-model experiment, the robot model progresses '
     'under its harness with reduced drive')

def sha(data):return hashlib.sha256(data).hexdigest()
def normalized(data):return data.replace(b'\r\n',b'\n')

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--apply',action='store_true')
    args=parser.parse_args()
    kept={name:sha((FINAL/name).read_bytes()) for name in KEEP}
    jobs=[]
    patches=[]
    with zipfile.ZipFile(SOURCE) as archive:
        for name in sorted(archive.namelist()):
            parts=PurePosixPath(name).parts
            if len(parts)!=2 or parts[0]!='chapters' or not name.endswith('.tex') or name in KEEP:continue
            target=FINAL/name
            if not target.resolve().is_relative_to(FINAL.resolve()):raise RuntimeError('Invalid chapter path')
            before=target.read_bytes()
            supplied=archive.read(name)
            after=supplied
            exception=None
            if name=='chapters/70-overall.tex':
                if supplied.count(OLD.encode())!=1:raise RuntimeError('Reconcile Chapter 7 drive wording')
                if b'the shared locomotor-drive input remains active' not in before:
                    raise RuntimeError('Expected verified current clarification not present')
                after=supplied.replace(OLD.encode(),NEW.encode())
                exception='Retain verified supplementary/shared descending-drive clarification.'
            backup=BACKUP/name
            if before!=after and backup.exists() and backup.read_bytes()!=before:
                raise RuntimeError(f'Preserved prior version already exists with different bytes: {backup}')
            jobs.append({'chapter':name,'before':before,'supplied':supplied,'after':after,'backup':backup,
                         'exception':exception})
            patches.extend(difflib.unified_diff(
                normalized(before).decode('utf-8').splitlines(True),
                normalized(after).decode('utf-8').splitlines(True),
                fromfile='before/'+name,tofile='Final/'+name))
    if args.apply:
        for job in jobs:
            target=FINAL/job['chapter']
            if job['before']!=job['after']:
                job['backup'].parent.mkdir(parents=True,exist_ok=True)
                if not job['backup'].exists():shutil.copy2(target,job['backup'])
                target.write_bytes(job['after'])
            if target.read_bytes()!=job['after']:raise RuntimeError('Imported chapter differs')
        for name,value in kept.items():
            if sha((FINAL/name).read_bytes())!=value:raise RuntimeError(f'Protected rewrite changed: {name}')
    records=[{'chapter':j['chapter'],'changed':j['before']!=j['after'],
              'text_changed':normalized(j['before'])!=normalized(j['after']),
              'before_sha256':sha(j['before']),'zip_sha256':sha(j['supplied']),
              'final_sha256':sha(j['after']),
              'zip_byte_identical':j['after']==j['supplied'],
              'exception':j['exception'],
              'backup':str(j['backup'].relative_to(ROOT)) if j['before']!=j['after'] else None}
             for j in jobs]
    report={'mode':'applied' if args.apply else 'plan','source':str(SOURCE.relative_to(ROOT)),
            'source_sha256':sha(SOURCE.read_bytes()),'chapters_synchronized':len(jobs),
            'preserved_rewritten_chapters':kept,'chapters':records,
            'main_tex':'Already textually identical to downloaded main.tex; no change required.',
            'bibliographies_and_figures':'Unchanged; independent dependency check found no additions required.'}
    (NOTES/('FINAL_SOURCE_SYNC.json' if args.apply else 'FINAL_SOURCE_SYNC_PLAN.json')).write_text(json.dumps(report,indent=2),encoding='utf-8')
    if args.apply:(NOTES/'OTHER_CHAPTERS_SYNC.diff').write_text(''.join(patches),encoding='utf-8')
    print(json.dumps({'mode':report['mode'],'chapters_synchronized':len(jobs),
                      'text_changes':[j['chapter'] for j in records if j['text_changed']],
                      'rewrites_preserved':sorted(KEEP),
                      'drive_exception':'chapters/70-overall.tex','backup':str(BACKUP)},indent=2))

if __name__=='__main__':main()
