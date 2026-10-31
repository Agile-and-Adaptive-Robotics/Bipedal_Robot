"""Preserve required bytes then remove live-code dependence on generated clones.
Does not delete any source folder or historical experiment data.
"""
import json,hashlib,shutil
from pathlib import Path
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[3];dest=HERE/'runs/source_inputs';records=HERE/'runs/archived_native_records';dest.mkdir(exist_ok=True);records.mkdir(exist_ok=True);rows=[];copied=[]
stems=['W2L_VirtualWalker_Guided','W2L_VirtualWalker_Guided30','W2L_VirtualWalker_GuidedT2','W2L_VirtualWalker_Newer_GUI1nA','W2L_VirtualWalker_Newer_SIContact','W2L_VirtualWalker_NewerV3Guided']
for stem in stems:
 folder=REPO/'Neuromechanical_Models'/f'{stem}_20261031';source=folder/f'{stem}_Standalone.asim';snapshot=dest/source.name;shutil.copy2(source,snapshot);sha=hashlib.sha256(source.read_bytes()).hexdigest();assert hashlib.sha256(snapshot.read_bytes()).hexdigest()==sha;rows.append(dict(original=str(source),snapshot=str(snapshot),sha256=sha))
 for file in folder.iterdir():
  if file.is_file() and (file.suffix.lower()=='.json' or file.name in ('EXPERIMENT_RECORD.md','README.md')):
   target=records/(stem+'__'+file.name);shutil.copy2(file,target);copied.append(dict(original=str(file),retained=str(target),sha256=hashlib.sha256(file.read_bytes()).hexdigest()))
archive=REPO/'tmp/w2l_generated_copies_before_cleanup_20261031.zip';archive_dest=HERE/'runs/archives'/archive.name;archive_dest.parent.mkdir(exist_ok=True);shutil.copy2(archive,archive_dest);assert hashlib.sha256(archive.read_bytes()).hexdigest()==hashlib.sha256(archive_dest.read_bytes()).hexdigest()
(dest/'MANIFEST.json').write_text(json.dumps(dict(machine='easteregg2',sources=rows,retained_records=copied,complete_original_archive=str(archive_dest),archive_sha256=hashlib.sha256(archive_dest.read_bytes()).hexdigest(),note='Historical immutable source inputs and evidence only; active W2L is canonical NEWER. No folder deletion performed.'),indent=2))
patched=[]
for file in HERE.glob('*.py'):
 if file.name in ('rehome_generated_sources.py','source_inputs.py'):continue
 before=file.read_text();text=before
 for stem in stems:
  text=text.replace("REPO/'Neuromechanical_Models/"+stem+"_20261031/"+stem+"_Standalone.asim'","historical_source('"+stem+"')")
 text=text.replace("REPO/'Neuromechanical_Models'/f'{stem}_20261031'/f'{stem}_Standalone.asim'","historical_source(stem)")
 if text!=before:
  text='from source_inputs import historical_source\n'+text;file.write_text(text);patched.append(str(file))
(dest/'REHOME_VERIFICATION.json').write_text(json.dumps(dict(sources=rows,patched_scripts=patched,all_snapshot_hashes_equal=True,full_archive_retained=str(archive_dest),deletion_performed=False),indent=2));print(json.dumps(dict(sources=len(rows),records=len(copied),patched=patched,archive=str(archive_dest)),indent=2))
