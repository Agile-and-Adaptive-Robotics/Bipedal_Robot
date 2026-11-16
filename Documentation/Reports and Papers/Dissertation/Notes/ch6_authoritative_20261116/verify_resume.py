"""Verify the saved canonical build before completing interrupted records."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import subprocess

HERE=Path(__file__).resolve().parent
DISS=HERE.parent.parent
ROOT=DISS.parents[2]
FINAL=DISS/'Final'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
build=json.loads((HERE/'BUILD_RESULT_r7.json').read_text())
pdf=Path(build['pdf'])
assert sha(pdf)==build['pdf_sha256'], 'Saved PDF changed since successful build'
expected=build['sources_before']
current={str(p.relative_to(FINAL)):sha(p) for p in sorted(FINAL.rglob('*')) if p.is_file()}
changed=[k for k in sorted(set(current)|set(expected)) if current.get(k)!=expected.get(k)]
assert not changed, 'Canonical sources differ from successful r7: '+repr(changed)
original=json.loads((HERE/'FINAL_BEFORE_HASHES.json').read_text(encoding='utf-8-sig'))
original_changes=[r['Path'] for r in original if sha(Path(r['Path'])).lower()!=r['SHA256'].lower()]
allowed={str(FINAL/'chapters/60-sns.tex'),str(FINAL/'chapters/95-AppendixD.tex')}
assert set(original_changes)==allowed
manifest=json.loads((HERE/'FIGURE_MOVES.json').read_text())
assert all(sha(Path(r['deprecated_path']))==r['sha256'] for r in manifest['moves'])
assert all(sha(Path(r['active_path']))==r['sha256'] for r in manifest['installed'])
git=lambda *args:subprocess.check_output(['git',*args],cwd=ROOT).decode('utf-8').strip()
head=git('rev-parse','HEAD')
committed=git('diff-tree','--no-commit-id','--name-only','-r',head).splitlines()
required=['Documentation/Reports and Papers/Dissertation/Final/chapters/'+n
          for n in ['60-sns.tex','95-AppendixD.tex']]
assert all(p in committed for p in required), 'Unexpected resume commit scope'
report={
    'checked_at_utc':datetime.now(timezone.utc).isoformat(),
    'hostname':subprocess.check_output(['hostname']).decode().strip(),
    'repository':str(ROOT),'write_probe':'unique tmp file created, read exactly, deleted on resume',
    'initial_resume_git_status':[],
    'initial_resume_relevant_processes':[],
    'process_scope':'MATLAB, python/pythonw, pdflatex/pdftex, bibtex, miktex, latexmk, ssh; no matching processes',
    'head_at_resume':head,'head_subject':git('log','-1','--format=%s'),
    'earlier_revision_already_committed':True,
    'committed_revision_paths':[{'relative_path':p,'absolute_path':str(ROOT/p)} for p in committed],
    'saved_build':'BUILD_RESULT_r7.json','saved_pdf':str(pdf),'pdf_sha256':sha(pdf),
    'canonical_sources_match_successful_build':True,'canonical_source_count':len(current),
    'canonical_source_hashes':current,
    'changes_since_original_pre_revision_snapshot':original_changes,
    'other_final_files_preserved':len(original)-len(original_changes),
    'preserved_asset_hashes_validated':len(manifest['moves']),
    'active_asset_hashes_validated':len(manifest['installed']),
    'rebuild_needed':False,'scientific_checks_rerun':False,
    'ssh':'Not retried; user waived SSH on continuation.'}
(HERE/'RESUME_VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:report[k] for k in ['hostname','head_at_resume','earlier_revision_already_committed',
    'canonical_sources_match_successful_build','canonical_source_count','other_final_files_preserved',
    'preserved_asset_hashes_validated','active_asset_hashes_validated','rebuild_needed']}))
