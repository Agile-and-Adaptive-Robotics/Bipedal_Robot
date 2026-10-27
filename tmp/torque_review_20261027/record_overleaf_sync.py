from pathlib import Path
import zipfile,difflib,json
R=Path(r'D:\GitHub\Bipedal_Robot')
P=R/'Documentation/Reports and Papers/Knee_Torque_Test'
A=R/'Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/94-AppendixC.tex'
O=R/'tmp/torque_review_20261027'
report=P/'CHATGPT_REVIEW_20261027.md'
t=report.read_text(encoding='utf-8')
t=t.replace('Nothing was committed, pushed, or uploaded to Overleaf.','Nothing was committed or pushed. The subsequent Overleaf upload is documented at the end of this report.')
t+='''

## Overleaf synchronization, 27 October 2026

Uploaded the revised sources and schematic PDFs to both user-selected projects and recompiled them:

- [Torque paper](https://www.overleaf.com/project/622e9536b8ea1e0a7f167cf5): Main.tex, frontiers_SupplementaryMaterial.tex, and both new schematic PDFs in Figures/. The manuscript compiles to 41 pages with **0 errors** and **0 undefined-reference warnings**. The 19 remaining warnings are font substitutions; typesetting messages also remain. Refreshed Main.aux and frontiers_SupplementaryMaterial.aux because older manually uploaded copies masked the new supplementary labels. The project's main document is left as Main.tex.
- [Dissertation](https://www.overleaf.com/project/644dc82bcd6a5481e3c3a8bb): chapters/15-background.tex, 30-torque.tex, 50-database.tex, 94-AppendixC.tex, and both new schematic PDFs in figs/30-torque/. Added the missing figs/30-torque/ fallback to Appendix C's graphicspath, both locally and on Overleaf, so its existing pulley and cam PDFs resolve there. The dissertation compiles with **0 errors** and no missing-figure errors. Its live project also contains Appendix E, which duplicates three database labels from Chapter 5: fig:sadb_explorer, tab:sadb_afferents, and tab:sadb_pathways. Overleaf reports these three conflicts plus their summary warning (4 warnings). Those other live appendix contents were left intact.

No optimization was run. No comments were resolved, and no existing tracked changes were accepted or rejected. The browser file helper lacked file-URL access, so file uploads used the standard Windows picker. Native mouse/keyboard control was stopped when Ben requested his screen back; subsequent checks used the Chrome extension.

Build screenshots are saved under tmp/torque_review_20261027/overleaf/. The changed-files ZIP and Appendix C session diff include the Overleaf path correction.
'''
report.write_text(t,encoding='utf-8')
old=(O/'94-AppendixC_before.tex').read_text(encoding='utf-8')
cur=A.read_text(encoding='utf-8')
(P/'Review_20261027/94-AppendixC.diff').write_text(''.join(difflib.unified_diff(old.splitlines(True),cur.splitlines(True),fromfile='94-AppendixC.tex (session start)',tofile='94-AppendixC.tex')),encoding='utf-8')
zp=P/'Review_20261027/Overleaf_changed_files_20261027.zip'
with zipfile.ZipFile(zp) as z: entries={n:z.read(n) for n in z.namelist()}
entries['dissertation/chapters/94-AppendixC.tex']=A.read_bytes()
entries['README.txt']+=b'\nSynchronization completed on 27 October 2026. Appendix C includes the figs/30-torque/ Overleaf fallback. See CHATGPT_REVIEW_20261027.md for live build results and remaining duplicate-label warnings.\n'
temp=zp.with_suffix('.new.zip')
with zipfile.ZipFile(temp,'w',zipfile.ZIP_DEFLATED) as z:
    for name,data in entries.items():z.writestr(name,data)
temp.replace(zp)
print('Recorded both Overleaf uploads and final build results; updated the source package.')
