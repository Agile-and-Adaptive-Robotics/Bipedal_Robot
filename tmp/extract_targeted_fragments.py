from pathlib import Path
E=Path(r'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261023\Revision_E');V=E/'review/targeted_reconciliation'
s=(E/'chapters/60-sns.tex').read_text(encoding='utf-8');a=(E/'chapters/95-AppendixD.tex').read_text(encoding='utf-8')
start=s.index(r'\begin{figure}[p]');end=s.index(r'\end{figure}',start)+len(r'\end{figure}')
(V/'P04_figure61_block.tex').write_text(s[start:end]+'\n',encoding='utf-8')
paragraph=next(p for p in s.split('\n\n') if p.startswith('Years later, I split'))
(V/'P05_controller_history.tex').write_text(paragraph+'\n',encoding='utf-8')
start=a.index(r'\section{OpenSim-to-MuJoCo Conversion Repairs}');end=a.index(r'\section{Model Construction and Recording}',start)
(V/'P06_appendix_conversion_section.tex').write_text(a[start:end],encoding='utf-8')
start=s.index('My first conversion attempt');end=s.index('For the BPA-routed robot model',start)
(V/'P07_conversion_summary.tex').write_text(s[start:end],encoding='utf-8')
with (V/'OVERLEAF_TARGETED_PATCHES.md').open('a',encoding='utf-8') as f:
 f.write('''
## Bounded fragments for the already requested figure/history/conversion edits

- **P04_figure61_block.tex:** replace only the figure environment containing `\\label{fig:animatlab_model_architecture}`, from its own `\\begin{figure}` through `\\end{figure}`. Upload the supplied `connor_ben_shared_rg_original.png` asset to the existing figure directory before compiling. Preserve the label. This is not a source-file replacement.
- **P05_controller_history.tex:** replace only the paragraph beginning “I subsequently expanded the neural model into left- and right-side circuits” and ending with the Appendix simulation-verification reference. The replacement preserves Ben’s later RG/PF and contact-feedback additions and corrects the chronology. If that paragraph has since changed, compare it before applying.
- **P06_appendix_conversion_section.tex:** insert this section into the supplementary-verification appendix before “Model Construction and Recording,” only if that appendix exists and `sec:app_conversion_repairs` is absent. Do not replace the appendix file or create a second copy of its chapter.
- **P07_conversion_summary.tex:** replace only the paragraph beginning “MyoConverter supplies the initial OpenSim-to-MuJoCo translation” and ending with the 2 ms time-step sentence. This expands that paragraph into the concrete conversion summary. It stops before the existing BPA-interface paragraph. Apply P06 with it, or omit the new appendix-reference sentence until the target section exists.

These fragments are prepared for targeted review/application; none has been applied to Overleaf. P04/P05 follow Ben’s explicit figure and chronology correction. P06/P07 carry the source-supported repair details. Preserve and recheck each affected advisor thread after editing its span. This still does not constitute approval or completion of the remaining larger clarity proposals.
''')
print('Extracted P04–P07 as bounded fragments; no chapter uploads needed.')
