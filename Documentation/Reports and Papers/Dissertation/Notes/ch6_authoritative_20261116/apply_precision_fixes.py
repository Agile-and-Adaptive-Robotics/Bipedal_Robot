"""Apply bounded reviewer wording fixes to the two authorized manuscript files."""
from pathlib import Path
HERE=Path(__file__).resolve().parent
FINAL=HERE.parent.parent/'Final/chapters'
fixes={
'60-sns.tex':[
    (r'Minimum anatomical knee flexion was \ang{2.087} on the right and \ang{1.946} on the left.',
     r'Minimum anatomical knee flexion over $1\leq t<10$\,s was \ang{2.087} on the right and \ang{1.946} on the left.'),
    (r'outgoing \qty{0}{\milli\volt} transmission threshold',r'outgoing rest-relative \qty{0}{\milli\volt} transmission threshold'),
    ('source-setting Hunt benchmark',r'source-setting benchmark of the archived \texttt{hunt} implementation'),
    ('source-setting Hunt controller',r'source-setting \texttt{hunt} implementation'),
    ('1,033-cell Hunt-PT+RC controller','1,033-cell implementation labelled Hunt-PT+RC'),
    ('1,033-cell Hunt-PT+RC circuit','1,033-cell implementation labelled Hunt-PT+RC'),
],
'95-AppendixD.tex':[
    (r'Minimum anatomical knee flexion was \ang{2.087} on the right and \ang{1.946} on the left.',
     r'Minimum anatomical knee flexion over $1\leq t<10$\,s was \ang{2.087} on the right and \ang{1.946} on the left.'),
    (r'Its MN difference is below \qty{6e-8}{\milli\volt}.',r'Its MN difference is below \qty{6e-8}{\milli\volt}.'),
    (r'onset $\leq t < $ onset $+\qty{3.4}{\second}$',r'the interval from onset, inclusive, to \qty{3.4}{\second} after onset, exclusive'),
    ('on the original calibration clock without shifting individual channels to improve alignment.',
     'on the original calibration clock without shifting individual channels to improve alignment. Here $V$ is rest-relative membrane voltage.'),
    ('source-setting 1,025-cell Hunt controller',r'source-setting 1,025-cell implementation selected by \texttt{hunt}'),
    ('electrical parity and whole-process timing are checked separately','numerical agreement and scoped evaluator timing are checked separately'),
]}
for name,changes in fixes.items():
    p=FINAL/name
    text=p.read_text(encoding='utf-8')
    for old,new in changes:
        if old==new:
            continue
        if old not in text:
            if new in text:
                continue
            raise ValueError((name,old))
        if text.count(old)!=1:
            raise ValueError(('ambiguous',name,old,text.count(old)))
        text=text.replace(old,new,1)
    p.write_text(text,encoding='utf-8',newline='')
print('Applied bounded reviewer wording fixes.')
