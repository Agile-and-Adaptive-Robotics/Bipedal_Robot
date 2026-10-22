from pathlib import Path
import difflib
H=Path(__file__).resolve().parent
for n in ['60','70']:
    file={'60':'60-sns.tex','70':'70-overall.tex'}[n]
    a=(H/'overleaf_expected'/file).read_text(encoding='utf-8');b=(H/f'live{n}.txt').read_text(encoding='utf-8')
    print(file,'expected',len(a),'actual',len(b),'same',a==b)
    print(''.join(difflib.unified_diff(a.splitlines(True),b.splitlines(True)))[:5000])
