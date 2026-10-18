import re, os, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
os.chdir(r"D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Knee_Torque_Test\ProofFinal")

main = open('Main.tex', encoding='utf-8').read()
keys = set()
for m in re.finditer(r'\\(citep|citet|cite)\s*(\[[^\]]*\])?\{([^}]*)\}', main):
    for k in m.group(3).split(','):
        if k.strip():
            keys.add(k.strip())
bib = open('BolenFrontiers22.bib', encoding='utf-8').read()
bibkeys = set(k.strip() for k in re.findall(r'@\w+\{([^,]+),', bib))
missing = sorted(k for k in keys if k not in bibkeys)
print('cite keys used:', len(keys), '| missing from BolenFrontiers22.bib:', missing if missing else 'NONE')

figs = ['testJigs.pdf', 'KneeICR.eps', 'xiFrameTransform.pdf', 'bktFrame.pdf',
        'FlxPin_group.eps', 'ExtPin_group.eps', 'Ext10mm_52cm.eps', 'FlxGrp.eps',
        'ExtensorRoute_redesigned.pdf']
for f in figs:
    print(('OK       ' if os.path.isfile('Figures/' + f) else 'MISSING  ') + f)
sup = open('frontiers_SupplementaryMaterial.tex', encoding='utf-8').read()
print('supplementary passThru label present:', 'fig:passThru' in sup)
