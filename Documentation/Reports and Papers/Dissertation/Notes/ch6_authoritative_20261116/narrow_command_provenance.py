"""Distinguish declared adapter commands from missing original input telemetry."""
from pathlib import Path
HERE=Path(__file__).resolve().parent
CH=HERE.parent.parent/'Final/chapters'
fixes={
'60-sns.tex':[
('Stage 1 disables sensory return and sustained flexor/extensor commands; stage 2 enables live sensory return and sets both normalized commands to 0.15 bilaterally. The period change therefore reflects different command and sensory conditions, rather than a feedback-only intervention.',
 'The retained adapter defines stage 1 with sensory return and sustained flexor/extensor commands disabled, and stage 2 with live sensory return and both normalized commands set to 0.15 bilaterally. Complete applied-command histories and the original launch adapter were not recovered, so these definitions cannot establish sample-by-sample command identity in the original runs. The period change therefore cannot be interpreted as a controlled feedback-only effect.')],
'95-AppendixD.tex':[
('Stage 1 disables sensory return and sets the sustained normalized flexor and extensor commands to zero. Stage 2 enables live sensory return and sets each command to 0.15 on both sides.',
 'The retained adapter defines stage 1 with sensory return disabled and sustained normalized flexor and extensor commands set to zero. It defines stage 2 with live sensory return enabled and each command set to 0.15 on both sides, consistent with the fitted command manifest.'),
('Thus these air stages are not a matched feedback-only comparison.',
 'The original launch adapter and complete consumed-command traces were not recovered; these are declared adapter settings rather than independently verified input samples from that launch. Thus these air stages cannot be treated as a matched feedback-only comparison.')]
}
fixes['95-AppendixD.tex'].extend([
('Both begin with the same initialization seed,', 'That adapter also defines the same initialization seed for both stages,'),
('The postural, balance, and functional-command channels are clamped off, direct tonic RG current is zero, and the legacy shared-drive input is disabled in this implementation.',
 'It clamps the postural, balance, and functional-command channels off, sets direct tonic RG current to zero, and disables the legacy shared-drive input.')])
for name, pairs in fixes.items():
    p=CH/name;text=p.read_text()
    for old,new in pairs:
        if old not in text:
            if new in text:continue
            raise ValueError((name,old))
        assert text.count(old)==1
        text=text.replace(old,new)
    p.write_text(text,encoding='utf-8',newline='')
print('Narrowed command provenance in both canonical files.')
