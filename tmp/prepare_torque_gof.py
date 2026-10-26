from pathlib import Path
import json
w=Path(r'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Knee_Torque_Test\review_audit_20261026')
p=(w/'paper_Main_live.tex').read_text(encoding='utf-8')
d=(w/'dissertation_30-torque_live.tex').read_text(encoding='utf-8')
paragraph=next(x for x in p.split('\n\n') if x.startswith('The correction terms enter the prediction'))
oldp=paragraph[:paragraph.index('The solver was')]
definitions=r'''The correction terms were identified by global optimization against measured isometric torque. Three goodness-of-fit (GoF) objectives were evaluated over the tested knee angles: root mean square error (RMSE), fraction of variance unexplained (FVU), and maximum absolute residual. RMSE summarizes the overall error magnitude. FVU divides the sum of squared residuals by the sum of squared deviations of the measured torque from its mean, expressing the error relative to the torque variation in that test. The maximum absolute residual records the largest pointwise error. RMSE and maximum residual have torque units, whereas FVU is dimensionless; smaller values indicate closer agreement for all three.

Together, these measures assess agreement in torque magnitude and variation across the range of motion. Offsets, amplitude differences, and shifts along the knee-angle axis can affect several metrics at once, so the metrics are interpreted jointly rather than assigning each one to a single type of disagreement. Across tests with different torque ranges, the absolute and variance-normalized errors also weight the tests differently. Multiobjective optimization retains these comparisons alongside the largest local error, making the trade-offs visible instead of allowing a favorable average score to conceal a poor fit elsewhere. The resulting Pareto front contains alternative compromises among the objectives, from which the adopted terms are selected using held-out validation.

'''
old_start=d.index('The goodness-of-fit metrics used to judge')
old_end=d.index('\\subsection{Optimization Algorithms',old_start)
oldd=d[old_start:old_end]
newd=definitions.replace('The correction terms were identified by global optimization against measured isometric torque. ','',1)
newd+=r'''This joint interpretation is consistent with describing curve disagreement in terms of offset, amplitude, and pattern \citep{taylor_summarizing_2001, iosa_assessment_2014}, but none of the three GoF measures identifies one component exclusively. In particular, FVU is based on the full squared residual, including any mean offset; it is the complement of the Nash--Sutcliffe model efficiency and is below one when the prediction outperforms a constant prediction equal to the measured mean \citep{nash_river_1970}. FVU provides the variance-normalized comparison across tests, while RMSE and maximum residual retain the absolute error information needed to judge torque predictions.

'''
patches=[{'document':'paper','old':oldp,'new':definitions,'reason':'Explain GoF jointly and the multiobjective trade-offs; retain the existing solver and selection procedure.'},{'document':'dissertation','old':oldd,'new':newd,'reason':'Replace the one-metric-per-error interpretation with a joint assessment and explicit multiobjective rationale.'}]
for c in patches:
    source=p if c['document']=='paper' else d
    assert source.count(c['old'])==1
(w/'gof_prose_patches.json').write_text(json.dumps(patches,indent=2),encoding='utf-8')
print('Prepared two bounded GoF prose replacements')
