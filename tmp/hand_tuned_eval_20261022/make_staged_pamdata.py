import io, os
base = r"D:\GitHub\Bipedal_Robot\tmp\hand_tuned_eval_20261022"
src = io.open(os.path.join(base, 'connor', 'PamData.m'), newline='').read()
nl = '\r\n' if '\r\n' in src else '\n'
assert src.count(';9hi ') == 1, 'unexpected 9hi count: %d' % src.count(';9hi ')
assert src.count('classdef PamData' + nl) == 1
assert src.count('function PD = PamData(muscle') == 1
out = src.replace(';9hi ', ';', 1)
out = out.replace('classdef PamData' + nl, 'classdef PamData_20261022' + nl, 1)
out = out.replace('function PD = PamData(muscle', 'function PD = PamData_20261022(muscle', 1)
hdr = ("% STAGED COPY 2026-10-22 (tmp/capability_table_20261022 on eb475ws4)." + nl +
       "% Source: repo Code/Matlab/Previous Optimization Code/PamData.m, VERBATIM except:" + nl +
       "%   (1) class renamed PamData -> PamData_20261022;" + nl +
       "%   (2) line 271 'size = delta/(1 - 1/def);9hi ' -> 'size = delta/(1 - 1/def);'" + nl +
       "%       (the '9hi' token is a syntax corruption; Ben's original Size.m has the" + nl +
       "%       clean line; the repo class FAILS TO PARSE as-is - verified 2026-10-22)." + nl +
       "% No design values changed. The repo file was NOT modified." + nl)
out = out.replace('classdef PamData_20261022', hdr + 'classdef PamData_20261022', 1)
io.open(os.path.join(base, 'PamData_20261022.m'), 'w', newline='').write(out)
print('written OK, %d bytes' % len(out))
