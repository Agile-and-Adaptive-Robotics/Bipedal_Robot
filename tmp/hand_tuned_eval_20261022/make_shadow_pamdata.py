import io, os
base = r"D:\GitHub\Bipedal_Robot\tmp\hand_tuned_eval_20261022"
src = io.open(os.path.join(base, 'PamData_20261022.m'), newline='').read()
nl = '\r\n' if '\r\n' in src else '\n'
assert src.count('classdef PamData_20261022' + nl) == 1
assert src.count('function PD = PamData_20261022(muscle') == 1
out = src.replace('classdef PamData_20261022' + nl, 'classdef PamData' + nl, 1)
out = out.replace('function PD = PamData_20261022(muscle', 'function PD = PamData(muscle', 1)
# replace the staging header with the shadow header
old_hdr_start = '% STAGED COPY 2026-10-22 (tmp/capability_table_20261022 on eb475ws4).'
i0 = out.index(old_hdr_start)
i1 = out.index('classdef PamData') - 1  # up to the newline before classdef... find line end
# find the end of the header block (last header line ends right before classdef line)
hdr = ('% SHADOW COPY 2026-10-22 (tmp/capability_table_20261022 on eb475ws4).' + nl +
       '% Purpose: let Connor''s section scripts run VERBATIM in this driver; the repo' + nl +
       '% PamData.m does not parse (9hi corruption, line 271). This shadow is the' + nl +
       '% staged PamData_20261022.m with the class name reverted to PamData - the' + nl +
       '% same single-token fix, no value changes. It lives ONLY in this tmp folder' + nl +
       '% and shadows the repo class solely via path order in this run.' + nl)
out = out[:i0] + hdr + out[i1:]
io.open(os.path.join(base, 'PamData.m'), 'w', newline='').write(out)
print('shadow written OK, %d bytes' % len(out))
