import os
from PIL import Image

figdir = r'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\30-torque'
names = [
    'Vas_Pam_20mm_cam_winner_1_r2.png',
    'Vas_Pam_20mm_cam_winner_2_r2.png',
    'Vas_Pam_20mm_cam_winner_3_r2.png',
    'route_check_run8_r2.png',
    'route_check_run8b_r2.png',
]
for n in names:
    p = os.path.join(figdir, n)
    im = Image.open(p)
    dpi = im.info.get('dpi', ('?', '?'))
    print('%-32s %sx%s px  dpi=%s  bytes=%d' % (n, im.size[0], im.size[1], dpi, os.path.getsize(p)))
    # print scale check at the QC's 6.2 in print width
    w_in = im.size[0] / float(dpi[0]) if dpi[0] != '?' else 0
    if w_in:
        print('   native %.2f in wide -> %.3fx at 6.2 in print; 13pt -> %.1f pt eff, 11pt -> %.1f pt eff'
              % (w_in, 6.2 / w_in, 13 * 6.2 / w_in, 11 * 6.2 / w_in))
