# gen_bpaWrapping_alpha_20261026.py -- P2-5 dissertation figure: the
# wrapping-geometry panel of the old composite fig:xiProjection as its own
# figure, with the bend angle alpha MARKED ON THE CIRCLE (Ben, 2026-10-26
# 13:38: "Should be two separate figures. alpha can be marked on the
# circle."). Drawing adapted verbatim from the established generator
# tmp/revision_c_user_assets.py (same palette/layout); the change is the
# alpha arc drawn at the joint circle between the vertical radius reference
# and the wrap-surface end radius, replacing the floating alpha annotation.
# Output: tmp/paper_bridge/figscripts/bpa_wrapping_loss_20261026.pdf (+png).
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Arc

blue = '#2A5DB0'; ink = '#20262E'; orange = '#D97E00'; grey = '#E9EEF2'

fig, ax = plt.subplots(figsize=(7.0, 2.4))
ax.set_aspect('equal'); ax.axis('off')
ax.set_xlim(-2.6, 2.9); ax.set_ylim(-.25, 1.65)

# joint circle
ax.add_patch(Circle((0, 0), .91, color=grey, ec='#80909A', lw=1.1))
# BPA following the joint surface (thick blue arc over 30..150 deg)
angles = np.deg2rad(np.linspace(30, 150, 120))
ax.plot(np.cos(angles), np.sin(angles), color=blue, lw=9, solid_capstyle='round')
# tangent BPA legs leaving the surface
for a in [30, 150]:
    th = np.deg2rad(a)
    pt = np.array([np.cos(th), np.sin(th)])
    direction = np.array([-np.sin(th), np.cos(th)]) * (-1 if a == 150 else 1)
    end = pt - direction * .8
    ax.plot([pt[0], end[0]], [pt[1], end[1]], color=blue, lw=9, solid_capstyle='round')
    ax.plot([0, pt[0]], [0, pt[1]], ls='--', color=ink, lw=.8)
# R arrow (radius reference, vertical)
ax.annotate('', (0, 1), (0, 0), arrowprops={'arrowstyle': '<->', 'lw': 1, 'color': orange})
ax.text(.08, .65, r'$R$', color=orange)
# wrap-surface angle sweep label (the arc between the radius lines)
ax.add_patch(Arc((0, 0), .75, .75, theta1=30, theta2=150, color=orange, lw=1.2))
ax.text(-.75, 1.28, 'BPA following joint surface', fontsize=10)
# ALPHA MARKED ON THE CIRCLE: arc between the vertical radius (R) and the
# 30-deg end radius, with the label on the arc.
ax.add_patch(Arc((0, 0), .34, .34, theta1=30, theta2=90, color=ink, lw=1.4))
ax.text(.24, .44, r'$\alpha$', fontsize=13, color=ink)
ax.text(1.05, .95, r'$\alpha=|\theta_{\mathrm{wrap}}-\theta_{k}|$', fontsize=11)
ax.text(1.1, .68, r'$\Delta l=\chi_3 R\alpha\zeta^2$', fontsize=11)
ax.text(1.1, .40, 'Fitted usable-length loss', fontsize=9)
ax.text(-.28, -.12, 'Joint', fontsize=9)
fig.subplots_adjust(left=.01, right=.99, top=.99, bottom=.01)
fig.savefig('bpa_wrapping_loss_20261026.pdf')
fig.savefig('bpa_wrapping_loss_20261026.png', dpi=200)
print('WROTE bpa_wrapping_loss_20261026.pdf/.png')
