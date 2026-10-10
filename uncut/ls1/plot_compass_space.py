"""Render LS1's exact coefficient surface and its noncollapsed zero seam."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
fig = plt.figure(figsize=(12, 6.5), facecolor='white')
fig.suptitle('Compass cut space', fontsize=22, fontweight='bold', y=.97)
fig.text(.5, .91, 'Two symmetric cut components (u, v) and a signed turn component (r)',
         ha='center', fontsize=12, color='#44546A')
ax = fig.add_subplot(121, projection='3d')
angles, depths = np.meshgrid(np.linspace(0, 2*np.pi, 80), np.linspace(-1.35, 1.35, 43))
u = np.sqrt(1+depths**2)*np.cos(angles)
v = np.sqrt(1+depths**2)*np.sin(angles)
ax.plot_surface(u, v, depths, color='#4C82B8', alpha=.16, edgecolor='none')
ax.plot_wireframe(u, v, depths, rstride=7, cstride=10, color='#4C82B8', linewidth=.55, alpha=.5)
angle = np.linspace(0, 2*np.pi, 400)
ax.plot(np.cos(angle), np.sin(angle), np.zeros_like(angle), color='#CE6419', linewidth=3.5)
ax.set_xlabel('u', labelpad=4)
ax.set_ylabel('v', labelpad=4)
ax.set_zlabel('r', labelpad=4)
ax.set_title(r'$u^2+v^2-r^2=1$', pad=3, fontsize=15)
ax.view_init(elev=19, azim=-56)
ax.set_box_aspect((1, 1, .92))
ax.set_xticks([-1, 0, 1]); ax.set_yticks([-1, 0, 1]); ax.set_zticks([-1, 0, 1])
ax.grid(False)
ax2 = fig.add_subplot(122)
ax2.plot(np.cos(angle), np.sin(angle), color='#CE6419', linewidth=3)
ax2.axhline(0, color='#D7DEE7', zorder=0, linewidth=1)
ax2.axvline(0, color='#D7DEE7', zorder=0, linewidth=1)
for x, y, label, tx, ty in ((1, 0, r'$K$', 1.16, 0), (0, 1, r'$RK$', 0, 1.17),
                           (-1, 0, r'$-K$', -1.2, 0), (0, -1, r'$-RK$', 0, -1.17)):
    ax2.scatter([x], [y], s=55, color='#293F62', zorder=4)
    ax2.text(tx, ty, label, ha='center', va='center', fontsize=14, color='#293F62')
ax2.text(0, .18, 'r = 0', ha='center', fontsize=18, color='#CE6419', fontweight='bold')
ax2.text(0, -.15, 'Distinct cuts\nSame return B = I', ha='center', va='center', fontsize=13, color='#44546A')
ax2.set_title('The zero seam is a circle', fontsize=15, pad=12)
ax2.set_aspect('equal'); ax2.set_xlim(-1.45, 1.45); ax2.set_ylim(-1.45, 1.45)
ax2.axis('off')
fig.text(.5, .06, 'The full cut retains orientation at zero. Keeping only B would collapse this circle to one point.',
         ha='center', fontsize=11, color='#44546A')
fig.subplots_adjust(left=.035, right=.965, bottom=.13, top=.82, wspace=.12)
fig.savefig(HERE/'LS1_compass_space.png', dpi=170, facecolor='white')
