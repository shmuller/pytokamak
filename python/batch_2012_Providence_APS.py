from matplotlib.patches import Rectangle
from matplotlib import rc

rc('lines', linewidth=1.5)

import probe_xpr

R0 = 164.5

XPR  = probe_xpr.ProbeXPR(shn=27692)
XPR2 = probe_xpr.ProbeXPR(shn=27695)

def make_fig(keys):
    kw = dict(keys=keys, xkey='R')

    fig = XPR.res.plot(inout='in', figsize=(10,5), **kw)
    fig = XPR.res.plot(inout='out', fig=fig, **kw)
    fig = XPR2.res.plot(inout='in', fig=fig, **kw)
    return fig

def draw_rect(ax, x, color='b'):
    ym, yM = ax.get_ylim()
    rect = Rectangle((R0 - x - 0.5, ym), width=1, height=yM-ym, 
                     edgecolor=None, facecolor=color, alpha=0.2)
    ax.add_patch(rect)

def draw_all_rect(fig):
    for ax in fig.axes:
        draw_rect(ax, 17, 'b')
        draw_rect(ax, 27, 'b')
        draw_rect(ax, 18, 'r')
        draw_rect(ax, 24, 'r')

keys1 = ('n', 'Mach'), ('Vf', 'v')
keys2 = ('Te', 'mnv'), ('Vp', 'pe')

fname1 = 'XPRres_high_n_vs_low_n_part1'
fname2 = 'XPRres_high_n_vs_low_n_part2'

# part1
fig = make_fig(keys1)
fig.axes[1].set_ylim((-1.5, 1.5))
fig.axes[3].set_ylim((-50, 50))

draw_all_rect(fig)
fig.axes[0].legend(loc="upper right")

fig.savefig(fname1 + '.pdf')
fig.savefig(fname1 + '.png')

# part2
fig = make_fig(keys2)
draw_all_rect(fig)
fig.savefig(fname2 + '.pdf')
fig.savefig(fname2 + '.png')



