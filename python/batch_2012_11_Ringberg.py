from matplotlib.patches import Rectangle
from matplotlib import rc

rc('lines', linewidth=1.5)

from LP import probe_xpr

R0 = 164.5

XPR  = probe_xpr.ProbeXPR(shn=27691)
XPR2 = probe_xpr.ProbeXPR(shn=28795)
XPR3 = probe_xpr.ProbeXPR(shn=28799)
XPR4 = probe_xpr.ProbeXPR(shn=28797)


keys = ('n', 'Mach'), ('Vf', 'pe')

kw = dict(keys=keys, xkey='R', inout='in')

fname = 'XPRres_Leena'

fig = XPR.res.plot(figsize=(10,5), **kw)
for ax in fig.axes:
    ax.set_xlim((130, 165))
fig.axes[0].set_ylim((0, 2.5e19))
fig.axes[1].set_ylim(-1.5, 2)

fig.savefig(fname + '1.pdf')
fig.savefig(fname + '1.png')

fig = XPR2.res.plot(fig=fig, **kw)
fig.savefig(fname + '2.pdf')
fig.savefig(fname + '2.png')

fig = XPR3.res.plot(fig=fig, **kw)
fig.savefig(fname + '3.pdf')
fig.savefig(fname + '3.png')

fig = XPR4.res.plot(fig=fig, **kw)
fig.savefig(fname + '4.pdf')
fig.savefig(fname + '4.png')




