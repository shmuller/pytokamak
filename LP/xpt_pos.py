import numpy as np

import probe_xpr

from sig import PositionSignal

from sm_pyplot.tight_figure import get_fig

R0, Z0 = 1.645, -0.966

def get_xpt_pos(shn):
    IO = probe_xpr.IOMdsAUG(shn=shn, diag='GQH')
    return IO.load(('Rxpu', 'Zxpu', 't'))


def compare(shn):
    xpt_pos = get_xpt_pos(shn)
    Rxpu = PositionSignal(100*xpt_pos['Rxpu'], xpt_pos['t'])
    Zxpu = PositionSignal(100*xpt_pos['Zxpu'], xpt_pos['t'])

    
    XPR = probe_xpr.ProbeXPR(shn=shn)

    Rs = XPR.S['Rs']
    R = PositionSignal(100*(R0 - Rs.x), Rs.t).interp(Rxpu.t)

    z = np.zeros_like(R.x)
    z.fill(100*Z0)
    Z = PositionSignal(z, R.t).interp(Zxpu.t)
    
    dR = R.x - Rxpu.x
    dZ = Z.x - Zxpu.x
    D = PositionSignal(np.sqrt(dR**2 + dZ**2), R.t)

    im = D.x.argmin()
    print (D.t[im], D.x[im])

    fig = get_fig(shape=(3, 1), xlab='t (s)', ylab=('R (cm)', 'Z (cm)', 'Distance (cm)'))

    ax = fig.axes[0]
    R.plot(ax=ax)
    Rxpu.plot(ax=ax)

    ax = fig.axes[1]
    Z.plot(ax=ax)
    Zxpu.plot(ax=ax)

    ax = fig.axes[2]
    D.plot(ax=ax)

    return fig

