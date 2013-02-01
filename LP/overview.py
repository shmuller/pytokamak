import numpy as np
import numpy.ma as ma

import probe_xpr
from sig import memoized_property

from sm_pyplot.tight_figure import get_fig, get_axes


class AUGOverview:
    def __init__(self, shn):
        self.shn = shn

        self.chn = dict(
            DCN = ('H-1', 'H-2', 'H-3', 'H-4', 'H-5'),
            ECS = ('PECRH',),
            NIS = ('PNI',),
            GQI = ('Wmhd',),
            MAC = ('Ipolsola', 'Ipolsoli'))

        self.XPR = probe_xpr.ProbeXPR(shn=shn)

        self.def_plots = ('power', 'density', 'XPR_I', 'XPR_R', 'Ipolsol')

    def load(self, diag):
        IO = probe_xpr.IOMdsAUG(shn=self.shn, diag=diag)
        chn = self.chn[diag] + ('t',)
        try:
            return IO.load(chn)
        except probe_xpr.probe.TdiError:
            return {c: 0. for c in chn}

    @memoized_property
    def S(self):
        return {k: self.load(k) for k in self.chn.keys()}

    def plot_power(self, ax):
        ax = get_axes(ax)
        ax.set_ylabel('Power (MW)')
        S_ECRH = self.S['ECS']
        S_NI   = self.S['NIS']
        S_WMHD = self.S['GQI']

        ax.plot(S_ECRH['t'], 1e-6*S_ECRH['PECRH'], label="ECRH")
        ax.plot(S_NI['t'], 1e-6*S_NI['PNI'], label="NBI")
        ax.plot(S_WMHD['t'], 1e-5*S_WMHD['Wmhd'], label="WMHD (x10)")
        ax.legend()
        return ax

    def plot_density(self, ax):
        ax = get_axes(ax)
        ax.set_ylabel('n (10$^{\mathdefault{19}}$ m$^{\mathdefault{-3}}$)')
        S_DCN  = self.S['DCN']
        
        t = S_DCN['t']
        for c in ('H-1', 'H-4', 'H-5'):
            n = S_DCN[c]
            n[n < 0] = np.nan
            ax.plot(t, 1e-19*n, label=c)
        ax.legend()
        return ax

    def plot_XPR_I(self, ax):
        ax = get_axes(ax)
        ax.set_ylabel('Current (A)')
        XPR = self.XPR

        t = XPR['tip1'].t
        I1, I2, I3 = XPR['tip1'].x, XPR['tip2'].x, XPR['tip3'].x

        #I1 = ma.masked_array(I1, I1 < 0)
        #I2 = ma.masked_array(I2, I2 < 0)

        ax.plot(t, I1, label='Mach tip 1')
        ax.plot(t, I2, label='Mach tip 2')
        ax.plot(t, I3, label='Single tip')
        ax.legend()
        return ax

    def plot_XPR_R(self, ax):
        ax = get_axes(ax)
        ax.set_ylabel('Pos (cm)')
        XPR = self.XPR
        ax.plot(XPR['R'].t, 100*XPR['R'].x)
        return ax

    def plot_Ipolsol(self, ax):
        ax = get_axes(ax)
        ax.set_ylabel('Current (kA)')
        S_MAC = self.S['MAC']

        t = S_MAC['t']
        Ia = ma.masked_array(S_MAC['Ipolsola'], t > 6.)
        Ii = ma.masked_array(S_MAC['Ipolsoli'], t > 6.)

        ax.plot(t, 1e-3*Ia, label='Ipolsola')
        ax.plot(t, 1e-3*Ii, label='Ipolsoli')
        ax.legend()
        return ax

    def plot(self, fig=None, plots=None):
        if plots is None:
            plots = self.def_plots

        fig = get_fig(fig, shape=(len(plots), 1), figsize=(6,6), xlab='t (s)')
        fig.axes[0].set_xlim((1,7))

        for p, ax in zip(plots, fig.axes):
            getattr(self, 'plot_' + p)(ax)
        return fig

