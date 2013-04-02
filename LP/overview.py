import numpy as np
import numpy.ma as ma

from sig import memoized_property, Digitizer

from sm_pyplot.tight_figure import get_fig, get_axes

from probe_xpr import TdiError, IOMdsAUG, IOFileAUG, ProbeXPR


class DigitizerAUG(Digitizer):
    def __init__(self, shn, diag, nodes, sock=None, name=""):
        Digitizer.__init__(self, shn, sock, name=name)

        self.IO_mds = IOMdsAUG(shn, sock, diag=diag, raw=False)
        self.IO_file = IOFileAUG(shn, suffix='_AUG', group=name + '/' + diag)
        self.nodes = nodes + ('t',)

    def load(self):
        try:
            return Digitizer.load(self)
        except TdiError:
            return {node: 0. for node in self.nodes}


AUG_diags = dict(
    dens = dict(diag='DCN', nodes=('H-1', 'H-2', 'H-3', 'H-4', 'H-5')),
    pnbi = dict(diag='NIS', nodes=('PNI',)),
    pech = dict(diag='ECS', nodes=('PECRH',)),
    wmhd = dict(diag='FPG', nodes=('Wmhd',)),
    isol = dict(diag='MAC', nodes=('Ipolsola', 'Ipolsoli')),
    tdiv = dict(diag='MAC', nodes=('Tdiv',)))


class AUGOverview:
    def __init__(self, shn):
        self.shn = shn

        self.XPR = ProbeXPR(shn=shn)

        self.def_plots = ('power', 'density', 'XPR_I', 'XPR_R', 'Ipolsol', 'Tdiv')

    @memoized_property
    def S(self):
        return {k: DigitizerAUG(self.shn, name=k, **v) for k, v in AUG_diags.iteritems()}

    def plot_power(self, ax):
        S = self.S
        ax = get_axes(ax)
        ax.set_ylabel('Power (MW)')
        
        ax.plot(S['pech']['t'], 1e-6*S['pech']['PECRH'], label="ECRH")
        ax.plot(S['pnbi']['t'], 1e-6*S['pnbi']['PNI'], label="NBI")
        ax.plot(S['wmhd']['t'], 1e-5*S['wmhd']['Wmhd'], label="WMHD (x10)")
        ax.legend()
        return ax

    def plot_density(self, ax):
        S = self.S
        ax = get_axes(ax)
        ax.set_ylabel('n (10$^{\mathdefault{19}}$ m$^{\mathdefault{-3}}$)')
        
        t = S['dens']['t']
        for c in ('H-1', 'H-4', 'H-5'):
            n = S['dens'][c]
            n[n < 0] = np.nan
            ax.plot(t, 1e-19*n, label=c)
        ax.legend()
        return ax

    def plot_XPR_I(self, ax, no_Mach=False, no_single=False):
        ax = get_axes(ax)
        ax.set_ylabel('Current (A)')
        XPR = self.XPR

        t = XPR['tip1'].t
        I1, I2, I3 = XPR['tip1'].x, XPR['tip2'].x, XPR['tip3'].x

        #I1 = ma.masked_array(I1, I1 < 0)
        #I2 = ma.masked_array(I2, I2 < 0)

        if not no_Mach:
            ax.plot(t, I1, label='Mach tip 1')
            ax.plot(t, I2, label='Mach tip 2')
        if not no_single:
            ax.plot(t, I3, label='Single tip')
        ax.legend()
        return ax

    def plot_XPR_V(self, ax, no_Mach=False, no_single=False):
        ax = get_axes(ax)
        ax.set_ylabel('Voltage (V)')
        XPR = self.XPR

        t = XPR['tip1'].t
        V1, V2, V3 = XPR['tip1'].V.x, XPR['tip2'].V.x, XPR['tip3'].V.x

        if not no_Mach:
            ax.plot(t, V1, label='Mach tip 1')
            ax.plot(t, V2, label='Mach tip 2')
        if not no_single:
            ax.plot(t, V3, label='Single tip')
        ax.legend()
        return ax

    def plot_XPR_R(self, ax):
        ax = get_axes(ax)
        ax.set_ylabel('Pos (cm)')
        XPR = self.XPR
        ax.plot(XPR['R'].t, 100*XPR['R'].x)
        return ax

    def plot_Ipolsol(self, ax):
        S = self.S
        ax = get_axes(ax)
        ax.set_ylabel('Current (kA)')

        t = S['isol']['t']
        Ia = ma.masked_array(S['isol']['Ipolsola'], t > 6.)
        Ii = ma.masked_array(S['isol']['Ipolsoli'], t > 6.)

        ax.plot(t, 1e-3*Ia, label='Ipolsola')
        ax.plot(t, 1e-3*Ii, label='Ipolsoli')
        ax.legend()
        return ax

    def plot_Tdiv(self, ax):
        S = self.S
        ax = get_axes(ax)
        ax.set_ylabel('Temp (eV)')

        t = S['tdiv']['t']
        ax.plot(t, S['tdiv']['Tdiv'], label='Tdiv')
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

