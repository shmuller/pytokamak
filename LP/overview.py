import numpy as np
import numpy.ma as ma

from sig import memoized_property, Digitizer, Signal

from sm_pyplot.tight_figure import get_tfig, get_axes

from probe_xpr import TdiError, IOMdsAUG, IOFileAUG, ProbeXPR, ShotNotFoundError


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
            return {node: np.zeros(1) for node in self.nodes}


AUG_diags = dict(
    dens = dict(diag='DCN', nodes=('H-1', 'H-2', 'H-3', 'H-4', 'H-5')),
    pnbi = dict(diag='NIS', nodes=('PNI',)),
    pech = dict(diag='ECS', nodes=('PECRH',)),
    wmhd = dict(diag='FPG', nodes=('Wmhd',)),
    isol = dict(diag='MAC', nodes=('Ipolsola', 'Ipolsoli')),
    tdiv = dict(diag='MAC', nodes=('Tdiv',)),
    elmh = dict(diag='POT', nodes=('ELMa-Han', 'ELMi-Han')),
    gasv = dict(diag='UVS', nodes=('D_tot',)),
    prad = dict(diag='BPD', nodes=('Pradtot',)))


class AUGOverview:
    def __init__(self, shn):
        self.shn = shn

        try:
            self.XPR = ProbeXPR(shn=shn)
            self.def_plots = ('power', 'density', 'XPR_I', 'XPR_R', 'Ipolsol')
        except ShotNotFoundError:
            self.def_plots = ('power', 'density', 'Ipolsol')

        self.all_plots = self.def_plots + ('Tdiv', 'Da', 'gas')
        
    @memoized_property
    def S(self):
        return {k: DigitizerAUG(self.shn, name=k, **v) for k, v in AUG_diags.iteritems()}

    def plot_power(self, ax):
        S = self.S
        ax = get_axes(ax)
        ax.set_ylabel('Power (MW)')
        
        Signal(1e-6*S['pech']['PECRH'], S['pech']['t'], name="ECRH").plot(ax)
        Signal(1e-6*S['pnbi']['PNI'], S['pnbi']['t'], name="NBI").plot(ax)
        Signal(1e-5*S['wmhd']['Wmhd'], S['wmhd']['t'], name="WMHD (x10)").plot(ax)

        t = S['prad']['t']
        mask = t > 6.
        Signal(1e-6*S['prad']['Pradtot'], t, name="Pradtot").masked(mask).plot(ax)
        ax.legend()
        return ax

    def plot_density(self, ax, chn=('H-1', 'H-4', 'H-5')):
        S = self.S['dens']
        ax = get_axes(ax)
        ax.set_ylabel('n (10$^{\mathdefault{19}}$ m$^{\mathdefault{-3}}$)')
        
        t = S['t']
        for c in chn:
            Signal(1e-19*S[c], t, name=c).masked(S[c] < 0).plot(ax)
        ax.legend()
        return ax

    def plot_H1(self, ax):
        return self.plot_density(ax, chn=('H-1',))

    def plot_XPR_I(self, ax, no_Mach=False, no_single=False):
        ax = get_axes(ax)
        ax.set_ylabel('Current (A)')
        
        try:
            XPR = self.XPR
        except AttributeError:
            return ax
                
        t = XPR['tip1'].t
        I1, I2, I3 = XPR['tip1'].x, XPR['tip2'].x, XPR['tip3'].x

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

        try:
            XPR = self.XPR
        except AttributeError:
            return ax

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

        try:
            XPR = self.XPR
        except AttributeError:
            return ax

        ax.plot(XPR['R'].t, 100*XPR['R'].x)
        return ax

    def plot_Ipolsol(self, ax):
        S = self.S['isol']
        ax = get_axes(ax)
        ax.set_ylabel('Current (kA)')

        t = S['t']
        m = t > 6.
        Signal(1e-3*S['Ipolsola'], t, name='Ipolsola').masked(m).plot(ax)
        Signal(1e-3*S['Ipolsoli'], t, name='Ipolsoli').masked(m).plot(ax)
        ax.legend()
        return ax

    def plot_Tdiv(self, ax):
        S = self.S['tdiv']
        ax = get_axes(ax)
        ax.set_ylabel('Temp (eV)')

        t = S['t']
        Signal(S['Tdiv'], t, name='Tdiv').plot(ax)
        ax.legend()
        return ax

    def plot_Da(self, ax):
        S = self.S['elmh']
        ax = get_axes(ax)
        ax.set_ylabel('Photons (au)')

        t = S['t']
        m = t > 6.
        Signal(S['ELMa-Han'], t, name='Da outer').masked(m).plot(ax)
        Signal(S['ELMi-Han'], t, name='Da inner').masked(m).plot(ax)
        ax.legend()
        return ax

    def plot_gas(self, ax):
        S = self.S
        ax = get_axes(ax)
        ax.set_ylabel('Gas (10$^{\mathdefault{21}}$ el s$^{\mathdefault{-1}}$)')

        t = S['gasv']['t']
        m = (t < 0.5) | (t > 6.)
        Signal(1e-21*S['gasv']['D_tot'], t, name='D total').masked(m).plot(ax)
        ax.legend()
        return ax

    def plot(self, plots=None, fig=None):
        if plots is None:
            plots = self.def_plots

        fig = get_tfig(fig, shape=(len(plots), 1), figsize=(6,6), xlab='t (s)')
        fig.axes[0].set_xlim((1,7))

        for p, ax in zip(plots, fig.axes):
            getattr(self, 'plot_' + p)(ax)
        return fig

    def plot_all(self, **kw):
        return self.plot(self.all_plots, **kw)

