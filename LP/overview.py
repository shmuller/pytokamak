import numpy as np
import numpy.ma as ma

from sig import memoized_property, Digitizer, Signal

from sm_pyplot.tight_figure import get_tfig, get_axes

from probe_xpr import TdiError, IOMdsAUG, IOFileAUG, ProbeXPR, ShotNotFoundError


class DigitizerAUG(Digitizer):
    def __init__(self, shn, diag, nodes, name=""):
        Digitizer.__init__(self, shn, name=name)

        self.IO_mds = IOMdsAUG(shn, diag=diag, raw=False)
        self.IO_file = IOFileAUG(shn, suffix='_AUG', group=name + '/' + diag)
        self.nodes = nodes + ('t',)

    def load(self, **kw):
        try:
            return Digitizer.load(self, **kw)
        except TdiError:
            return {node: np.zeros(1) for node in self.nodes}


AUG_diags = dict(
    dens = dict(diag='DCN', nodes=('H-1', 'H-2', 'H-3', 'H-4', 'H-5')),
    ncor = dict(diag='TOT', nodes=('H-1_corr', 'H-2_corr', 'H-3_corr', 'H-4_corr', 'H-5_corr')),
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
        
        (S['pech']['PECRH']*1e-6).plot(ax, label="ECRH")
        (S['pnbi']['PNI']*1e-6).plot(ax, label="NBI")
        (S['wmhd']['Wmhd']*1e-5).plot(ax, label="WMHD (x10)")
        ax.legend()
        return ax

    def plot_rad(self, ax):
        ax = get_axes(ax)
        ax.set_ylabel('Power (MW)')

        S = self.S['prad']['Pradtot']*1e-6
        S.masked(S.t > 6.).plot(ax)
        ax.legend()
        return ax

    def plot_power_rad(self, ax):
        ax = self.plot_power(ax)
        return self.plot_rad(ax)

    def plot_density(self, ax, chn=('H-1', 'H-4', 'H-5')):
        S = self.S['dens']
        ax = get_axes(ax)
        ax.set_ylabel('n (10$^{\mathdefault{19}}$ m$^{\mathdefault{-3}}$)')
        
        for c in chn:
            (S[c]*1e-19).masked(S[c] < 0).plot(ax)
        ax.legend()
        return ax

    def plot_n(self, ax, chn=('H-1_corr', 'H-4_corr', 'H-5_corr')):
        S = self.S['ncor']
        ax = get_axes(ax)
        ax.set_ylabel('n (10$^{\mathdefault{19}}$ m$^{\mathdefault{-3}}$)')
        
        for c in chn:
            (S[c]*1e-19).masked(S[c] < 0).plot(ax)
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
        
        I1, I2, I3 = XPR['tip1'], XPR['tip2'], XPR['tip3']

        if not no_Mach:
            I1.plot(ax, label='Mach tip 1')
            I2.plot(ax, label='Mach tip 2')
        if not no_single:
            I3.plot(ax, label='Single tip')
        ax.legend()
        return ax

    def plot_XPR_I_Mach(self, ax):
        return self.plot_XPR_I(ax, no_single=True)

    def plot_XPR_I_single(self, ax):
        return self.plot_XPR_I(ax, no_Mach=True)

    def plot_XPR_V(self, ax, no_Mach=False, no_single=False):
        ax = get_axes(ax)
        ax.set_ylabel('Voltage (V)')

        try:
            XPR = self.XPR
        except AttributeError:
            return ax

        V1, V2, V3 = XPR['tip1'].V, XPR['tip2'].V, XPR['tip3'].V

        if not no_Mach:
            V1.plot(ax, label='Mach tip 1')
            V2.plot(ax, label='Mach tip 2')
        if not no_single:
            V3.plot(ax, label='Single tip')
        ax.legend()
        return ax

    def plot_XPR_V_Mach(self, ax):
        return self.plot_XPR_V(ax, no_single=True)

    def plot_XPR_V_single(self, ax):
        return self.plot_XPR_V(ax, no_Mach=True)

    def plot_XPR_R(self, ax):
        ax = get_axes(ax)
        ax.set_ylabel('Pos (cm)')

        try:
            XPR = self.XPR
        except AttributeError:
            return ax

        (XPR['R']*100).plot(ax)
        return ax

    def plot_Ipolsol(self, ax):
        S = self.S['isol']
        ax = get_axes(ax)
        ax.set_ylabel('Div. cur. (kA)')

        Sa, Si = S['Ipolsola']*1e-3, S['Ipolsoli']*1e-3
        Sa.masked(Sa.t > 6.).plot(ax, label='Outer')
        Si.masked(Si.t > 6.).plot(ax, label='Inner')
        ax.legend()
        return ax

    def plot_Tdiv(self, ax):
        S = self.S['tdiv']['Tdiv']
        ax = get_axes(ax)
        ax.set_ylabel('Temp (eV)')

        S.plot(ax)
        ax.legend()
        return ax

    def plot_Da(self, ax):
        S = self.S['elmh']
        ax = get_axes(ax)
        ax.set_ylabel('Photons (au)')

        Sa, Si = S['ELMa-Han'], S['ELMi-Han']
        Sa.masked(Sa.t > 6.).plot(ax, label='Da outer')
        Si.masked(Si.t > 6.).plot(ax, label='Da inner')
        ax.legend()
        return ax

    def plot_gas(self, ax):
        S = self.S['gasv']['D_tot']*1e-21
        ax = get_axes(ax)
        ax.set_ylabel('Gas (10$^{\mathdefault{21}}$ el s$^{\mathdefault{-1}}$)')

        S.masked((S.t < 0.5) | (S.t > 6.)).plot(ax, label='D total')
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

