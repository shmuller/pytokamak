import numpy as np
import numpy.ma as ma

from sm_pyplot.tight_figure import get_tfig, get_axes

from LP.sig import memoized_property, Digitizer, Signal
from LP.probe_xpr import TdiError, IOMdsAUG, IOFileAUG, ProbeXPR, ShotNotFoundError

from equilibrium import Eqi

class DigitizerAUG(Digitizer):
    def __init__(self, shn, diag, nodes, subgrp='', **kw):
        Digitizer.__init__(self, shn, name=diag, **kw)

        self.IO_mds = IOMdsAUG(shn, diag=diag, raw=False)
        self.IO_file = IOFileAUG(shn, suffix='_AUG', group=diag + '/' + subgrp)
        if self.tnode not in nodes:
            nodes += (self.tnode,)
        self.nodes = nodes

    def load(self, **kw):
        try:
            return Digitizer.load(self, **kw)
        except TdiError:
            return {node: np.zeros(1) for node in self.nodes}

    def calib(self):
        for node, x in self.x.iteritems():
            if x.ndim > 1:
                self.x[node] = x.transpose(np.roll(np.arange(x.ndim), 1))
        Digitizer.calib(self)


class DigitizerAUGMAC(DigitizerAUG):
    def __init__(self, shn):
        DigitizerAUG.__init__(self, shn, diag='MAC', nodes=('Ipolsola', 'Ipolsoli'))
        self.dig_Tdiv = DigitizerAUG(shn, diag='MAC', subgrp='Tdiv', nodes=('Tdiv',))

    def __getitem__(self, indx):
        try:
            return DigitizerAUG.__getitem__(self, indx)
        except KeyError:
            return self.dig_Tdiv[indx]


class DigitizerAUGEQI(DigitizerAUG):
    def __init__(self, shn):
        DigitizerAUG.__init__(self, shn, diag='EQI', 
                              nodes=('Ri', 'Zj', 'PFM', 'ikCAT', 'RPFx', 'zPFx', 'PFxx'))
        self.mapspec = dict(magnaxis=0, xpoint=1, innerlim=2, xpoint2=3, outerlim=4)
   
    def get_R_z_psi(self):
        return self.x['Ri'][0], self.x['Zj'][0], self['PFM']

    def get_R_z_psi_special(self, spec):
        i = self.mapspec[spec]
        return self['RPFx'][:,i], self['zPFx'][:,i], self['PFxx'][:,i]


AUG_diags = dict(
    DCN = dict(nodes=('H-1', 'H-2', 'H-3', 'H-4', 'H-5')),
    TOT = dict(nodes=('H-1_corr', 'H-2_corr', 'H-3_corr', 'H-4_corr', 'H-5_corr')),
    NIS = dict(nodes=('PNI',)),
    ECS = dict(nodes=('PECRH',)),
    FPG = dict(nodes=('Wmhd',)),
    POT = dict(nodes=('ELMa-Han', 'ELMi-Han')),
    UVS = dict(nodes=('D_tot',)),
    BPD = dict(nodes=('Pradtot',)),
    MAG = dict(nodes=('Ipa', 'ULid12')),
    MHE = dict(nodes=('C09-23',), s=slice(None, None, 4)),
    CEZ = dict(nodes=('R', 'z', 'phi', 'vrot', 'Ti', 'inte', 
                      'err_vrot', 'err_Ti', 'err_inte')))


class AUGOverview:
    def __init__(self, shn):
        self.shn = shn

        try:
            self.XPR = ProbeXPR(shn=shn)
            self.def_plots = ('power', 'density', 'XPR_I', 'XPR_R', 'Ipolsol')
        except ShotNotFoundError:
            self.def_plots = ('power', 'density', 'Ipolsol')

        self.all_plots = self.def_plots + ('Tdiv', 'Da', 'gas')
        
    def __getitem__(self, indx):
        return self.S[indx]

    @memoized_property
    def S(self):
        S = {k: DigitizerAUG(self.shn, diag=k, **v) for k, v in AUG_diags.iteritems()}
        S['MAC'] = DigitizerAUGMAC(self.shn)
        S['EQI'] = DigitizerAUGEQI(self.shn)
        return S

    @memoized_property
    def eqi(self):
        return Eqi(self.S['EQI'])

    def plot_eqi(self, ax=None, Lvls=np.linspace(0., 1., 10)):
        ax = get_axes(ax)
        ax = self.eqi.get_flux_surf(3.45, Lvls).plot(ax)
        ax = self.eqi.get_separatrix(3.45).plot(ax, color='b', linewidth=2)
        return ax

    def plot_power(self, ax=None):
        S = self.S
        ax = get_axes(ax)
        ax.set_ylabel('Power (MW)')
        
        (S['ECS']['PECRH']*1e-6).plot(ax, label="ECRH")
        (S['NIS']['PNI']*1e-6).plot(ax, label="NBI")
        (S['FPG']['Wmhd']*1e-5).plot(ax, label="WMHD (x10)")
        ax.legend()
        return ax

    def plot_rad(self, ax=None):
        ax = get_axes(ax)
        ax.set_ylabel('Power (MW)')

        S = self.S['BPD']['Pradtot']*1e-6
        S.masked(S.t > 6.).plot(ax)
        ax.legend()
        return ax

    def plot_power_rad(self, ax=None):
        ax = self.plot_power(ax)
        return self.plot_rad(ax)

    def plot_density(self, ax=None, chn=('H-1', 'H-4', 'H-5')):
        S = self.S['DCN']
        ax = get_axes(ax)
        ax.set_ylabel('n (10$^{\mathdefault{19}}$ m$^{\mathdefault{-3}}$)')
        
        for c in chn:
            (S[c]*1e-19).masked(S[c] < 0).plot(ax)
        ax.legend()
        return ax

    def plot_n(self, ax=None, chn=('H-1_corr', 'H-4_corr', 'H-5_corr')):
        S = self.S['TOT']
        ax = get_axes(ax)
        ax.set_ylabel('n (10$^{\mathdefault{19}}$ m$^{\mathdefault{-3}}$)')
        
        for c in chn:
            (S[c]*1e-19).masked(S[c] < 0).plot(ax)
        ax.legend()
        return ax

    def plot_H1(self, ax=None):
        return self.plot_density(ax, chn=('H-1',))

    def plot_XPR_I(self, ax=None, no_Mach=False, no_single=False):
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

    def plot_XPR_I_Mach(self, ax=None):
        return self.plot_XPR_I(ax, no_single=True)

    def plot_XPR_I_single(self, ax=None):
        return self.plot_XPR_I(ax, no_Mach=True)

    def plot_XPR_V(self, ax=None, no_Mach=False, no_single=False):
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

    def plot_XPR_V_Mach(self, ax=None):
        return self.plot_XPR_V(ax, no_single=True)

    def plot_XPR_V_single(self, ax=None):
        return self.plot_XPR_V(ax, no_Mach=True)

    def plot_XPR_R(self, ax=None):
        ax = get_axes(ax)
        ax.set_ylabel('Pos (cm)')

        try:
            XPR = self.XPR
        except AttributeError:
            return ax

        (XPR['R']*100).plot(ax)
        return ax

    def plot_Ipolsol(self, ax):
        S = self.S['MAC']
        ax = get_axes(ax)
        ax.set_ylabel('Div. cur. (kA)')

        Sa, Si = S['Ipolsola']*1e-3, S['Ipolsoli']*1e-3
        Sa.masked(Sa.t > 6.).plot(ax, label='Outer')
        Si.masked(Si.t > 6.).plot(ax, label='Inner')
        ax.legend()
        return ax

    def plot_Tdiv(self, ax=None):
        S = self.S['MAC']['Tdiv']
        ax = get_axes(ax)
        ax.set_ylabel('Temp (eV)')

        S.plot(ax)
        ax.legend()
        return ax

    def plot_Da(self, ax=None):
        S = self.S['POT']
        ax = get_axes(ax)
        ax.set_ylabel('Photons (au)')

        Sa, Si = S['ELMa-Han'], S['ELMi-Han']
        Sa.masked(Sa.t > 6.).plot(ax, label='Da outer')
        Si.masked(Si.t > 6.).plot(ax, label='Da inner')
        ax.legend()
        return ax

    def plot_gas(self, ax=None):
        S = self.S['UVS']['D_tot']*1e-21
        ax = get_axes(ax)
        ax.set_ylabel('Gas (10$^{\mathdefault{21}}$ el s$^{\mathdefault{-1}}$)')

        S.masked((S.t < 0.5) | (S.t > 6.)).plot(ax, label='D total')
        ax.legend()
        return ax

    def plot_ipvl(self, ax=None):
        S = self.S['MAG']
        ax = get_axes(ax)
        Ip, Vl = S['Ipa']*1e-6, S['ULid12']

        Ip.plot(ax, label='Ip (MA)')
        Vl.masked((Vl.t < 0.5) | (Vl.t > 6.)).plot(ax, label='V loop (V)')
        ax.legend()
        return ax

    def plot_mirn(self, ax=None):
        S = self.S['MHE']['C09-23']
        ax = get_axes(ax)
        S.plot(ax)
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

