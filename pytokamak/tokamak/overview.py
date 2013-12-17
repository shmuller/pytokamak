import numpy as np
import numpy.ma as ma

from pytokamak.utils.utils import memoized_property, BoundingBox, Timer, tic, toc
from pytokamak.utils.sig import Signal, get_tfig, get_axes, show
from pytokamak.LP.probe_xpr import ProbeXPR, ShotNotFoundError

from digitizer_aug import DigitizerAUG, DigitizerAUGMAC, DigitizerAUGMIR, eqi_digitizers
from diaggeom_aug import Vessel
from equilibrium import Eqi2 as Eqi
from equilibrium import EqiViewer
from diags_aug import DCN

from sm_pyplot.observer_viewer import ToggleViewer, ToggleViewerVtk

aug_diags = dict(
    DCN = dict(nodes=('H-1', 'H-2', 'H-3', 'H-4', 'H-5')),
    TOT = dict(nodes=('H-1_corr', 'H-2_corr', 'H-3_corr', 'H-4_corr', 'H-5_corr')),
    NIS = dict(nodes=('PNI',)),
    ECS = dict(nodes=('PECRH',)),
    FPG = dict(nodes=('Wmhd',)),
    POT = dict(nodes=('ELMa-Han', 'ELMi-Han')),
    UVS = dict(nodes=('D_tot',)),
    BPD = dict(nodes=('Pradtot',)),
    MAG = dict(nodes=('Ipa', 'ULid12')),
    CEZ = dict(nodes=('vrot', 'Ti', 'inte', 'err_vrot', 'err_Ti', 'err_inte',
                      'R', 'z', 'phi'), tnode='time'),
    XVS = dict(nodes=('S2L1A10',)))


class EqiViewerAUG(EqiViewer):
    def viewer(self, event=None, **kw):
        return EqiViewer.viewer(self, event, 
                                pos=(50, 150), figsize=(4.5, 6), 
                                xlab="R (m)", ylab="z (m)", **kw)


class EqiViewerAUGXPR(EqiViewerAUG):
    def __init__(self, eqi, XPR):
        EqiViewerAUG.__init__(self, eqi)
        self.XPR = XPR
    
    def viewer(self, event=None, **kw):
        # cache 'pos' when viewer is active
        #self.pos = self.XPR.pos
        self.pos = self.XPR.pos.as_spline(k=2)
        return EqiViewerAUG.viewer(self, event, **kw)
    
    def on_close(self, event=None):
        # remove 'pos' when viewer is inactive
        del self.pos
        return EqiViewerAUG.on_close(self, event)

    def plotfun(self, event):
        artists = EqiViewerAUG.plotfun(self, event)
        
        t_event = event.xdata
        ax = self.ax
        #R, z = self.pos(t_event, masked=True).x[0]
        R, z = self.pos(t_event).ravel()
        pp_head = self.XPR.head.as_path_patch(R, z)
        ax.add_patch(pp_head)
        artists.append(pp_head)

        if R < self.XPR.R0:
            y0 = np.array([R, z, 0.])
            fl = self.fli.solve_bdry(y0)
            print fl
            pp_fl = fl.as_path_patch()
            ax.add_patch(pp_fl)
            artists.append(pp_fl)
        return artists


class ProfViewerAUG(ToggleViewer):
    def __init__(self, x, y, menu_entry, **kw):
        self.x, self.y, self.kw = x, y, kw
        ToggleViewer.__init__(self, menu_entry=menu_entry)

    def plotfun(self, event):
        t_event = event.xdata
        x, y = self.x, self.y(t_event).x[0]
        return self.ax.plot(x, y, 'b')

    def viewer(self, event=None, **kw):
        kw2 = self.kw.copy()
        kw2.update(kw)
        fig = get_tfig(pos=(950, 150), figsize=(5, 5), **kw2)
        self.ax = ax = fig.axes[0]
        ax.set_xlim((1.5, 2.5))
        ax.set_ylim((-50, 100))
        return ax


class MIRViewerAUG(ToggleViewer):
    def __init__(self, dig):
        self.x = np.array([int(node[4:]) for node in dig.nodes])
        t = dig[dig.nodes[0]].t
        y = np.zeros((t.size, self.x.size), 'f')
        for i, node in enumerate(dig.nodes):
            y[:,i] = dig[node].x
        self.y = Signal(y, t)
        ToggleViewer.__init__(self, menu_entry='MIR viewer')

    def plotfun(self, event):
        t_event = event.xdata
        x, y = self.x, self.y(t_event).x[0]
        return self.ax.plot(x, y, 'b-+')

    def viewer(self, event=None, **kw):
        fig = get_tfig(pos=(950, 150), figsize=(5, 5), xlab="Number", **kw)
        self.ax = ax = fig.axes[0]
        ax.set_xlim((0, 33))
        ax.set_ylim((-2, 2))
        return ax


class EqiViewerAUGVtk(ToggleViewerVtk):
    def __init__(self, eqi):
        self.eqi = eqi
        ToggleViewerVtk.__init__(self, menu_entry='VTK viewer')

    def viewer(self, event=None, **kw):
        self.refine = kw.pop('refine', 1)
        self.ax = ax = self.eqi.vessel.render()
        return ax

    def plotfun(self, event):
        t_event = event.xdata
        win = self.ax
        FS = self.eqi.get_flux_surf(t_event, norm=True, refine=self.refine)
        FS.render(win=win)
        return [win.ren.GetActors().GetLastActor()]


class EqiAUG(Eqi):
    def get_viewers(self, XPR=None):
        if XPR is None:
            self.viewers = [EqiViewerAUG(self), EqiViewerAUGVtk(self)]
        else:
            self.viewers = [EqiViewerAUGXPR(self, XPR), EqiViewerAUGVtk(self)]
        return self.viewers


class AUGOverview:
    def __init__(self, shn, eqi_dig='EQI'):
        self.shn, self.eqi_dig = shn, eqi_dig

        self.ves = Vessel()
        self.eqi = EqiAUG(self.S['EQI'], self.ves)
        self.DCN = DCN(shn=shn, eqi=self.eqi)
        try:
            self.XPR = ProbeXPR(shn=shn, eqi=self.eqi)
            self.def_plots = ('power', 'density', 'XPR_I_Mach', 'XPR_R', 'Ipolsol')
        except ShotNotFoundError:
            self.XPR = None
            self.def_plots = ('power', 'density', 'Ipolsol')
        self.eqi.get_viewers(self.XPR)

    def __getitem__(self, indx):
        return self.S[indx]

    @memoized_property
    def S(self):
        S = {k: DigitizerAUG(self.shn, diag=k, **v) for k, v in aug_diags.iteritems()}
        S['MAC'] = DigitizerAUGMAC(self.shn)
        S['MIR'] = DigitizerAUGMIR(self.shn)
        S['EQI'] = eqi_digitizers[self.eqi_dig](self.shn)
        return S

    def plot_power(self, ax=None, **kw):
        S = self.S
        ax = get_axes(ax)
        ax.set_ylabel('Power (MW)')
        
        (S['ECS']['PECRH']*1e-6).plot(ax, label="ECRH")
        (S['NIS']['PNI']*1e-6).plot(ax, label="NBI")
        (S['FPG']['Wmhd']*1e-5).plot(ax, label="WMHD (x10)")
        ax.legend()
        return ax

    def plot_rad(self, ax=None, **kw):
        S = self.S['BPD']
        ax = get_axes(ax)
        ax.set_ylabel('Power (MW)')

        (S['Pradtot']*1e-6).t_lt(6.).plot(ax)
        ax.legend()
        return ax

    def plot_power_rad(self, ax=None, **kw):
        ax = self.plot_power(ax, **kw)
        return self.plot_rad(ax, **kw)

    def plot_density(self, ax=None, chn=('H-1', 'H-4', 'H-5'), **kw):
        S = self.DCN
        ax = get_axes(ax)
        ax.set_ylabel(r'$\int$n dl (10$^{\text{19}}$ m$^{\text{-2}}$)')
        
        for c in chn:
            (S[c]*1e-19).plot(ax, **kw)
        ax.legend()
        return ax

    def plot_DCN_sep(self, ax=None, chn=('H-1', 'H-4', 'H-5'), field='d', **kw):
        S = self.DCN
        ax = S.plot_separatrix(ax, chn=chn, field=field)
        ax.set_ylabel(S.sep[chn[0]][field].ylab)
        ax.legend()
        return ax

    def plot_DCN_Rsep(self, ax=None, **kw):
        return self.plot_DCN_sep(ax, field='R', **kw)

    def plot_DCN_zsep(self, ax=None, **kw):
        return self.plot_DCN_sep(ax, field='z', **kw)

    def plot_DCN_dsep(self, ax=None, **kw):
        return self.plot_DCN_sep(ax, field='d', **kw)

    def plot_n(self, ax=None, chn=('H-1_corr', 'H-4_corr', 'H-5_corr'), **kw):
        S = self.S['TOT']
        ax = get_axes(ax)
        ax.set_ylabel(r'n (10$^{\text{19}}$ m$^{\text{-3}}$)')
        
        for c in chn:
            (S[c]*1e-19).nonneg().plot(ax, **kw)
        ax.legend()
        return ax

    def plot_H1(self, ax=None, **kw):
        return self.plot_density(ax, chn=('H-1',), **kw)

    def plot_XPR_I(self, ax=None, **kw):
        return self.XPR.plot_I(ax, **kw)

    def plot_XPR_I_Mach(self, ax=None, **kw):
        return self.XPR.plot_I_Mach(ax, **kw)

    def plot_XPR_I_single(self, ax=None, **kw):
        return self.XPR.plot_I_single(ax, **kw)

    def plot_XPR_V(self, ax=None, **kw):
        return self.XPR.plot_V(ax, **kw)

    def plot_XPR_V_Mach(self, ax=None, **kw):
        return self.XPR.plot_V_Mach(ax, **kw)

    def plot_XPR_V_single(self, ax=None, **kw):
        return self.XPR.plot_V_single(ax, **kw)

    def plot_XPR_R(self, ax=None, legend_loc=None, **kw):
        return self.XPR.plot_R(ax, legend_loc=legend_loc, **kw)

    def plot_Ipolsol(self, ax, **kw):
        S = self.S['MAC']
        ax = get_axes(ax)
        ax.set_ylabel('Div. cur. (kA)')

        (S['Ipolsola']*1e-3).t_lt(6.).plot(ax, label='Outer')
        (S['Ipolsoli']*1e-3).t_lt(6.).plot(ax, label='Inner')
        ax.legend()
        return ax

    def plot_Tdiv(self, ax=None, **kw):
        S = self.S['MAC']['Tdiv']
        ax = get_axes(ax)
        ax.set_ylabel('Temp (eV)')

        S.plot(ax)
        ax.legend()
        return ax

    def plot_Da(self, ax=None, **kw):
        S = self.S['POT']
        ax = get_axes(ax)
        ax.set_ylabel('Photons (au)')

        S['ELMa-Han'].t_lt(6.).plot(ax, label='Da outer')
        S['ELMi-Han'].t_lt(6.).plot(ax, label='Da inner')
        ax.legend()
        return ax

    def plot_gas(self, ax=None, **kw):
        S = self.S['UVS']
        ax = get_axes(ax)
        ax.set_ylabel(r'Gas (10$^{\text{21}}$ el s$^{\text{-1}}$)')

        (S['D_tot']*1e-21).t_between(0.5, 6.).plot(ax, label='D total')
        ax.legend()
        return ax

    def plot_ipvl(self, ax=None, **kw):
        S = self.S['MAG']
        ax = get_axes(ax)

        (S['Ipa']*1e-6).plot(ax, label='Ip (MA)')
        S['ULid12'].t_between(0.5, 6.).plot(ax, label='V loop (V)')
        ax.legend()
        return ax

    def plot_mirn(self, ax=None, **kw):
        S = self.S['MIR']['C09-23']
        ax = get_axes(ax)
        S.plot(ax)
        ax.legend()
        return ax

    def plot_CER(self, ax=None, name='vrot', 
                 ylab=r'vrot (km s$^{\text{-1}}$)', **kw):
        S = self.S['CEZ']
        ax = get_axes(ax)
        ax.set_ylabel(ylab)
        
        i = [0, 14, 23]
        R = S['R'].x[i]
        S = S[name][:,i]*1e-3
        S.plot(ax)
        ax.legend(['R = %.2f m' % x for x in R])
        return ax

    def plot_vrot(self, ax=None, **kw):
        return self.plot_CER(ax, 'vrot', r'vrot (km s$^{\text{-1}}$)')

    def plot_Ti(self, ax=None, **kw):
        return self.plot_CER(ax, 'Ti', 'Ti (keV)')

    def plot_AXUV(self, ax=None, **kw):
        S = self.S['XVS']['S2L1A10']
        ax = get_axes(ax)
        S.plot(ax)
        ax.legend()
        return ax

    def specgram(self, **kw):
        S_list = [self.XPR['tip1+tip2'],
                  self['MAC']['Ipolsoli'],
                  self['MIR']['C09-23'],
                  self['XVS']['S2L1A10']]

        fig = get_tfig(shape=(5, 1), figsize=(6, 8), xlab='t (s)')
        self.plot_XPR_I_Mach(ax=fig.axes[0])

        for ax, S in zip(fig.axes[1:], S_list):
            S.specgram(ax=ax, **kw)
            ax.set_ylabel('%s f (kHz)' % S.name)
        return fig

    def figure(self, fig=None, shape=(5, 1), figsize=(6, 6), 
               CEZ_viewer=False, MIR_viewer=False, **kw):
        self.viewers = self.eqi.viewers[:]
        if CEZ_viewer:
            S = self.S['CEZ']
            y, x = S['vrot'][:,:24]*1e-3, S['R'].x[:24]
            viewer_kw = dict(xlab='R (m)', ylab=r'vrot (km s$^{\text{-1}}$)')
            CEZ_viewer = ProfViewerAUG(x, y, menu_entry='CEZ viewer', **viewer_kw)
            self.viewers.append(CEZ_viewer)

        if MIR_viewer:
            MIR_viewer = MIRViewerAUG(self.S['MIR'])
            self.viewers.append(MIR_viewer)

        fig = get_tfig(fig, pos=(450, 150), shape=shape, figsize=figsize, 
                       xlab='t (s)', viewers=self.viewers, **kw)
        fig.axes[0].set_xlim((1, 7))
        return fig

    def plot(self, plots=None, fig=None, kw_plots=None, **kw):
        if plots is None:
            plots = self.def_plots
        if kw_plots is None:
            kw_plots = dict()

        fig = self.figure(fig, shape=(len(plots), 1), **kw)

        for p, ax in zip(plots, fig.axes):
            getattr(self, 'plot_' + p)(ax, **kw_plots)
        return fig


if __name__ == "__main__":
    AUG = AUGOverview(shn=30017)
    FS = AUG.eqi.get_flux_surf(2.5)
    fli = AUG.eqi.get_field_line_integrator(2.5)
        
    R0, l = fli.test()
    
    fl = fli(1.475, -0.966)
    fl2 = fli(1.5, 0.)

    ax = AUG.ves.plot()
    FS.plot(ax=ax)
    fl.plot(ax=ax, linewidth=2)
    fl2.plot(ax=ax, linewidth=2)
    ax.figure.show()

    """
    win = AUG.ves.prerender()
    FS.prerender(win=win)
    fl.prerender(win=win)
    fl2.prerender(win=win, color=(1., 0., 0.))
    win.render()
    win.start()
    """
