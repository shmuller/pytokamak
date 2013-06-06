import numpy as np
import numpy.ma as ma

from sm_pyplot.tight_figure import get_tfig, get_axes, show

from LP.sig import memoized_property, BoundingBox
from LP.probe_xpr import ProbeXPR, ShotNotFoundError

from digitizer_aug import DigitizerAUG, DigitizerAUGMAC, eqi_digitizers, dig_YGC
from equilibrium import Eqi, EqiViewer

from vtk_aug import VtkProxy, VtkWindow, VtkRotatingPolygon

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
    MHE = dict(nodes=('C09-23',), s=slice(None, None, 4)),
    CEZ = dict(nodes=('vrot', 'Ti', 'inte', 'err_vrot', 'err_Ti', 'err_inte',
                      'R', 'z', 'phi')))


class Vessel(VtkProxy):
    def __init__(self, digitizer):
        self.digitizer = digitizer

    @memoized_property
    def bdry(self):
        return np.concatenate(self.digitizer.xy_bdry)

    @memoized_property
    def vtk(self):
        kw = dict(phi1=3./2.*np.pi)
        return [VtkRotatingPolygon(*xy.T, **kw) for xy in self.digitizer.xy]

    def get_bbox(self):
        return BoundingBox(self.bdry.min(axis=0), self.bdry.max(axis=0))

    def plot(self, *args, **kw):
        return self.digitizer.plot(*args, **kw)

    def prerender(self, **kw):
        alpha = np.ones(len(self.vtk))
        alpha[[0, 1]] = 0.2
        alpha[[2, 4, 5]] = 0

        win = VtkWindow(**kw)
        for rpoly, a in zip(self.vtk, alpha):
            if a > 0:
                win = rpoly.prerender(win=win, alpha=a)
        return win


class EqiViewerAUG(EqiViewer):
    def viewer(self, event):
        EqiViewer.viewer(self, event, 
                         pos=(50, 150), figsize=(4.5, 6), xlab="R (m)", ylab="z (m)")


class EqiViewerAUGXPR(EqiViewerAUG):
    def __init__(self, eqi, XPR):
        EqiViewerAUG.__init__(self, eqi)
        self.XPR = XPR

    def plotfun(self, event):
        artists = EqiViewerAUG.plotfun(self, event)

        t_event = event.xdata
        ax = self.ax
        R, z = self.XPR.pos(t_event, masked=True).x[0]
        pp_head = self.XPR.head.as_path_patch(R, z)
        ax.add_patch(pp_head)
        artists.append(pp_head)

        if R < self.XPR.R0:
            fli = self.eqi.get_field_line_integrator(t_event)
            y0 = np.array([R, z, 0.])
            fl = fli.solve_bdry(y0)
            print fl
            pp_fl = fl.as_path_patch()
            ax.add_patch(pp_fl)
            artists.append(pp_fl)
        return artists


class ProfViewerAUG(ToggleViewer):
    def __init__(self, x, y):
        self.x, self.y = x, y

        ToggleViewer.__init__(self, menu_entry='Prof viewer')

    def plotfun(self, event):
        t_event = event.xdata
        x, y = self.x, self.y(t_event).x[0]
        return self.ax.plot(x, y, 'b')

    def viewer(self, event):
        fig = get_tfig(pos=(950, 150), figsize=(5, 5), 
                xlab="R (m)", ylab="vrot (km s$^{\mathdefault{-1}}$)")
        self.ax = fig.axes[0]
        self.ax.set_xlim((1.5, 2.5))
        self.ax.set_ylim((-50, 100))


class EqiViewerAUGVtk(ToggleViewerVtk):
    def __init__(self, eqi):
        self.eqi = eqi
        ToggleViewerVtk.__init__(self, menu_entry='VTK viewer')

    def viewer(self, event=None):
        self.ax = self.eqi.vessel.render()

    def plotfun(self, event):
        t_event = event.xdata
        win = self.ax
        FS = self.eqi.get_flux_surf(t_event)
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

        try:
            ProbeXPR.find_shot(shn=shn)
            self.def_plots = ('power', 'density', 'XPR_I', 'XPR_R', 'Ipolsol')
        except ShotNotFoundError:
            self.def_plots = ('power', 'density', 'Ipolsol')

        self.all_plots = self.def_plots + ('Tdiv', 'Da', 'gas')
        
    def __getitem__(self, indx):
        return self.S[indx]

    @memoized_property
    def S(self):
        S = {k: DigitizerAUG(self.shn, diag=k, **v) for k, v in aug_diags.iteritems()}
        S['MAC'] = DigitizerAUGMAC(self.shn)
        S['EQI'] = eqi_digitizers[self.eqi_dig](self.shn)
        S['YGC'] = dig_YGC
        return S

    @memoized_property
    def ves(self):
        return Vessel(self.S['YGC'])

    @memoized_property
    def eqi(self):
        return EqiAUG(self.S['EQI'], self.ves)

    @memoized_property
    def XPR(self):
        return ProbeXPR(shn=self.shn, eqi=self.eqi)

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

    def plot_CER(self, ax=None, name='vrot', ylab='vrot (km s$^{\mathdefault{-1}}$)'):
        S = self.S['CEZ']
        ax = get_axes(ax)
        ax.set_ylabel(ylab)
        
        i = [0, 14, 23]
        R = S['R'].x[i]
        S = S[name][:,i]*1e-3
        S.plot(ax)
        ax.legend(['R = %.2f m' % x for x in R])
        return ax

    def plot_vrot(self, ax=None):
        return self.plot_CER(ax, 'vrot', 'vrot (km s$^{\mathdefault{-1}}$)')

    def plot_Ti(self, ax=None):
        return self.plot_CER(ax, 'Ti', 'Ti (keV)')

    def plot(self, plots=None, fig=None):
        if plots is None:
            plots = self.def_plots

        try:
            self.viewers = [EqiViewerAUGXPR(self.eqi, self.XPR)]
        except AttributeError:
            self.viewers = [EqiViewerAUG(self.eqi)]

        try:
            self.viewers += [EqiViewerAUGVtk(self.eqi)]
        except AttributeError:
            pass

        try:
            S = self.S['CEZ']
            self.viewers += [ProfViewerAUG(S['R'].x[:24], S['vrot'][:,:24]*1e-3)]
        except:
            pass

        fig = get_tfig(fig, pos=(450, 150), figsize=(6,6), shape=(len(plots), 1), 
                       xlab='t (s)', viewers=self.viewers)
        fig.axes[0].set_xlim((1,7))

        for p, ax in zip(plots, fig.axes):
            getattr(self, 'plot_' + p)(ax)
        return fig

    def plot_all(self, **kw):
        return self.plot(self.all_plots, **kw)


if __name__ == "__main__":
    AUG = AUGOverview(shn=30017)
    FS = AUG.eqi.get_flux_surf(2.5)
    fli = AUG.eqi.get_field_line_integrator(2.5)
        
    R0, l = fli.test()
    
    fl = fli(1.475, -0.966)
    fl2 = fli(1.5, 0.)

    #ax = AUG.ves.plot()
    #FS.plot(ax=ax)
    #fl.plot(ax=ax, linewidth=2)
    #fl2.plot(ax=ax, linewidth=2)

    win = AUG.ves.prerender()
    FS.prerender(win=win)
    fl.prerender(win=win)
    fl2.prerender(win=win, color=(1., 0., 0.))
    win.render()

    win.show()
    show()
