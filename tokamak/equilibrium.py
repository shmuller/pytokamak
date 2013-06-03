import math
import numpy as np

from scipy.ndimage import map_coordinates

from LP.sig import memoized_property
from LP.splines import Spline, Spline2D

from sm_pyplot.tight_figure import get_tfig, get_axes, show
from sm_pyplot.observer_viewer import ToggleViewer, ToggleViewerIntegrated, \
                                      ToggleViewerVtk
from vtk_aug import VtkProxy, VtkContour, VtkPolyline

from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.collections import PathCollection

from fieldline._fieldline import solve_bdry

class Interpolator:
    def __init__(self, t, x, y, f):
        self.t, self.x, self.y, self.f = t, x, y, f
        self.l, self.m, self.n = self.shape = f.shape


class Interpolator3DSlab(Interpolator):    
    def eval(self, t, x, y):
        coo = self.as_coo(t, self.t), self.as_coo(y, self.y), self.as_coo(x, self.x)
        return map_coordinates(self.f, coo)

    @staticmethod
    def as_coo(X, x):
        return (X - x[0]) / (x[-1] - x[0]) * (x.size - 1)


class InterpolatorSlice(Interpolator):
    @memoized_property
    def _splines(self):
        return np.array(None).repeat(self.l)

    def get_spline(self, i):
        sp = self._splines[i]
        if sp is None:
            sp = self._splines[i] = Spline2D(self.y, self.x, self.f[i])
        return sp

    def eval_path_on_slices(self, x, y):
        fs = np.empty((self.l, x.size))
        for i in xrange(self.l):
            sp = self.get_spline(i)
            fs[i] = sp.ev(y, x)
        return fs

    def get_path_spline(self, s, x, y):
        fs = self.eval_path_on_slices(x, y)
        return Spline2D(self.t, s, fs)


class FluxSurf(VtkContour):
    def __init__(self, x, y, f, lvls):
        VtkContour.__init__(self, x, y, f, lvls, plane='xz')

    def plot(self, ax=None, **kw):
        kw.setdefault('edgecolors', 'r')
        ax = get_axes(ax)
        pc = self.as_path_collection(**kw)
        ax.add_collection(pc)
        return ax


class FieldLine(VtkProxy):
    def __init__(self, y, t):
        self.y, self.t = y, t

        self.start, self.end = y[[0, -1], :2]
        self.n_turns = t[-1]
        self.length = y[-1, 2]

    @memoized_property
    def vtk(self):
        R, z = self.y[:,:2].T
        phi = 2.*np.pi*self.t
        return VtkPolyline(R*np.cos(phi), R*np.sin(phi), z, 
                           color=(0., 0., 1.), linewidth=1.5)

    def plot(self, ax=None, **kw):
        ax = get_axes(ax)
        ax.plot(self.y[:, 0], self.y[:, 1], **kw)
        return ax

    def as_path_collection(self, **kw):
        kw.setdefault('facecolors', 'none')
        kw.setdefault('edgecolors', 'b')
        kw.setdefault('linewidth', 2)
        paths = [Path(self.y[:,:2])]
        return PathCollection(paths, **kw)

    def as_path_patch(self, **kw):
        kw.setdefault('facecolor', 'none')
        kw.setdefault('edgecolor', 'b')
        kw.setdefault('linewidth', 2)
        path = Path(self.y[:,:2])
        return PathPatch(path, **kw)

    def __repr__(self):
        return "%s at (%.2f,% .2f) m" % (self.__class__.__name__, 
                                          self.start[0], self.start[1])

    def __str__(self):
        return self.__repr__() + ": " + \
               "n_turns={s.n_turns:5.2f}, l={s.length:6.2f} m".format(s=self)


class OpenFieldLine(FieldLine):
    def __init__(self, co, ctr):
        self.co, self.ctr = co, ctr
        self.y = np.concatenate((ctr.y[:0:-1], co.y))
        self.t = np.concatenate((ctr.t[:0:-1], co.t))

        self.start = co.start
        self.n_turns = co.n_turns - ctr.n_turns
        self.length = co.length - ctr.length

    def __str__(self):
        return self.__repr__() + ": " + \
               "n_turns=(%5.2f,%5.2f), l=(%6.2f,%6.2f) m" % \
               (self.co.n_turns, self.ctr.n_turns, self.co.length, self.ctr.length)


class FieldLineIntegrator:
    def __init__(self, splR, splz, bdry=None, bbox=None):
        if bdry is None:
            if bbox is None:
                bbox = splR.get_bbox()
                z0, R0 = bbox.x0
                z1, R1 = bbox.x1
            else:
                R0, z0 = bbox.x0
                R1, z1 = bbox.x1
            bdry = np.array([(R0, z0), (R1, z0), (R1, z1), (R0, z1)])

        self.splR, self.splz = splR, splz
        self.bdry = np.ascontiguousarray(bdry.T.ravel(), 'd')
    
    def __call__(self, R0, z0, **kw):
        y0 = np.array((R0, z0, 0.), np.float64)
        return self.solve_bdry(y0, **kw)

    def _solve_bdry(self, y0, t):
        y = np.zeros((t.size, y0.size))
        y[0] = y0
        points_done = solve_bdry(self.splR.astuple(), self.splz.astuple(), 
                                 y, t, self.bdry)
        return FieldLine(y[:points_done], t[:points_done])

    def solve_bdry(self, y0, n_turns=20., n=1000):
        t = np.linspace(0., n_turns, n)
        fl = self._solve_bdry(y0, t)
        if fl.t.size < n:
            # field line hit the wall, so go the other way too
            t = np.linspace(0., -n_turns, n)
            fl_ctr = self._solve_bdry(y0, t)
            fl = OpenFieldLine(fl, fl_ctr)
        return fl

    def test(self, n_turns=20., n=1000):
        R0 = np.linspace(1.29, 1.64, 100)
        y0 = np.array([0., -0.966, 0.])
        l = np.zeros((R0.size, 2))
        solver = self.solve_bdry

        for i in xrange(R0.size):
            y0[0] = R0[i]
            fl = solver(y0, n_turns, n)
            l[i, 0] = fl.co.length
            l[i, 1] = fl.ctr.length
        return R0, l


class Eqi:
    def __init__(self, digitizer, vessel=None):
        self.digitizer, self.vessel = digitizer, vessel

    def __getitem__(self, indx):
        return self.digitizer[indx]

    @memoized_property
    def R_z_psi(self):
        return self.digitizer.get_R_z_psi()

    @memoized_property
    def R(self):
        return self.R_z_psi[0].astype(np.float64)

    @memoized_property
    def z(self):
        return self.R_z_psi[1].astype(np.float64)

    @memoized_property
    def psi(self):
        return self.R_z_psi[2].astype(np.float64)

    @memoized_property
    def t(self):
        return self.psi.t

    @memoized_property
    def shape(self):
        return self.psi.shape

    @memoized_property
    def psi0(self):
        return self.digitizer.get_R_z_psi_special('magnaxis')[2]

    @memoized_property
    def psi1(self):
        return self.digitizer.get_R_z_psi_special('xpoint')[2]

    @memoized_property
    def psi_n(self):
        dpsi = self.psi1 - self.psi0
        return (self.psi - self.psi0[:, None, None]) / dpsi[:, None, None]

    def _property_factory(name):
        @memoized_property
        def prop(self):
            return self[name]
        return prop

    psii = _property_factory('psii')

    @memoized_property
    def psii_n(self):
        dpsi = self.psi1 - self.psi0
        return (self.psii - self.psi0[:, None]) / dpsi[:, None]

    def _get_flux_surf(self, f, lvls=None):
        return FluxSurf(self.R, self.z, f, lvls)

    def get_flux_surf(self, ti, lvls=None):
        if lvls is None:
            lvls = np.linspace(0., 1., 20)
        return self._get_flux_surf(self.psi_n(ti).x[0], lvls)

    def get_flux_surf_unnorm(self, ti, lvls=None):
        return self._get_flux_surf(self.psi(ti).x[0], lvls)

    def get_separatrix(self, ti):
        return self.get_flux_surf(ti, lvls=np.array([1.]))

    @memoized_property
    def interpolator_3d_slab(self):
        return Interpolator3DSlab(self.t, self.R, self.z, self.psi_n.x)

    def eval(self, t, R, z):
        return self.interpolator_3d_slab.eval(t, R, z)

    @memoized_property
    def interpolator_slice(self):
        return InterpolatorSlice(self.t, self.R, self.z, self.psi_n.x)
    
    def get_path_spline(self, s, R, z):
        return self.interpolator_slice.get_path_spline(s, R, z)

    def _psi(self, ti):
        R, z, psi = self.R, self.z, self.psi(ti).x[0]
        return Spline2D(z, R, psi)

    def _grad_psi(self, ti):
        psi = self._psi(ti)
        return psi.derivy(), psi.derivx()

    def get_psi(self, ti, R, z):
        return self._psi(ti).ev(z, R)
    
    def _check_grid(self, R, z):
        if R is None:
            R = self.R
        if z is None:
            z = self.z
        return R, z

    def get_psi_grid(self, ti, R=None, z=None):
        R, z = self._check_grid(R, z)
        return self._psi(ti).eval(z, R)

    def get_Bpol(self, ti, R, z):
        dRpsi, dzpsi = self._grad_psi(ti)
        fact = 1./R
        BR = -dzpsi.ev(z, R) * fact
        Bz =  dRpsi.ev(z, R) * fact
        return BR, Bz

    def get_Bpol_grid(self, ti, R=None, z=None):
        R, z = self._check_grid(R, z)
        dRpsi, dzpsi = self._grad_psi(ti)
        fact = (1./R)[None]
        BR = -dzpsi.eval(z, R) * fact
        Bz =  dRpsi.eval(z, R) * fact
        return BR, Bz

    def get_Bpol_spline(self, ti, R=None, z=None):
        R, z = self._check_grid(R, z)
        Bpol = self.get_Bpol_grid(R, z)
        splR = Spline2D(z, R, Bpol[0])
        splz = Spline2D(z, R, Bpol[1])
        return splR, splz

    def _f(self, ti):
        psii, f = self['psii'](ti).x[0] , self['f'](ti).x[0]
        cnd = (f != 0)
        psii, f = psii[cnd], f[cnd]
        perm = psii.argsort()
        return Spline(psii[perm], f[perm])

    def get_Bphi(self, ti, R, z=None, psi=None):
        if psi is None:
            psi = self.get_psi(ti, R, z)
        spl = self._f(ti)
        fill = spl.assignal().x[0]
        return spl.eval(psi, fill=fill) / R

    def get_Bphi_grid(self, ti, R=None, z=None, psi=None):
        R, z = self._check_grid(R, z)
        if psi is None:
            psi = self.get_psi_grid(ti, R, z)
        spl = self._f(ti)
        fill = spl.assignal().x[0]
        return spl.eval(psi, fill=fill) / R[None]

    def get_Bphi_spline(self, ti, R=None, z=None):
        R, z = self._check_grid(R, z)
        Bphi = self.get_Bphi_grid(ti, R, z)
        splp = Spline2D(z, R, Bphi)
        return splp

    def get_B_ratios_spline(self, ti, R=None, z=None):
        R, z = self._check_grid(R, z)
        Bpol = self.get_Bpol_grid(ti, R, z)
        Bphi = self.get_Bphi_grid(ti, R, z)
        splR = Spline2D(z, R, Bpol[0] / Bphi)
        splz = Spline2D(z, R, Bpol[1] / Bphi)
        return splR, splz

    def get_field_line_integrator(self, ti, bdry=None, bbox=None):
        splR, splz = self.get_B_ratios_spline(ti)
        if bdry is None:
            try:
                bdry = self.vessel.bdry
            except AttributeError:
                pass
        if bbox is None:
            try:
                bbox = self.vessel.get_bbox()
            except AttributeError:
                pass
        return FieldLineIntegrator(splR, splz, bdry, bbox)


class FieldLineViewer(ToggleViewerIntegrated):
    def __init__(self, eqi):
        self.eqi = eqi
        self.fli = None
        ToggleViewerIntegrated.__init__(self, menu_entry='Field lines')

    def set_fli(self, fli):
        # Update field line integrator when t_event has changed. Since this
        # means that the flux surfaces shown in the background have changed 
        # as well, invalidate cache of MouseMotionObserverViewer.
        self.fli = fli
        self.invalidate()

    def plotfun(self, event):
        fl = self.fli(event.xdata, event.ydata)
        print fl
        patch = fl.as_path_patch()
        self.ax.add_patch(patch)
        return [patch]


class FieldLineViewerVtk(ToggleViewerVtk):
    def __init__(self, eqi):
        self.eqi = eqi
        self.fli = None
        ToggleViewerVtk.__init__(self, menu_entry='VTK field lines')

    def set_fli(self, fli):
        # Update field line integrator when t_event has changed. Since this
        # means that the flux surfaces shown in the background have changed 
        # as well, invalidate cache of MouseMotionObserverViewer.
        self.fli = fli
        self.invalidate()

    def viewer(self, event=None):
        self.ax = self.eqi.vessel.render()

    def plotfun(self, event):
        fl = self.fli(event.xdata, event.ydata)
        print fl
        win = self.ax
        fl.render(win=win)
        return [win.ren.GetActors().GetLastActor()]


class EqiViewer(ToggleViewer):
    def __init__(self, eqi):
        self.eqi = eqi
        self.lvls = np.linspace(0., 1., 10)
        #self.lvls = np.linspace(eqi.psi.x.min(), eqi.psi.x.max(), 50)

        ToggleViewer.__init__(self, menu_entry='Eqi viewer')

    def viewer(self, event, **kw):
        self.flv = FieldLineViewer(self.eqi)
        self.flv_vtk = FieldLineViewerVtk(self.eqi)
        
        fig = get_tfig(menu_entries_ax=self.flv.menu_entries_ax + \
                                       self.flv_vtk.menu_entries_ax, **kw)
        self.ax = fig.axes[0]
        try:
            self.eqi.vessel.plot(self.ax)
        except AttributeError:
            pass

    def plotfun(self, event):
        t_event = event.xdata
        ax = self.ax
        FS = self.eqi.get_flux_surf(t_event, lvls=self.lvls)
        #FS = self.eqi.get_flux_surf_unnorm(t_event, lvls=self.lvls)
        FS.plot(ax)

        # assign field line integrator for t_event to FieldLineViewer
        fli = self.eqi.get_field_line_integrator(t_event)
        self.flv.set_fli(fli)
        self.flv_vtk.set_fli(fli)
        return ax.collections[-1:]


if __name__ == "__main__":
    from digitizer_aug import DigitizerAUGEQI
    from digitizer_d3d import DigitizerD3DEFIT

    ax = None
    lvls = np.linspace(0., 1., 10)

    dig_AUG = DigitizerAUGEQI(shn=30017)
    eqi_AUG = Eqi(dig_AUG)
    ax = eqi_AUG.get_flux_surf(3.45, lvls).plot(ax, edgecolors='b')
    ax = eqi_AUG.get_separatrix(3.45).plot(ax, edgecolors='b', linewidth=2)

    dig_D3D = DigitizerD3DEFIT(shn=141451)
    eqi_D3D = Eqi(dig_D3D)
    ax = eqi_D3D.get_flux_surf(1.65, lvls).plot(ax, edgecolors='r')
    ax = eqi_D3D.get_separatrix(1.65).plot(ax, edgecolors='r', linewidth=2)
    ax.axis('equal')

    show()


