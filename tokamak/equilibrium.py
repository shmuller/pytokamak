import numpy as np

from scipy.ndimage import map_coordinates
from scipy.integrate import odeint

from LP.sig import memoized_property
from LP.splines import Spline, Spline2D

from sm_pyplot.tight_figure import get_tfig, get_axes, show
from sm_pyplot.observer_viewer import ToggleViewer
from sm_pyplot.vtk_contour import VtkContour


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


class FluxSurf:
    def __init__(self, vtk_ctr):
        self.vtk_ctr = vtk_ctr

    def plot(self, ax=None, **kw):
        kw.setdefault('edgecolors', 'r')
        ax = get_axes(ax)
        pc = self.vtk_ctr.as_path_collection(**kw)
        ax.add_collection(pc)
        return ax


class FieldLineIntegrator:
    def __init__(self, splR, splz):
        out = np.zeros(3)
        twopi = 2.*np.pi

        def dy_dt(y, t):
            R, z, l = y
            fact = twopi*R
            BR_Bphi = splR.eval(z, R)
            Bz_Bphi = splz.eval(z, R)
            out[0] = fact*BR_Bphi
            out[1] = fact*Bz_Bphi
            out[2] = fact*np.sqrt(1. + BR_Bphi**2 + Bz_Bphi**2)
            return out
        
        self.dy_dt = dy_dt

    def solve(self, y0, t):
        return odeint(self.dy_dt, y0, t)


class Eqi:
    def __init__(self, digitizer):
        self.digitizer = digitizer

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

    def get_flux_surf(self, ti, Lvls=None):
        x = self.R
        y = self.z
        z = np.zeros(1, np.float64)
        f = self.psi_n(ti).x[0]
        vtk_ctr = VtkContour(x, y, z, f, Lvls)
        vtk_ctr.contour()
        return FluxSurf(vtk_ctr)

    def get_separatrix(self, ti):
        return self.get_flux_surf(ti, Lvls=np.array([1.]))

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

    def get_B_ratios_spline(self, ti, R=None, z=None):
        R, z = self._check_grid(R, z)
        Bpol = self.get_Bpol_grid(ti, R, z)
        Bphi = self.get_Bphi_grid(ti, R, z)
        splR = Spline2D(z, R, Bpol[0] / Bphi)
        splz = Spline2D(z, R, Bpol[1] / Bphi)
        return splR, splz

    def get_field_line_integrator(self, ti):
        return FieldLineIntegrator(*self.get_B_ratios_spline(ti))


class EqiViewer(ToggleViewer):
    def __init__(self, eqi):
        self.eqi = eqi
        self.Lvls = np.linspace(0., 1., 10)

        ToggleViewer.__init__(self, 'Eqi viewer')

    def plotfun(self, event):
        t_event = event.xdata
        ax = self.ax
        self.clear()
        FS = self.eqi.get_flux_surf(t_event, Lvls=self.Lvls)
        FS.plot(ax)
        return ax.collections[-1:]
                

if __name__ == "__main__":
    from digitizer_aug import DigitizerAUGEQI
    from digitizer_d3d import DigitizerD3DEFIT

    Lvls = np.linspace(0., 1., 10)

    dig_AUG = DigitizerAUGEQI(shn=30017)
    eqi_AUG = Eqi(dig_AUG)
    ax = eqi_AUG.get_flux_surf(3.45, Lvls).plot()
    ax = eqi_AUG.get_separatrix(3.45).plot(ax, color='b', linewidth=2)

    dig_D3D = DigitizerD3DEFIT(shn=141451)
    eqi_D3D = Eqi(dig_D3D)
    ax = eqi_D3D.get_flux_surf(1.65, Lvls).plot(ax=ax)
    ax = eqi_D3D.get_separatrix(1.65).plot(ax, color='r', linewidth=2)
    ax.axis('equal')

    show()


