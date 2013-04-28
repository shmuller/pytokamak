import numpy as np

from LP.sig import memoized_property

from sm_pyplot.tight_figure import get_tfig, get_axes, show
from sm_pyplot.observer_viewer import ToggleViewer
from sm_pyplot.vtk_contour import VtkContour

class FluxSurf:
    def __init__(self, vtk_ctr):
        self.vtk_ctr = vtk_ctr

    def plot(self, ax=None, **kw):
        ax = get_axes(ax)
        ax.add_collection(self.vtk_ctr.as_path_collection())
        return ax


class Eqi:
    def __init__(self, digitizer):
        self.digitizer = digitizer

    @memoized_property
    def R_z_psi(self):
        return self.digitizer.get_R_z_psi()

    @memoized_property
    def R(self):
        return self.R_z_psi[0]

    @memoized_property
    def z(self):
        return self.R_z_psi[1]

    @memoized_property
    def psi(self):
        return self.R_z_psi[2]

    @memoized_property
    def psi_n(self):
        psi0 = self.digitizer.get_R_z_psi_special('magnaxis')[2]
        psi1 = self.digitizer.get_R_z_psi_special('xpoint')[2]
        dpsi = psi1 - psi0
        return (self.psi - psi0[:, None, None]) / dpsi[:, None, None]

    def get_flux_surf(self, ti, Lvls=None):
        x = self.R.astype('d')
        y = self.z.astype('d')
        z = np.zeros(1, 'd')
        f = self.psi_n(ti).x[0].astype('d')
        vtk_ctr = VtkContour(x, y, z, f, Lvls)
        vtk_ctr.contour()
        return FluxSurf(vtk_ctr)

    def get_separatrix(self, ti):
        return self.get_flux_surf(ti, Lvls=np.array([1.]))


class EqiViewer(ToggleViewer):
    def __init__(self, eqi):
        self.eqi = eqi
        self.Lvls = np.linspace(0., 1., 10)

        ToggleViewer.__init__(self, 'Eqi viewer')

    def plotfun(self, event):
        t_event = event.xdata
        FS = self.eqi.get_flux_surf(t_event, Lvls=self.Lvls)
        ax = self.ax
        ax.collections = []
        FS.plot(ax)
        return ax.collections

    def viewer(self, event):
        fig = get_tfig(figsize=(4,4), xlab="R (m)", ylab="z (m)")
        self.ax = fig.axes[0]
        self.ax.set_xlim((0., 2.5))
        self.ax.set_ylim((-1.5, 1.))


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


