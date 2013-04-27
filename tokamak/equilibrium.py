import numpy as np

from LP.sig import memoized_property, Digitizer, Amp
from LP.probe_xpr import IOMdsAUG, IOFileAUG
from LP.probe_rcp import IOMdsD3D, IOFileD3D

from sm_pyplot.tight_figure import get_axes, show
from sm_pyplot.vtk_contour import VtkContour


class DigitizerEQIAUG(Digitizer):
    def __init__(self, shn):
        Digitizer.__init__(self, shn, name='EQI')
        self.IO_mds = IOMdsAUG(shn, diag='EQI')
        self.IO_file = IOFileAUG(shn, suffix='_AUG', group='EQI')
        self.tnode, self.tunits = 'time', 's'
        self.nodes = ('time', 'Ri', 'Zj', 'PFM', 'ikCAT', 'RPFx', 'zPFx', 'PFxx')
        self.mapspec = dict(magnaxis=0, xpoint=1, innerlim=2, xpoint2=3, outerlim=4)
   
    def get_R_z_psi(self):
        return self.x['Ri'][0], self.x['Zj'][0], self['PFM']

    def get_R_z_psi_special(self, spec):
        i = self.mapspec[spec]
        return self['RPFx'][:,i], self['zPFx'][:,i], self['PFxx'][:,i]

    def calib(self):
        for node, x in self.x.iteritems():
            if x.ndim > 1:
                self.x[node] = x.transpose(np.roll(np.arange(x.ndim), 1))


class DigitizerEFITD3D(Digitizer):
    def __init__(self, shn):
        Digitizer.__init__(self, shn, name='EFIT')
        self.IO_mds = IOMdsD3D(shn, diag='EFIT01')
        self.IO_file = IOFileD3D(shn, suffix='_D3D', group='EFIT01')
        self.tnode, self.tunits = 'GTIME', 's'
        self.nodes = ('GTIME', 'R', 'Z', 'PSIRZ', 
                      'RMAXIS', 'ZMAXIS', 'SSIMAG', 'SSIBRY', 'BDRY')

        self.amp = dict(GTIME=Amp(1e-3))

    def get_R_z_psi(self):
        return self.x['R'], self.x['Z'], self['PSIRZ']

    def get_R_z_psi_special(self, spec):
        if spec == 'magnaxis':
            return self['RMAXIS'], self['ZMAXIS'], self['SSIMAG']
        elif spec == 'xpoint':
            return None, None, self['SSIBRY']


class FluxSurf:
    def __init__(self, vtk_ctr):
        self.vtk_ctr = vtk_ctr

    def plot(self, ax=None, **kw):
        ax = get_axes(ax)
        self.vtk_ctr.plot(ax=ax, **kw)
        return ax


class EQI:
    def __init__(self, shn, DigitizerClass):
        self.digitizer = DigitizerClass(shn=shn)

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


if __name__ == "__main__":
    Lvls = np.linspace(0., 1., 10)

    eqi_AUG = EQI(shn=30017, DigitizerClass=DigitizerEQIAUG)
    ax = eqi_AUG.get_flux_surf(3.45, Lvls).plot()
    ax = eqi_AUG.get_separatrix(3.45).plot(ax, color='b', linewidth=2)

    eqi_D3D = EQI(shn=141451, DigitizerClass=DigitizerEFITD3D)
    ax = eqi_D3D.get_flux_surf(1.65, Lvls).plot(ax=ax)
    ax = eqi_D3D.get_separatrix(1.65).plot(ax, color='r', linewidth=2)
    ax.axis('equal')

    show()


