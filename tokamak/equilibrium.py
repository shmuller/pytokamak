import numpy as np

from LP.sig import memoized_property, Digitizer
from LP.probe_xpr import IOMdsAUG, IOFileAUG

from sm_pyplot.tight_figure import get_axes, show
from sm_pyplot.vtk_contour import VtkContour


class DigitizerEQI(Digitizer):
    def __init__(self, shn, sock=None):
        Digitizer.__init__(self, shn, sock, name='EQI')
        self.tnode, self.tunits = 'time', 's'

        self.IO_mds = IOMdsAUG(shn, sock, diag='EQI')
        self.IO_file = IOFileAUG(shn, suffix='_EQI')
        self.nodes = ('time', 'Ri', 'Zj', 'PFM', 'ikCAT', 'RPFx', 'zPFx', 'PFxx')

        self.mapspec = dict(magnaxis=0, xpoint=1, innerlim=2, xpoint2=3, outerlim=4)

    def __getitem__(self, indx):
        t = self.x['time']
        x = self.x[indx][:t.shape[0]]
        if indx == 'PFM':
            x = x[:, :self.x['Zj'].shape[1], :self.x['Ri'].shape[1]]
        return self.assignal(x, t, name=indx)

    def get_R_z_psi(self):
        return self['Ri'], self['Zj'], self['PFM']

    def get_R_z_psi_special(self, spec):
        i = self.mapspec[spec]
        return self['RPFx'][:,i], self['zPFx'][:,i], self['PFxx'][:,i]


class FluxSurf:
    def __init__(self, vtk_ctr):
        self.vtk_ctr = vtk_ctr

    def plot(self, ax=None):
        ax = get_axes(ax)
        self.vtk_ctr.plot(ax=ax)
        return ax


class EQI:
    def __init__(self, shn):
        self.digitizer = None

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
        x = self.R(ti).x[0].astype('d')
        y = self.z(ti).x[0].astype('d')
        z = np.zeros(1, 'd')
        f = self.psi_n(ti).x[0].astype('d')
        vtk_ctr = VtkContour(x, y, z, f, Lvls)
        vtk_ctr.contour()
        return FluxSurf(vtk_ctr)


class EQIAUG(EQI):
    def __init__(self, shn):
        self.digitizer = DigitizerEQI(shn=shn)


if __name__ == "__main__":
    eqi = EQIAUG(shn=30017)
    FS = eqi.get_flux_surf(3.45)
    ax = FS.plot()
    show()



