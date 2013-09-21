import math
import numpy as np
from utils.utils import memoized_property
from utils.sig import Signal2D
from vtk_aug import VtkContour

class DCN:
    def __init__(self, shn, eqi=None):
        from digitizer_aug import DigitizerAUGDCR
        from diaggeom_aug import DCNGeom

        self.digitizer = DigitizerAUGDCR(shn)
        self.geom = DCNGeom()
        self.eqi = eqi
   
    def __getitem__(self, indx):
        return self.digitizer[indx].t_gt(0.).nonneg()

    @memoized_property
    def _psi_spl(self):
        (Rm, zm), (RM, zM) = self.geom.get_rays('H-5')
        
        ni = 100
        si = np.linspace(0., 1., ni)
        Ri = np.linspace(Rm, RM, ni)
        zi = np.linspace(zm, zM, ni)
        return self.eqi.get_path_spline(si, Ri, zi)

    @memoized_property
    def psi(self):
        S = self['H-5'].compressed()
        t = S.t
        s = np.linspace(0., 1., 100)
        return Signal2D(t, s, self._psi_spl(t, s).T, type='Normalized flux',
                        xtype='t', xunits=S.tunits, ytype='Ray coordinate', yunits='')

    def get_separatrix(self):
        t = np.ascontiguousarray(self.psi.x, np.float64)
        s = np.ascontiguousarray(self.psi.y, np.float64)
        psi = np.ascontiguousarray(self.psi.Z.base, np.float64)
        return VtkContour(s, t, psi, lvls=1.)

