import math
import numpy as np
from utils.utils import memoized_property

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

