import numpy as np

from probe import PositionSignal

from pytokamak.tokamak.digitizer_aug import DigitizerAUG

class DigitizerMEM(DigitizerAUG):
    def __init__(self, shn, raw=False, s=slice(None, None, 4), **kw):
        DigitizerAUG.__init__(self, shn, diag='MHC', suffix='_MEM', raw=raw, s=s,
            nodes=('Isat_m05', 'Isat_m10'), **kw)

        #self.nodes = ('Usat_m09', 'Isat_m09', 'I_m07', 
        #              'Ufl_m03', 'Ufl_m08', 'I_0', 'I_1', 't')

        #self.nodes = ('Usat_m05', 'Isat_m05', 'Usat_m10', 'Isat_m10', 
        #              'Ufl_m02', 'Ufl_m04', 'U_m03', 'I_m03', 't')


class DigitizerMEMPos(DigitizerMEM):
    def __init__(self, shn, **kw):
        DigitizerMEM.__init__(self, shn, **kw)

        self.dig_lsm = DigitizerAUG(shn, diag='LSM', suffix='_MEM', nodes=('S-posi',))

    def load_raw_mds(self, **kw):
        R = -self.dig_lsm['S-posi'].astype(np.float64)
        i0, iM, i1 = PositionSignal(R.x, R.t).t_ind
        t0, t1 = R.t[i0[0]] - 0.2, R.t[i1[-1]] + 0.2

        kw.setdefault('t0', t0)
        kw.setdefault('t1', t1)
        return DigitizerMEM.load_raw_mds(self, **kw)

    def calib(self):
        DigitizerMEM.calib(self)
        R = self.dig_lsm['S-posi'].astype(np.float64)
        self.x['XPOS'] = R(self.x['t']).x.astype(np.float32)

