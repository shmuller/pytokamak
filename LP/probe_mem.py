import numpy as np

from sig import PositionSignal

from tokamak.digitizer_aug import DigitizerAUG

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

    def _load_raw_factory(name):
        def load_raw(self, **kw):
            x = getattr(self.dig_lsm, name)()
            R = PositionSignal(-x['S-posi'].astype('d'), x['t'].astype('d'))

            i0, iM, i1 = R.t_ind
            t0, t1 = R.t[i0[0]] - 0.2, R.t[i1[-1]] + 0.2

            kw.setdefault('t0', t0)
            kw.setdefault('t1', t1)

            getattr(DigitizerMEM, name)(self, **kw)
            self.x['S-posi'] = R(self.x['t']).x.astype(np.float32)
            return self.x
        return load_raw

    load_raw      = _load_raw_factory('load_raw')
    load_raw_mds  = _load_raw_factory('load_raw_mds')
    load_raw_file = _load_raw_factory('load_raw_file')


