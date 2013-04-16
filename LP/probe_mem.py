import numpy as np
import numpy.ma as ma

from sig import memoized_property, Digitizer, PositionSignal

from sm_pyplot.tight_figure import get_fig, get_axes

from probe_xpr import TdiError, IOMdsAUG, IOFileAUG

class DigitizerLSM(Digitizer):
    def __init__(self, shn, sock=None):
        Digitizer.__init__(self, shn, sock, name='LSM')

        self.IO_mds = IOMdsAUG(shn, sock, diag='LSM')
        self.IO_file = IOFileAUG(shn, suffix='_MEM_POS')
        self.nodes = ('S-posi', 't')


class DigitizerMEM(Digitizer):
    def __init__(self, shn, sock=None, raw=False):
        Digitizer.__init__(self, shn, sock, name='MHC')

        self.IO_mds = IOMdsAUG(shn, sock, diag='MHC', raw=raw)
        self.IO_file = IOFileAUG(shn, suffix='_MEM')
        
        self.nodes = ('Isat_m05', 'Isat_m10', 't')

        #self.nodes = ('Usat_m09', 'Isat_m09', 'I_m07', 
        #              'Ufl_m03', 'Ufl_m08', 'I_0', 'I_1', 't')


class DigitizerMEMPos(DigitizerMEM):
    def __init__(self, shn, sock=None, raw=False):
        DigitizerMEM.__init__(self, shn, sock)

        self.dig_lsm = DigitizerLSM(self.shn, self.sock)

        self.more_nodes = ('S-posi',)

    def _load_raw_factory(name):
        def load_raw(self, **kw):
            x = getattr(self.dig_lsm, name)()
            R = PositionSignal(-x['S-posi'].astype('d'), x['t'].astype('d'))

            i0, iM, i1 = R.t_ind
            t0, t1 = R.t[i0[0]] - 0.2, R.t[i1[-1]] + 0.2

            kw.setdefault('t0', t0)
            kw.setdefault('t1', t1)

            getattr(DigitizerMEM, name)(self, **kw)
            self.x['S-posi'] = R(self.x['t'])
            return self.x
        return load_raw

    load_raw      = _load_raw_factory('load_raw')
    load_raw_mds  = _load_raw_factory('load_raw_mds')
    load_raw_file = _load_raw_factory('load_raw_file')

    def save(self):
        DigitizerMEM.save(self)
        self.dig_lsm.save()


class DigitizerMEMCombi(Digitizer):
    def __init__(self, shn, sock=None, raw=False):
        Digitizer.__init__(self, shn, sock, name='MHC')

        self.IO_mds = IOMdsAUG(shn, sock, diag='MHC', raw=raw)
        self.IO_file = IOFileAUG(shn, suffix='_MEM_combi')
        self.nodes = ('Usat_m05', 'Isat_m05', 'Usat_m10', 'Isat_m10', 
                      'Ufl_m02', 'Ufl_m04', 'U_m03', 'I_m03', 't')


