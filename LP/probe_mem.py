import numpy as np
import numpy.ma as ma

from sig import memoized_property, Digitizer, PositionSignal

from sm_pyplot.tight_figure import get_fig, get_axes

from probe_xpr import TdiError, IOMdsAUG, IOFileAUG

class DigitizerMEMPos(Digitizer):
    def __init__(self, shn, sock=None):
        Digitizer.__init__(self, shn, sock, name='MEM_POS')

        self.IO_mds = IOMdsAUG(shn, sock, diag='LSM')
        self.IO_file = IOFileAUG(shn, suffix='_MEM_POS')
        self.nodes = ('S-posi', 't')


class DigitizerMEM(Digitizer):
    def __init__(self, shn, sock=None, raw=False):
        Digitizer.__init__(self, shn, sock, name='MHC')

        self.IO_mds = IOMdsAUG(shn, sock, diag='MHC', raw=raw)
        self.IO_file = IOFileAUG(shn, suffix='_MEM')
        self.nodes = ('Usat_m09', 'Isat_m09', 'I_m07', 
                      'Ufl_m03', 'Ufl_m08', 'I_0', 'I_1', 't')


class DigitizerMEMCombi(Digitizer):
    def __init__(self, shn, sock=None, raw=False):
        Digitizer.__init__(self, shn, sock, name='MHC')

        self.IO_mds = IOMdsAUG(shn, sock, diag='MHC', raw=raw)
        self.IO_file = IOFileAUG(shn, suffix='_MEM_combi')
        self.nodes = ('Usat_m05', 'Isat_m05', 'Usat_m10', 'Isat_m10', 
                      'Ufl_m02', 'Ufl_m04', 'U_m03', 'I_m03', 't')


