import numpy as np
import numpy.ma as ma

from sig import memoized_property, Digitizer

from sm_pyplot.tight_figure import get_fig, get_axes

from probe_xpr import TdiError, IOMdsAUG, IOFileAUG


class DigitizerMEM(Digitizer):
    def __init__(self, shn, sock=None, raw=False):
        Digitizer.__init__(self, shn, sock, name='MHC')

        self.IO_mds = IOMdsAUG(shn, sock, diag='MHC', raw=raw)
        self.IO_file = IOFileAUG(shn, suffix='_MEM')
        self.nodes = ('Usat_m09', 'Isat_m09', 'I_m07', 
                      'Ufl_m03', 'Ufl_m08', 'I_0', 'I_1', 't')


