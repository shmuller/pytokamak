import numpy as np

from pytokamak.tokamak.digitizer import Digitizer
from pytokamak.tokamak.digitizer_d3d import IOMdsD3D, IOFileD3D

from probe import Probe

class DigitizerRCP(Digitizer):
    def __init__(self, shn, plunge=1):
        Digitizer.__init__(self, shn, name='RCP')

        self.IO_mds = IOMdsD3D(shn, diag='RCP', suffix="_%d" % plunge)
        self.IO_file = IOFileD3D(shn, suffix="_RCP", group='plunge_%d' % plunge)

        self.nodes = ('FPCALPOS', 'ISATV', 'ISATM1', 'ISATM2', 
            'VF1', 'VF2', 'TE1V', 'TE2I', 't')
        

class ProbeRCP(Probe):
    def __init__(self, shn, plunge=1):    
        digitizer = DigitizerRCP(shn, plunge)
        Probe.__init__(self, digitizer)


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    shn = 141451
    plunge = 2

    RCP = ProbeRCP(shn=shn, plunge=plunge)

    fig = RCP.digitizer.plot()

    plt.show()
