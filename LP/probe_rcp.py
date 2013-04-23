import numpy as np

import probe
#reload(probe)

#import config
#reload(config)

IOMds = probe.IOMds
IOFile = probe.IOFile
Digitizer = probe.Digitizer
Amp = probe.Amp
Probe = probe.Probe

PositionSignal = probe.PositionSignal
VoltageSignal = probe.VoltageSignal
CurrentSignal = probe.CurrentSignal

ampUnity = Amp(fact=1., offs=0.)
ampInv   = Amp(fact=-1., offs=0.)
amp12Bit = Amp(fact=10./4095, offs=-5.)
amp14Bit = Amp(fact=20./16383, offs=-10.)


class IOMdsD3D(IOMds):
    def __init__(self, *args, **kw):
        diag = kw.pop('diag', '')
        suffix = kw.pop('suffix', '')

        IOMds.__init__(self, *args, **kw)
        self.mdsport = "8020"
        self.mdstree = diag
        self.mdsfmt = '\%s' + suffix

class IOFileD3D(IOFile):
    def __init__(self, *args, **kw):
        kw.setdefault('subdir', "D3D")
        IOFile.__init__(self, *args, **kw)


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
