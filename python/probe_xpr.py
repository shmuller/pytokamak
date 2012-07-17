import numpy as np

import scipy.interpolate as interp

import probe
reload(probe)

class IOMdsXPR(probe.IOMds):
    def __init__(self, *args):
        probe.IOMds.__init__(self, *args)
        self.mdsport = "8001"
        self.mdsfmt = 'augsignal(%d,"XPR","%s","AUGD")'


class IOFileXPR(probe.IOFile):
    def __init__(self, shn):
        probe.IOFile.__init__(self, shn=shn, suffix="_XPR", subdir="AUG")


class ProbeXPR(probe.Probe):
    def __init__(self, shn, sock=None):
        probe.Probe.__init__(self, shn, sock)

        self.IO_mds = IOMdsXPR(shn, sock)
        self.IO_file = IOFileXPR(shn)
        self.nodes = ('S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 't')

        self.Rcal = (442./5.6779, -106.6648)
        self.Vcal = (5.0952, -190.8193)
        self.I1cal = (1., 0.)
        self.I2cal = (0.5559, 0.)

    def mapsig(self):
        x = self.x

        t = x['t']
        R = x['S5']

        V = x['S1']
        I1 = x['S4']
        I2 = x['S2']

        M = 5000
        R[:] -= R[-M:].mean()
        I1[:] -= I1[:M].mean()
        I2[:] -= I2[:M].mean()

        R = probe.PositionSignal(t, R, 'R')
        V = probe.VoltageSignal(t, V, 'V')

        self.S = {'R': R,
                  'V': V,
                  'I1': probe.CurrentSignal(t, I1, 'I1', V),
                  'I2': probe.CurrentSignal(t, I2, 'I2', V)}

    def calib(self):
        R, V, I1, I2 = self['R', 'V', 'I1', 'I2']
        
        """
        R.x[:] = self.Rcal[0]*R.x + self.Rcal[1]
        V.x[:] = self.Vcal[0]*V.x + self.Vcal[1]

        I1.x[:] = self.I1cal[0]*I1.x + self.I1cal[1]
        I2.x[:] = self.I2cal[0]*I2.x + self.I2cal[1]
        """

if __name__ == "__main__":
    shn = 28435

    XPR = ProbeXPR(shn=shn)
    XPR.load()

    XPR.plot_raw()
    show()

