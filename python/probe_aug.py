import numpy as np

import probe
reload(probe)

class IOMdsAUG(probe.IOMds):
    def __init__(self, *args):
        probe.IOMds.__init__(self, *args)
        self.mdsport = "8001"
        self.mdsfmt = 'augsignal(%d,"LPS","%s","AUGD")'


class IOFileAUG(probe.IOFile):
    def __init__(self, shn=27695):
        probe.IOFile.__init__(self, shn=shn, subdir="AUG")


class ProbeAUG(probe.Probe):
    def __init__(self, shn=27695, sock=None):
        probe.Probe.__init__(self, shn, sock)

        self.IO_mds = IOMdsAUG(shn, sock)
        self.IO_file = IOFileAUG(shn)
        self.nodes = ('VOL3', 'VOL1', 'CUR1', 'CUR2', 't')

        self.Vcal = (5.0952, -190.8193)
        self.I1cal = (1., 0.)
        self.I2cal = (0.5559, 0.)

    def mapsig(self):
        s = slice(2048, None)
        x = self.x.view(s)

        t = x['t']
        R = x['VOL3']

        V = x['VOL1']
        I1 = x['CUR1']
        I2 = x['CUR2']

        R[:] -= R[0]
        V[:] = self.Vcal[0]*(-V) + self.Vcal[1]

        I1[:] = I1[:1000].mean() - I1
        I2[:] = I2[:1000].mean() - I2

        if self.shn < 27687:
            I2[:] = -I2

        I1[:] = self.I1cal[0]*I1 + self.I1cal[1]
        I2[:] = self.I2cal[0]*I2 + self.I2cal[1]

        R = probe.PositionSignal(t, R, 'R')
        V = probe.VoltageSignal(t, V, 'V')

        self.S = {'R': R,
                  'V': V,
                  'I1': probe.CurrentSignal(t, I1, 'I1', V),
                  'I2': probe.CurrentSignal(t, I2, 'I2', V)}

    def current_calib(self):
        I1, I2 = self.S['I1'].x, self.S['I2'].x
        cnd1 = (I1 != I1.min()) & (I1 != I1.max())
        cnd2 = (I2 != I2.min()) & (I2 != I2.max())
        cnd = cnd1 & cnd2
        I1, I2 = I1[cnd], I2[cnd]

        fact = I1.dot(I2)/I2.dot(I2)
        return fact


if __name__ == "__main__":
    shn = 27695

    XPR = ProbeAUG(shn=shn)
    XPR.load()

    XPR.plot()
    show()

