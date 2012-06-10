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

    def mapsig(self):
        s = slice(2048, None)
        x = self.x.view(s)

        t = x['t']
        R = x['VOL3']

        V = x['VOL1']
        I1 = x['CUR1']
        I2 = x['CUR2']

        R[:] -= R[0]
        V[:] = -V

        I1[:] = I1[:100].mean() - I1
        I2[:] = I2[:100].mean() - I2

        R = probe.PositionSignal(t, R, 'R')
        V = probe.VoltageSignal(t, V, 'V')

        self.S = {'R': R,
                  'V': V,
                  'I1': probe.CurrentSignal(t, I1, 'I1', V),
                  'I2': probe.CurrentSignal(t, I2, 'I2', V)}



if __name__ == "__main__":
    shn = 27695

    XPR = ProbeAUG(shn=shn)
    XPR.load()

    XPR.plot()
    show()

