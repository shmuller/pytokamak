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

        self.t = x['t']
        self.R = x['VOL3']

        self.V = x['VOL1']
        self.I = np.r_[x['CUR1'], x['CUR2']].reshape((2,-1))

        self.I_offs = np.mean(self.I[:,:100],1)

        self.R -= self.R[0]
        self.V = -self.V
        self.I = self.I_offs[:,None] - self.I

        R = probe.PositionSignal(self.t, self.R, 'R')
        V = probe.VoltageSignal(self.t, self.V, 'V')

        self.S = {'R': R,
                  'V': V,
                  'I1': probe.CurrentSignal(self.t, self.I[0], 'I1', V),
                  'I2': probe.CurrentSignal(self.t, self.I[1], 'I2', V)}



if __name__ == "__main__":
    shn = 27695

    XPR = ProbeAUG(shn=shn)
    XPR.load()

    XPR.plot()
    show()

