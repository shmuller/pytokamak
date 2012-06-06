import numpy as np

from probe import IOMds, IOFile, Probe

class IOMdsAUG(IOMds):
    def __init__(self, *args):
        IOMds.__init__(self, *args)
        self.mdsport = "8001"
        self.mdsfmt = 'augsignal(%d,"LPS","%s","AUGD")'


class IOFileAUG(IOFile):
    def __init__(self, shn=27695):
        IOFile.__init__(self, "AUG", shn)


class ProbeAUG(Probe):
    def __init__(self, shn=27695, sock=None):
        Probe.__init__(self, shn, sock)

        self.IO_mds = IOMdsAUG(shn, sock)
        self.IO_file = IOFileAUG(shn)
        self.nodes = ('VOL3', 'VOL1', 'CUR1', 'CUR2', 'VOL2')

    def mapsig(self):
        self.R  = self.x['VOL3']
        self.Vb = self.x['VOL1']
        self.I1 = self.x['CUR1']
        self.I2 = self.x['CUR2']
        self.V  = self.x['VOL2']

        self.I1m = np.mean(self.I1[:100])
        self.I2m = np.mean(self.I2[:100])

        self.Vb = -self.Vb
        self.I1 = self.I1m-self.I1
        self.I2 = self.I2m-self.I2



if __name__ == "__main__":
    shn = 27695

    XPR = ProbeAUG(shn=shn)
    XPR.load()

    XPR.plot()
    show()

