import numpy as np

from probe import IOMds, IOFile, Probe

class IOMdsAUG(IOMds):
    def __init__(self, *args):
        IOMds.__init__(self, *args)
        self.mdsport = "8001"
        self.mdsfmt = 'augsignal(%d,"LPS","%s","AUGD")'


class IOFileAUG(IOFile):
    def __init__(self, shn=27695):
        IOFile.__init__(self, shn=shn, subdir="AUG")


class ProbeAUG(Probe):
    def __init__(self, shn=27695, sock=None):
        Probe.__init__(self, shn, sock)

        self.IO_mds = IOMdsAUG(shn, sock)
        self.IO_file = IOFileAUG(shn)
        self.nodes = ('VOL3', 'VOL1', 'CUR1', 'CUR2', 't')

    def mapsig(self):
        s = slice(2048, None)
        x = self.x.view(s)

        self.t = x['t']
        self.R = x['VOL3']

        self.V = x['VOL1']
        self.I = np.c_[x['CUR1'], x['CUR2']]

        self.I_offs = np.mean(self.I[:100,:],0)

        self.V = -self.V
        self.I = self.I_offs - self.I



if __name__ == "__main__":
    shn = 27695

    XPR = ProbeAUG(shn=shn)
    XPR.load()

    XPR.plot()
    show()

