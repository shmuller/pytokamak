import numpy as np

import scipy.interpolate as interp

import probe
reload(probe)

class IOMdsLPS(probe.IOMds):
    def __init__(self, *args):
        probe.IOMds.__init__(self, *args)
        self.mdsport = "8001"
        self.mdsfmt = 'augsignal(%d,"LPS","%s","AUGD")'


class IOFileLPS(probe.IOFile):
    def __init__(self, shn=27695):
        probe.IOFile.__init__(self, shn=shn, suffix="_LPS", subdir="AUG")


class ProbeLPS(probe.Probe):
    def __init__(self, shn=27695, sock=None):
        probe.Probe.__init__(self, shn, sock)

        self.IO_mds = IOMdsLPS(shn, sock)
        self.IO_file = IOFileLPS(shn)
        self.nodes = ('VOL3', 'VOL1', 'CUR1', 'CUR2', 't')

        self.Rcal = (442./5.6779, -106.6648)
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

        M = 5000
        R[:] -= R[-M:].mean()
        V[:] = -V

        I1[:] = I1[:M].mean() - I1
        I2[:] = I2[:M].mean() - I2
        
        if self.shn < 27687:
            I2[:] = -I2
        
        It = I1 + I2

        R = probe.PositionSignal(t, R, name='R')
        V = probe.VoltageSignal(t, V, name='V')

        self.S = {'R': R,
                  'V': V,
                  'I1': probe.CurrentSignal(t, I1, V, name='I1'),
                  'I2': probe.CurrentSignal(t, I2, V, name='I2'),
                  'It': probe.CurrentSignal(t, It, V, name='It')}

    def calib(self):
        R, V, I1, I2, It = self['R', 'V', 'I1', 'I2', 'It']
        R.x[:] = self.Rcal[0]*R.x + self.Rcal[1]
        V.x[:] = self.Vcal[0]*V.x + self.Vcal[1]

        I1.x[:] = self.I1cal[0]*I1.x + self.I1cal[1]
        I2.x[:] = self.I2cal[0]*I2.x + self.I2cal[1]
        It.x[:] = I1.x + I2.x

    def current_calib(self):
        self.load(trim=True, calib=False)
        I1, I2 = self.S['I1'].x, self.S['I2'].x
        cnd1 = (I1 != I1.min()) & (I1 != I1.max())
        cnd2 = (I2 != I2.min()) & (I2 != I2.max())
        cnd = cnd1 & cnd2
        I1, I2 = I1[cnd], I2[cnd]

        fact = I1.dot(I2)/I2.dot(I2)
        return fact

    def position_calib(self):
        self.load(trim=True, calib=False)
        R = self['R']
        sp = interp.UnivariateSpline(R.t, R.x, s=10)
        Rs = sp(R.t)
        return Rs.max()


if __name__ == "__main__":
    shn = 27695

    LPS = ProbeLPS(shn=shn)
    LPS.load()

    LPS.plot_raw()
    show()

