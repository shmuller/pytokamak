import numpy as np

import scipy.interpolate as interp

if __name__ == "__main__":
    import matplotlib
    #matplotlib.use('TkAgg')
    matplotlib.use('Qt4Agg')
    import matplotlib.pyplot as plt

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

        self.Rcal = ((170+72)/3.3630, -72)
        self.Vcal = (100., -69.5700)
        self.I1cal = (0.5*20./10, 0.)  # mA/mV = A/V (0.5 from missing 50 Ohm term.)
        self.I2cal = (0.5*20./10, 0.)

    def mapsig(self):
        x = self.x

        t = x['t'].astype('d')
        R = x['S5'].astype('d')

        V = x['S1'].astype('d')
        I1 = x['S4'].astype('d')
        I2 = x['S2'].astype('d')

        M = 5000
        R[:] -= R[-M:].mean()
        I1[:] -= I1[:M].mean()
        I2[:] -= I2[:M].mean()

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

    def position_calib(self, plunge=0):
        self.load(plunge=plunge, trim=True, calib=False)
        R = self['R']
        sp = interp.UnivariateSpline(R.t, R.x, s=10)
        Rs = sp(R.t)
        return Rs

    def voltage_calib(self):
        self.load(trim=False, calib=False)
        return self['V'].x.mean()


if __name__ == "__main__":
    shn = 28469

    XPR = ProbeXPR(shn=shn)
    XPR.analyze(plunge=0)

    XPR.plot()
    plt.show()



