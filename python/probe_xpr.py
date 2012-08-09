import numpy as np

import scipy.interpolate as interp

if __name__ == "__main__":
    import matplotlib
    #matplotlib.use('TkAgg')
    matplotlib.use('Qt4Agg')
    import matplotlib.pyplot as plt

import probe
#reload(probe)


class IOMdsXPR(probe.IOMds):
    def __init__(self, *args):
        # augsignal(_shot, _diag, _signame, _experiment, _edition, 
        #   _t1, _t2, _oshot, _oedition, _qual)

        probe.IOMds.__init__(self, *args)
        self.mdsport = "8001"
        self.mdsfmt = '_s = augsignal(%d, "XPR", "%s", "AUGD", *, *, *, *, *, "raw")'

        self.datadeco = '%s; word_unsigned(data(_s))'
        self.timedeco = '%s; dim_of(_s)'
        self.sizedeco = '%s; size(_s)'


class IOFileXPR(probe.IOFile):
    def __init__(self, shn):
        probe.IOFile.__init__(self, shn=shn, suffix="_XPR", subdir="AUG")


class ProbeXPR(probe.Probe):
    def __init__(self, shn, sock=None):
        probe.Probe.__init__(self, shn, sock)

        self.IO_mds = IOMdsXPR(shn, sock)
        self.IO_file = IOFileXPR(shn)
        self.nodes = ('S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 't')

    def mapsig(self):
        b2V = probe.Amp(fact=20./16383, offs=-10.)
        x = self.x

        t  = x['t']
        R  = b2V*x['S5'].astype('d')
        V  = b2V*x['S1'].astype('d')
        I1 = b2V*x['S4'].astype('d')
        I2 = b2V*x['S2'].astype('d')
        VF = b2V*x['S6'].astype('d')

        R = probe.PositionSignal(t, R, name='R')
        V = probe.VoltageSignal(t, V, name='V')

        self.S = {'R': R,
                  'V': V,
                  'I1': probe.CurrentSignal(t, I1, V, name='I1'),
                  'I2': probe.CurrentSignal(t, I2, V, name='I2'),
                  'VF': probe.VoltageSignal(t, VF, name='VF')}

    def calib(self):
        #Rcal = ((170+72)/3.3630, -72)
        #Vcal = (100., -69.2227)
        #I1cal = (0.5*20./10, 0.)  
        #I2cal = (0.5*20./10, 0.)
        #VFcal = (100., 0.)

        fixpoints = (3.6812, -72), (7.0382, 170)
        ampR = probe.Amp(fixpoints=fixpoints)
        ampV = probe.Amp(fact=100., offs=-69.2227)

        ampI1 = probe.Amp(fact=0.5*20./10)  # mA/mV = A/V (0.5 from missing 50 Ohm term.)
        ampI2 = probe.Amp(fact=0.5*20./10)
        ampVF = probe.Amp(fact=100.)

        R, V, I1, I2, VF = self['R', 'V', 'I1', 'I2', 'VF']

        M = 5000
        #R.x[:] -= R.x[-M:].mean()
        #R.x[:] = Rcal[0]*R.x + Rcal[1]

        R *= ampR
        V *= ampV

        I1 *= ampI1
        I2 *= ampI2
        VF *= ampVF

        I1.x[:] -= I1.x[-M:].mean()
        I2.x[:] -= I2.x[-M:].mean()
        VF.x[:] -= VF.x[-M:].mean()

        self.S['It'] = probe.CurrentSignal(I1.t, I1.x + I2.x, V, name='It')

    def position_calib(self):
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



