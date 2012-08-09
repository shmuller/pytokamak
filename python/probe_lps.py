import numpy as np

import scipy.interpolate as interp

import probe
#reload(probe)

class IOMdsLPS(probe.IOMds):
    def __init__(self, *args):
        # augsignal(_shot, _diag, _signame, _experiment, _edition, 
        #   _t1, _t2, _oshot, _oedition, _qual)

        probe.IOMds.__init__(self, *args)
        self.mdsport = "8001"
        self.mdsfmt = '_s = augsignal(%d, "LPS", "%s", "AUGD", *, *, *, *, *, "raw")'

        self.datadeco = '%s; word_unsigned(data(_s))'
        self.timedeco = '%s; dim_of(_s)'
        self.sizedeco = '%s; size(_s)'


class IOFileLPS(probe.IOFile):
    def __init__(self, shn=27695):
        probe.IOFile.__init__(self, shn=shn, suffix="_LPS", subdir="AUG")


class ProbeLPS(probe.Probe):
    def __init__(self, shn=27695, sock=None):
        probe.Probe.__init__(self, shn, sock)

        self.IO_mds = IOMdsLPS(shn, sock)
        self.IO_file = IOFileLPS(shn)
        self.nodes = ('CUR1', 'VOL1', 'CUR2', 'VOL2', 'VOL3', 'VOL4', 't')

    def mapsig(self):
        s = slice(2048, None)
        x = self.x
        t = x['t'][s]
        
        b2V = probe.Amp(fact=10./4095, offs=-5.)
        inv = probe.Amp(fact=-1.)

        ampR  = inv*b2V
        ampV  = b2V
        ampI1 = inv*b2V
        ampI2 = b2V
        ampVF = b2V

        R  = ampR*x['VOL3'][s].astype('d')
        V  = ampV*x['VOL1'][s].astype('d')
        I1 = ampI1*x['CUR1'][s].astype('d')
        I2 = ampI2*x['CUR2'][s].astype('d')
        VF = ampVF*x['VOL2'][s].astype('d')
        
        R = probe.PositionSignal(t, R, name='R')
        V = probe.VoltageSignal(t, V, name='V')

        self.S = {'R': R,
                  'V': V,
                  'I1': probe.CurrentSignal(t, I1, V, name='I1'),
                  'I2': probe.CurrentSignal(t, I2, V, name='I2'), 
                  'VF': probe.VoltageSignal(t, VF, name='VF')}

    def calib(self):
        fixpoints = (-1.8767, -106), (3.8011, 336)
        ampR = probe.Amp(fixpoints=fixpoints)
        ampV = probe.Amp(fact=100., offs=-183.76)

        ampI1 = probe.Amp(fact=0.5*20./10)  # mA/mV = A/V (0.5 from missing 50 Ohm term.)
        ampI2 = probe.Amp(fact=0.5*20./10)
        ampVF = probe.Amp(fact=100.)

        inv = probe.Amp(fact=-1.)
        amp1x5 = probe.Amp(fact=2.58)
        amp2x5 = probe.Amp(fact=4.84)

        ampI1 *= amp1x5.inv()
        ampI2 *= amp2x5.inv()

        if self.shn < 28426:
            ampI1 *= inv

        if self.shn < 27687:
            ampI2 *= inv
    
        R, V, I1, I2, VF = self['R', 'V', 'I1', 'I2', 'VF']

        R  *= ampR
        V  *= ampV
        I1 *= ampI1
        I2 *= ampI2
        VF *= ampVF

        s = slice(-5000, None)
        I1.norm_to_region(s)
        I2.norm_to_region(s)
        VF.norm_to_region(s)

        self.S['It'] = probe.CurrentSignal(I1.t, I1.x + I2.x, V, name='It')
        
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

