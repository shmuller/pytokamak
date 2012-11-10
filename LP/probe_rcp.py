import copy
import numpy as np

import scipy.interpolate as interp

if __name__ == "__main__":
    import matplotlib
    #matplotlib.use('TkAgg')
    matplotlib.use('Qt4Agg')
    import matplotlib.pyplot as plt

import probe
#reload(probe)

#import config
#reload(config)

IOMds = probe.IOMds
IOFile = probe.IOFile
Digitizer = probe.Digitizer
Amp = probe.Amp
Probe = probe.Probe

PositionSignal = probe.PositionSignal
VoltageSignal = probe.VoltageSignal
CurrentSignal = probe.CurrentSignal

ampUnity = Amp(fact=1., offs=0.)
ampInv   = Amp(fact=-1., offs=0.)
amp12Bit = Amp(fact=10./4095, offs=-5.)
amp14Bit = Amp(fact=20./16383, offs=-10.)


class IOMdsD3D(IOMds):
    def __init__(self, *args, **kw):
        diag = kw.pop('diag', 'RCP')
        plunge = kw.pop('plunge', 1)

        IOMds.__init__(self, *args, **kw)
        self.mdsport = "8020"
        self.mdstree = diag
        self.mdsfmt = '\%%s_%d' % plunge

class IOFileD3D(IOFile):
    def __init__(self, shn, diag='RPC', plunge=1):
        IOFile.__init__(self, shn=shn, suffix="_%s_%d" % (diag, plunge), subdir="D3D")


class DigitizerRCP(Digitizer):
    def __init__(self, shn, sock=None, plunge=1):
        Digitizer.__init__(self, shn, sock, name='RCP')

        self.IO_mds = IOMdsD3D(shn, sock, diag='RCP', plunge=plunge)
        self.IO_file = IOFileD3D(shn, diag='RCP', plunge=plunge)

        self.nodes = ('FPCALPOS', 'ISATV', 'ISATM1', 'ISATM2', 
            'VF1', 'VF2', 'TE1V', 'TE2I', 't')
        
        self.amp = {node: ampUnity.copy() for node in self.nodes}
        self.amp['t'] = ampUnity.copy()
        
    def calib(self):
        Digitizer.calib(self)


class ProbeRCP(Probe):
    def __init__(self, shn, sock=None, plunge=1):
        
        digitizer = DigitizerRCP(shn, sock, plunge)
        Probe.__init__(self, digitizer)

    def mapsig(self):
        """
        mapping = self.config.mapping[self.digitizer.name]

        x = self.x[self.digitizer.window]

        t = x['t']
        R = x[mapping.R].astype('d')
        self.S = dict(R=PositionSignal(R, t, name='R'))

        unique_sigs = {k: x[k].astype('d') for k in mapping.unique('VI')}

        for i, (mapV, mapI) in enumerate(mapping.VI, start=1):
            if mapV is None:
                V = None
            else:
                V = VoltageSignal(unique_sigs[mapV], t, name='V%d' % i)
            if mapI is None:
                self.S[i] = V
            else:
                self.S[i] = CurrentSignal(unique_sigs[mapI], t, V, name='I%d' % i)
        """

    def calib(self):
        """
        amp = self.config.amp[self.digitizer.name]
        
        s = slice(-5000, None)
        for i, (ampV, ampI) in enumerate(amp.VI, start=1):
            S = self.S[i]
            if isinstance(S, CurrentSignal):
                S *= ampI
                V = S.V
            else:
                V = S
            if V is not None and ampV is not None:
                V *= ampV
            S.norm_to_region(s)

        self.S['R'] *= amp.R
        self.S['V'] = self.S[1].V
        self.S['It'] = self.S[1] + self.S[2]
        """

    def position_calib(self):
        """
        R = self['R']
        sp = interp.UnivariateSpline(R.t, R.x, s=10)
        Rs = sp(R.t)
        return Rs
        """

    def voltage_calib(self):
        """
        self.load(trim=False, calib=False)
        return self['V'].x.mean()
        """
   

if __name__ == "__main__":
    shn = 141451
    plunge = 2

    RCP = ProbeRCP(shn=shn, plunge=plunge)

    #XPR.plot()
    #plt.show()



