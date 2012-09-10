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

import config
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


class IOMdsAUG(IOMds):
    def __init__(self, *args, **kw):
        # augsignal(_shot, _diag, _signame, _experiment, _edition, 
        #   _t1, _t2, _oshot, _oedition, _qual)

        diag = kw.pop('diag', 'XPR')

        IOMds.__init__(self, *args, **kw)
        self.mdsport = "8001"
        self.mdsfmt = '_s = augsignal(%d,"%s","%%s","AUGD",*,*,*,*,*,"raw")' % (self.shn, diag)

        self.datadeco = '%s; word_unsigned(data(_s))'
        self.timedeco = '%s; dim_of(_s)'
        self.sizedeco = '%s; size(_s)'


class IOFileAUG(IOFile):
    def __init__(self, shn, diag='XPR'):
        IOFile.__init__(self, shn=shn, suffix="_"+diag, subdir="AUG")


class DigitizerXPR(Digitizer):
    def __init__(self, shn, sock=None):
        Digitizer.__init__(self, shn, sock, name='XPR')

        self.IO_mds = IOMdsAUG(shn, sock, diag='XPR')
        self.IO_file = IOFileAUG(shn, diag='XPR')
        self.nodes = ('S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 't')

        #f = 20./16283
        #offs = (8759, 8830, 8756, 8770, 8759, 8747, 8746, 8752, 0)
        #self.amp = {node: Amp(fact=f, offs=-f*o) for node, o in zip(self.nodes, offs)}
        self.amp = {node: amp14Bit.copy() for node in self.nodes}
        self.amp['t'] = ampUnity.copy()

    def update_window(self):
        """Detect failing data readout and cut window
        """
        t = self.x['t']
        dt = np.abs(np.diff(t))
        cnd = dt > 2*dt[0]
        iM = np.argmax(cnd)

        stop = None
        if cnd[iM] == True:
            stop = iM+1
        self.window = slice(None, stop)
        
    def calib(self):
        self.update_window()
        Digitizer.calib(self)


class DigitizerLPS(Digitizer):
    def __init__(self, shn, sock=None):
        Digitizer.__init__(self, shn, sock, name='LPS')

        self.IO_mds = IOMdsAUG(shn, sock, diag='LPS')
        self.IO_file = IOFileAUG(shn, diag='LPS')
        self.nodes = ('CUR1', 'VOL1', 'CUR2', 'VOL2', 'VOL3', 'VOL4', 't')

        self.window = slice(2048, None)

        #f = 10./4095
        #offs = (0,0,0,0,0,0,0)
        #self.amp = {node: Amp(fact=f, offs=-f*o) for node, o in zip(self.nodes, offs)}
        self.amp = {node: amp12Bit.copy() for node in self.nodes}
        self.amp['t'] = ampUnity.copy()

        for node in ('CUR1', 'VOL3'):
            self.amp[node] *= ampInv


class ProbeXPR(Probe):
    def __init__(self, shn, sock=None, dig=None):
        self.config = config.campaign.find_shot(shn)
        if dig is None:
            dig = self.config.dig
        
        if dig == 'LPS':
            DigitizerClass = DigitizerLPS
        else:
            DigitizerClass = DigitizerXPR

        digitizer = DigitizerClass(shn, sock)
        Probe.__init__(self, digitizer)

    def mapsig(self):
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

    def calib(self):
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



