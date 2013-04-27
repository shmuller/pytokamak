import numpy as np

import scipy.interpolate as interp

import config_xpr as config
ShotNotFoundError = config.ShotNotFoundError

from tokamak.digitizer import Digitizer
from tokamak.digitizer_aug import IOMdsAUG, IOFileAUG

from sig import Amp, Signal
from probe import Probe, PhysicalResults

ampUnity = Amp(fact=1., offs=0.)
ampInv   = Amp(fact=-1., offs=0.)
amp12Bit = Amp(fact=10./4095, offs=-5.)
amp14Bit = Amp(fact=20./16383, offs=-10.)


class DigitizerXPR(Digitizer):
    def __init__(self, shn, raw=False):
        Digitizer.__init__(self, shn, name='XPR')

        self.IO_mds = IOMdsAUG(shn, diag='XPR', raw=raw)
        self.IO_file = IOFileAUG(shn, suffix='_XPR')
        self.nodes = ('S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 't')

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


class DigitizerXPRRaw(DigitizerXPR):
    def __init__(self, shn):
        DigitizerXPR.__init__(self, shn, raw=True)

        self.amp = {node: amp14Bit.copy() for node in self.nodes[:-1]}


class DigitizerXPRRawPos(DigitizerXPRRaw):
    def __init__(self, shn):
        DigitizerXPRRaw.__init__(self, shn)
        self.nodes += ('PosL', 'PosH')

    def _load_raw_factory(name):
        def load_raw(self, **kw):
            getattr(DigitizerXPRRaw, name)(self, **kw)

            PosL = self.x['PosL'].astype(np.int16).view(np.uint16)
            PosH = self.x['PosH'].astype(np.int16).view(np.uint16) & 0b0011111111111111
        
            Pos  = np.c_[PosL, np.r_[PosH[0], PosH[:-1]]]
        
            self.x['Pos'] = Pos.copy().view(np.uint32)[:,0].astype(np.float32)
            return self.x
        return load_raw

    load_raw      = _load_raw_factory('load_raw')
    load_raw_mds  = _load_raw_factory('load_raw_mds')
    load_raw_file = _load_raw_factory('load_raw_file')


class DigitizerLPS(Digitizer):
    def __init__(self, shn, raw=False):
        Digitizer.__init__(self, shn, name='LPS')

        self.IO_mds = IOMdsAUG(shn, diag='LPS', raw=raw)
        self.IO_file = IOFileAUG(shn, suffix='_LPS')
        self.nodes = ('CUR1', 'VOL1', 'CUR2', 'VOL2', 'VOL3', 'VOL4', 't')

        self.window = slice(2048, None)


class DigitizerLPSRaw(DigitizerLPS):
    def __init__(self, shn):
        DigitizerLPS.__init__(self, shn, raw=True)

        self.amp = {node: amp12Bit.copy() for node in self.nodes[:-1]}

        for node in ('CUR1', 'VOL3'):
            self.amp[node] *= ampInv


class DigitizerXPOS(Digitizer):
    def __init__(self, shn):
        Digitizer.__init__(self, shn, name='XPOS')

        self.IO_mds = IOMdsAUG(shn, diag='LPS')
        self.IO_file = IOFileAUG(shn, suffix='_LPS_XPOS')
        self.nodes = ('XPOS', 't')


class DigitizerLPSOld(DigitizerLPS):
    def __init__(self, shn):
        DigitizerLPS.__init__(self, shn)

        self.dig_xpos = DigitizerXPOS(shn)

        for node in ('CUR1', 'CUR2', 'VOL3'):
            self.amp[node] = ampInv.copy()

    def _load_raw_factory(name):
        def load_raw(self, **kw):
            x = getattr(self.dig_xpos, name)()
            R = Signal(x['XPOS'].astype('d'), x['t'].astype('d'))
            R.despike()

            getattr(DigitizerLPS, name)(self, **kw)
            self.x['XPOS'] = R(self.x['t']).x.astype(np.float32)
            return self.x
        return load_raw

    load_raw      = _load_raw_factory('load_raw')
    load_raw_mds  = _load_raw_factory('load_raw_mds')
    load_raw_file = _load_raw_factory('load_raw_file')

    def save(self):
        DigitizerLPS.save(self)
        self.dig_xpos.save()


DigitizerClasses = dict(
        LPS = DigitizerLPSRaw,
        XPR = DigitizerXPRRaw,
        XPR_pos = DigitizerXPRRawPos,
        LPS_old = DigitizerLPSOld)


class ProbeXPR(Probe):
    def __init__(self, shn, dig=None):
        
        reload(config)
        
        self.config = config.campaign.find_shot(shn)
        if dig is None:
            dig = self.config.dig

        digitizer = DigitizerClasses[dig](shn)
        Probe.__init__(self, digitizer)

    def mapsig(self):
        w = self.digitizer.window
        x = {k: v[w] for k, v in self.x.iteritems()}

        self.S = self.config.mapsig(x, self.digitizer.name)

    def calib(self):
        self.config.calib(self.digitizer.name)

    def calc_res(self, ID='IV'):
        tips = self.config.head.tips

        A = (0.5*tips[0].area, 0.5*tips[1].area, tips[2].area)

        try:
            PP_j  = self.IV['tip1+tip2'].PP[ID][:3].copy()
            PP_jp = self.IV['tip1'].PP[ID][0].copy()
            PP_jm = self.IV['tip2'].PP[ID][0].copy()
            PP_jt = PP_j[0].copy()
            A_j = A[0] + A[1]
        except KeyError:
            PP_j  = self.IV['tip3'].PP[ID][:3].copy()
            PP_jp = self.S['tip1'].as_PP(PP_j)
            PP_jm = self.S['tip2'].as_PP(PP_j)
            PP_jt = self.S['tip1+tip2'].as_PP(PP_j)
            A_j = A[2]

        PP_j[0].c /= A_j
        PP_jp.c /= A[0]
        PP_jm.c /= A[1]
        PP_jt.c /= A[0] + A[1]

        keys = ['j', 'Vf', 'Te', 'jp', 'jm', 'jt']
        c = np.dstack((PP_j.c, PP_jp.c, PP_jm.c, PP_jt.c))

        try:
            PP_j3 = self.IV['tip3'].PP[ID][0].copy()
            PP_j3.c /= A[2]

            keys.append('j3')
            c = np.dstack((c, PP_j3.c))
        except KeyError:
            pass

        self.PP_meas = PP_j.__class__(c, PP_j.x, **PP_j.kw)
        i, meas = self.PP_meas.eval()

        dtype = zip(keys, [np.double]*len(keys))
        meas = meas.view(dtype).reshape(-1).view(np.recarray)

        shn = self.digitizer.shn
        return PhysicalResults(shn, self['Rs'], i, meas)

    def position_calib(self):
        R = self['R']
        sp = interp.UnivariateSpline(R.t, R.x, s=10)
        Rs = sp(R.t)
        return Rs

    def voltage_calib(self):
        self.load(trim=False, calib=False)
        return self['V'].x.mean()
   

def get_dwell_params():
    shots = config.campaign.shots_with_min_stars('*')
    for shot in shots:
        XPR = ProbeXPR(shn=shot)
        tM, RM = XPR.get_dwell_params()
        print "%d: tM = %s s, RM = %s m" % (shot, tM, RM)


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    shn = 29859

    XPR = ProbeXPR(shn=shn)
    fig = XPR.res.plot_R(plunge=1)
    plt.show()



