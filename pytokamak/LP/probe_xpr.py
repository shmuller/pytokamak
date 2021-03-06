import numpy as np
from warnings import warn

from sm_pyplot.annotations import vrect

import config_xpr as config
ShotNotFoundError = config.ShotNotFoundError

from pytokamak.tokamak.digitizer_aug import DigitizerAUG
from pytokamak.tokamak.equilibrium import NormalizedFluxSignal

from pytokamak.utils.utils import memoized_property
from pytokamak.utils.sig import Amp, amp_inv, get_fig, get_tfig
from probe import PositionSignal, Probe, PhysicalResults

amp_12bit = Amp(fact=10./4095, offs=-5.)
amp_14bit = Amp(fact=20./16383, offs=-10.)


class DigitizerXPR(DigitizerAUG):
    def __init__(self, shn, raw=False, **kw):
        DigitizerAUG.__init__(self, shn, diag='XPR', suffix='_XPR', group='', raw=raw,
            nodes=('S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8'), **kw)

    @memoized_property
    def window(self):
        '''Detect failing data readout and cut window
        '''
        t = self.x[self.tnode]
        dt = np.abs(np.diff(t))
        cnd = dt > 2*dt[0]
        iM = np.argmax(cnd)

        stop = None
        if cnd[iM] == True:
            stop = iM+1
        return slice(None, stop)


class DigitizerXPRRaw(DigitizerXPR):
    def __init__(self, shn):
        DigitizerXPR.__init__(self, shn, raw=True)

        self.amp = {node: amp_14bit.copy() for node in self.nodes}


class DigitizerXPRRawPos(DigitizerXPRRaw):
    def __init__(self, shn):
        DigitizerXPRRaw.__init__(self, shn)
        self.nodes += ('PosL', 'PosH')

    def __getitem__(self, node):
        try:
            return DigitizerXPRRaw.__getitem__(self, node)
        except KeyError:
            if node == 'Pos':
                L = self.x['PosL'].astype(np.int16).view(np.uint16)
                H = self.x['PosH'].astype(np.int16).view(np.uint16) & 0b0011111111111111
                x = np.c_[L, np.r_[H[0], H[:-1]]].copy().view(np.uint32)[:,0]
                return self.make_sig(node, x, self.x[self.tnode])
            else:
                raise
            

class DigitizerLPS(DigitizerAUG):
    def __init__(self, shn, raw=False):
        DigitizerAUG.__init__(self, shn, diag='LPS', suffix='_LPS', group='', raw=raw,
            nodes=('CUR1', 'VOL1', 'CUR2', 'VOL2', 'VOL3', 'VOL4'))

        self.window = slice(2048, None)


class DigitizerLPSRaw(DigitizerLPS):
    def __init__(self, shn):
        DigitizerLPS.__init__(self, shn, raw=True)

        self.amp = {node: amp_12bit.copy() for node in self.nodes}

        for node in ('CUR1', 'VOL3'):
            self.amp[node] *= amp_inv


class DigitizerLPSOld(DigitizerLPS):
    def __init__(self, shn):
        DigitizerLPS.__init__(self, shn)
        self.tnode = 'TIME2'

        for node in ('CUR1', 'CUR2', 'VOL3'):
            self.amp[node] = amp_inv.copy()

        self.dig_xpos = DigitizerAUG(shn, diag='LPS', suffix='_LPS', group='XPOS',
                nodes=('XPOS',))

    def __getitem__(self, node):
        try:
            return DigitizerLPS.__getitem__(self, node)
        except KeyError:
            if node == 'XPOS':
                R = self.dig_xpos['XPOS'].despike()
                x = R(self.x[self.tnode]).x.astype(np.float32)
                return self.make_sig(node, x, self.x[self.tnode])
            else:
                raise


DigitizerClasses = dict(
        LPS = DigitizerLPSRaw,
        XPR = DigitizerXPRRaw,
        XPR_pos = DigitizerXPRRawPos,
        LPS_old = DigitizerLPSOld)


class ProbeXPR(Probe):
    def __init__(self, shn, shot=None, head=None, dig=None, eqi=None):
        if shot is None:
            shot = self.find_shot(shn)
        if head is None:
            head = shot.head
        if dig is None:
            dig = shot.dig

        self.shot = shot
        digitizer = DigitizerClasses[dig](shn)
        try:
            viewers = eqi.get_viewers(self)
        except AttributeError:
            viewers = ()
        
        Probe.__init__(self, head, digitizer, R0=1.645, z0=-0.966, 
                eqi=eqi, viewers=viewers)

    @classmethod
    def find_shot(cls, shn):
        reload(config)
        return config.campaign.find_shot(shn)

    @property
    def pos(self):
        t = self.S['R'].t
        Rz = np.empty((t.size, 2))
        Rz[:,0] = self.R0 - self.S['R'].x
        Rz[:,1] = self.z0
        return PositionSignal(Rz, t, name='Rz')

    """
    @memoized_property
    def _psi_spl(self):
        R, z = self.pos.x.T
        Rm, RM = R.min(), R.max()
        if Rm == RM: Rm -= 0.01
        Ri = np.linspace(Rm, RM, 20)
        zi = np.array([self.z0]).repeat(Ri.size)
        return self.eqi.get_path_spline(Ri, Ri, zi)

    @memoized_property
    def psi(self):
        t = self.pos.t
        R, z = self.pos.x.T
        psi = self._psi_spl.ev(t, R)
        return NormalizedFluxSignal(psi, t)
    """
    def get_keys(self, name):
        return self.shot.tipmap[name]

    def get_mapping(self, key):
        return self.shot.get(self.digitizer.name, 'mapping', key)

    def get_amp(self, key):
        return self.shot.get(self.digitizer.name, 'amp', key)

    def calib(self):
        Probe.calib(self)

        s = slice(5000)
        self.norm_to_region(s)
        
        S = self.S
        S['Rs'] = S['R'].smooth(100, mode='gaussian')
        S['tip1+tip2'] = S['tip1'] + S['tip2']
        S['tip1+tip2'].update(label='Mach tips sum')
        
        if not S['tip1+tip2'].V.is_swept:
            mask = (S['tip1+tip2'].V > -150.) | (S['tip1+tip2'] > 1.9)
            S['tip1+tip2'] = S['tip1+tip2'].masked(mask)
            S['tip1'] = S['tip1'].masked(mask)
            S['tip2'] = S['tip2'].masked(mask)

        tips = self.head.tips
        A = [0.5*tips[0].area, 0.5*tips[1].area, tips[2].area]
        A.append(A[0] + A[1])
        for tip, area in zip(('tip1', 'tip2', 'tip3', 'tip1+tip2'), A):
            S[tip + 'j'] = S[tip] / area
            S[tip + 'j'].update(type='Current density', units='A / m**2')

    def calc_res(self, ID='IV'):
        tips = self.head.tips

        A = (0.5*tips[0].area, 0.5*tips[1].area, tips[2].area)

        self.PP = PP = dict()
        if self.IV.x.has_key('tip1+tip2'):
            PP['j'], PP['Vf'], PP['Te'] = self.IV['tip1+tip2'].PP[ID][:3].copy()
            PP['jp'] = self.IV['tip1'].PP[ID][0].copy()
            PP['jm'] = self.IV['tip2'].PP[ID][0].copy()
            PP['jt'] = PP['j'].copy()
            A_j = A[0] + A[1]
        else:
            PP['j'], PP['Vf'], PP['Te'] = self.IV['tip3'].PP[ID][:3].copy()
            PP['jp'] = self.S['tip1'].as_PP(PP['j'])
            PP['jm'] = self.S['tip2'].as_PP(PP['j'])
            PP['jt'] = self.S['tip1+tip2'].as_PP(PP['j'])
            A_j = A[2]

        PP['j'].c /= A_j
        PP['jp'].c /= A[0]
        PP['jm'].c /= A[1]
        PP['jt'].c /= A[0] + A[1]

        i = PP['j'].eval()[0]
        meas = {k: v.eval()[1] for k, v in PP.iteritems()}
        
        if self.IV.x.has_key('tip3'):
            PP['j3'] = self.IV['tip3'].PP[ID][0].copy()
            PP['j3'].c /= A[2]
            i3, meas_j3 = PP['j3'].eval()
            if np.all(i == i3):
                meas['j3'] = meas_j3
            else:
                meas['j3'] = PP['j3'](PP['j'].x[i])

        return PhysicalResults(self, i, meas)    

    def _plot(self, *args, **kw):
        kw.setdefault('sepmode', 'patch')
        return Probe._plot(self, *args, **kw)

    def plot(self, *args, **kw):
        kw.setdefault('sepmode', 'patch')
        return Probe.plot(self, *args, **kw)

    def plot_I_Mach(self, *args, **kw):
        return self._plot([self.S['tip1'], self.S['tip2']], *args, **kw)

    def plot_I_single(self, *args, **kw):
        return self._plot([self.S['tip3']], *args, **kw)

    def plot_V_Mach(self, *args, **kw):
        return self._plot([self.S['tip1'].V, self.S['tip2'].V], *args, **kw)

    def plot_V_single(self, *args, **kw):
        return self._plot([self.S['tip3'].V], *args, **kw)

    def plot_separatrix_crossings(self, axes, sepmode='patch', **kw):
        if sepmode == 'patch':
            try:
                xsep = self._get_xsep(**kw).reshape((-1, 2))
                for ax in axes:
                    for xlim in xsep:
                        vrect(ax, xlim, **kw)
            except:
                warn("Could not generate patches, falling back to lines")
                Probe.plot_separatrix_crossings(self, axes, **kw)
        elif sepmode == 'lines':
            Probe.plot_separatrix_crossings(self, axes, **kw)
        else:
            raise ValueError("Unknown sepmode: '%s'" % sepmode)

    def mk_R_axis(self, axes, Rj=None, tlab=None, xlab='Probe R (cm)', **kw):
        if tlab is None:
            tlab = np.r_[np.arange(131, 135), np.arange(135, 160, 5)]
        if Rj is None:
            Rj = 0.01*tlab
        return Probe.mk_R_axis(self, axes, Rj, tlab=tlab, xlab=xlab, **kw)


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



