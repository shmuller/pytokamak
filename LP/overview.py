import numpy as np
import numpy.ma as ma

from LP import probe_xpr
from LP.sig import memoized_property

from sm_pyplot import tight_figure

figure = tight_figure.pickable_linked_lod_figure


class AUGOverview:
    def __init__(self, shn):
        self.shn = shn

        self.chn = dict(
            DCN = ('H-1', 'H-2', 'H-3', 'H-4', 'H-5'),
            ECS = ('PECRH',),
            NIS = ('PNI',),
            GQI = ('Wmhd',),
            MAC = ('Ipolsola', 'Ipolsoli'))

        self.XPR = probe_xpr.ProbeXPR(shn=shn)

    def load(self, diag):
        IO = probe_xpr.IOMdsAUG(shn=self.shn, diag=diag)
        chn = self.chn[diag] + ('t',)
        try:
            return IO.load(chn)
        except probe_xpr.probe.TdiError:
            return {c: 0. for c in chn}

    @memoized_property
    def S(self):
        return {k: self.load(k) for k in self.chn.keys()}

    def plot(self):
        S, XPR = self.S, self.XPR

        S_DCN  = S['DCN']
        S_ECRH = S['ECS']
        S_NI   = S['NIS']
        S_WMHD = S['GQI']
        S_MAC  = S['MAC']

        fig = figure(figsize=(6,6))

        ax = fig.add_subplot(511)

        ax.set_ylabel('Power (MW)')
        ax.plot(S_ECRH['t'], 1e-6*S_ECRH['PECRH'], label="ECRH")
        ax.plot(S_NI['t'], 1e-6*S_NI['PNI'], label="NBI")
        ax.plot(S_WMHD['t'], 1e-5*S_WMHD['Wmhd'], label="WMHD (x10)")

        ax.legend()

        ax = fig.add_subplot(512)
        ax.set_ylabel('n (10$^{\mathdefault{19}}$ m$^{\mathdefault{-3}}$)')

        t = S_DCN['t']
        for c in ('H-1', 'H-4', 'H-5'):
            n = S_DCN[c]
            n[n < 0] = np.nan
            ax.plot(t, 1e-19*n, label=c)

        ax.legend()

        ax = fig.add_subplot(513)
        ax.set_ylabel('Current (A)')

        t = XPR['tip1'].t
        I1, I2, I3 = XPR['tip1'].x, XPR['tip2'].x, XPR['tip3'].x

        #I1 = ma.masked_array(I1, I1 < 0)
        #I2 = ma.masked_array(I2, I2 < 0)

        ax.plot(t, I1, label='Mach tip 1')
        ax.plot(t, I2, label='Mach tip 2')
        ax.plot(t, I3, label='Single tip')

        #XPR['tip3'].V.plot(ax=ax)

        ax.legend()

        #tips = 'tip1', 'tip2', 'tip3'
        #
        #for tip in tips:
        #    XPR[tip].plot(ax=ax)

        ax = fig.add_subplot(514)
        ax.set_ylabel('Pos (cm)')

        ax.plot(XPR['R'].t, 100*XPR['R'].x)

        ax = fig.add_subplot(515)
        ax.set_ylabel('Current (kA)')

        t = S_MAC['t']
        Ia = ma.masked_array(S_MAC['Ipolsola'], t > 6.)
        Ii = ma.masked_array(S_MAC['Ipolsoli'], t > 6.)

        ax.plot(t, 1e-3*Ia, label='Ipolsola')
        ax.plot(t, 1e-3*Ii, label='Ipolsoli')
        ax.legend()

        for ax in fig.axes:
            ax.set_xlim((1,7))
            ax.grid(True)

        for ax in fig.axes[:-1]:
            ax.set_xticklabels('')

        fig.axes[-1].set_xlabel('t (s)')

        fig.tight_layout(pad=0.2)


