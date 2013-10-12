import math
import numpy as np
import numpy.ma as ma
from itertools import izip
from pytokamak.utils.utils import memoized_property, GeneratorDict
from pytokamak.utils.sig import AmpSignal

class DCN:
    def __init__(self, shn, eqi=None, **kw):
        from digitizer_aug import DigitizerAUGDCR
        from diaggeom_aug import DCNGeom

        self.digitizer = DigitizerAUGDCR(shn)
        self.geom = DCNGeom()
        self.eqi = eqi
        self.kw = kw

        self.names = self.geom.names
        self._psi_spl = GeneratorDict(self._get_psi_spl)
        self.sep = GeneratorDict(self.get_separatrix)
  
    def set_kw(self, **kw):
        self.kw = kw

    def __getitem__(self, indx):
        return self.digitizer[indx].t_gt(0.).gt(1e18).compressed()

    def _get_psi_spl(self, chn='H-5'):
        (Rm, zm), (RM, zM) = self.geom.get_rays(chn)
        
        ni = 100
        si = np.linspace(0., 1., ni)
        Ri = np.linspace(Rm, RM, ni)
        zi = np.linspace(zm, zM, ni)
        return self.eqi.get_path_spline(si, Ri, zi, **self.kw)

    def get_separatrix(self, chn='H-5'):
        S = self[chn]
        t = np.asanyarray(S.t, np.float64)
        Spl = (self._psi_spl[chn] - 1.).eval_x(t)

        sep = np.zeros((t.size, 2))
        sep.fill(np.nan)
        for spl, s in izip(Spl, sep):
            spl.roots(s)
        
        kw = dict(units='m', label=chn)
        sep = ma.masked_invalid(sep)
        (Rm, zm), (RM, zM) = self.geom.get_rays(chn)
        dR, dz = RM - Rm, zM - zm
        d = np.sqrt(dR*dR + dz*dz)

        R = AmpSignal(Rm + dR * sep, t, type='Rsep', **kw)
        z = AmpSignal(zm + dz * sep, t, type='zsep', **kw)
        d = AmpSignal(d * (sep[:,1] - sep[:,0]), t, type='dsep', **kw)
        sep = dict(R=R, z=z, d=d)
        return sep

    def get_n(self, chn):
        return self[chn] / self.sep[chn]['d']

    def plot(self, ax=None, chn=None):
        if chn is None:
            chn = self.names
        for c in chn:
            ax = (self[c]*1e-19).plot(ax)
        return ax

    def plot_n(self, ax=None, chn=None):
        if chn is None:
            chn = self.names
        for c in chn:
            ax = (self.get_n(c)*1e-19).plot(ax)
        return ax

    def plot_separatrix(self, ax=None, chn=None, field='d'):
        if chn is None:
            chn = self.names
        for c in chn:
            ax = self.sep[c][field].plot(ax, color='next')
        return ax

    



