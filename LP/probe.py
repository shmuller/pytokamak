import numpy as np
import numpy.ma as ma

from warnings import warn

import logging
reload(logging)
logging.basicConfig(level=logging.WARN)
logger = logging

from collections import Iterable, OrderedDict

from tokamak.digitizer import IOH5

try:
    import h5py as H5
except ImportError:
    pass

import fitter_IV
reload(fitter_IV)

FitterIV = fitter_IV.FitterIV
IVContainer = fitter_IV.IVContainer

from sm_pyplot.tight_figure import get_fig, get_tfig, get_axes
from sm_pyplot.annotations import vlines, vrect

from utils.utils import memoized_property, dict_pop, DictView
from utils.sig import Signal, amp_unity, AmpSignal, PeriodPhaseFinder

class PositionSignal(AmpSignal):
    def __init__(self, x, t, amp=None, dtype=np.float64, **kw):
        kw.setdefault('type', 'Position')
        kw.setdefault('units', 'm')

        AmpSignal.__init__(self, x, t, amp, dtype, **kw)
    
        self.baseline_slice = kw.get('baseline_slice', slice(None, 1000))
        self.lvl_fact = kw.get('lvl_fact', 0.1)
        self.dist_threshold = kw.get('dist_threshold', 1000)

    def to_cm(self):
        x = self*100
        x.type = 'Pos'
        x.units = 'cm'
        return x

    def get_baseline(self):
        return self.x[self.baseline_slice]

    def get_crossings(self):
        x = self.get_baseline()
        x0, xM = x.mean(), max(self.x.max(), 0.)
        lvl = x0 + self.lvl_fact*(xM - x0)
        return self.crossings(lvl, self.dist_threshold)

    @memoized_property
    def t_ind(self):
        ind0, ind1, is_rising = self.get_crossings()
        i0, i1 = ind0[is_rising], ind1[~is_rising]
        iM = self.local_argmax(i0, i1)
        return i0, iM, i1

    @memoized_property
    def Nplunges(self):
        return len(self.t_ind[1])

    def tM(self, plunge=None):
        i0, iM, i1 = self.t_ind
        if plunge is not None:
            if not isinstance(plunge, Iterable):
                plunge = [plunge]
            iM = iM[plunge]
        return self.t[iM]

    def region_boundaries(self, inout=None):
        i0, iM, i1 = self.t_ind
        if inout == 'in':
            return np.array([i0, iM])
        elif inout == 'out':
            return np.array([iM, i1])
        else:
            return np.array([i0, i1])

    def plunges(self, plunge=None, inout=None):
        w = self.region_boundaries(inout)
        if plunge is not None:
            if not isinstance(plunge, Iterable):
                plunge = [plunge]
            w = w[:,plunge]
        return w

    def plunge_mask(self, i, plunge=None, inout=None):
        w = self.plunges(plunge, inout)
        if w.shape[1] == 0:
            return np.arange(0)
        else:
            ind0, ind1 = np.searchsorted(i, w)
            return np.concatenate(map(np.arange, ind0, ind1))

    def plot_plunges(self, ax=None, **kw):
        ax = Signal.plot(self, ax, **kw)
        i0, iM, i1 = self.t_ind
        ax.plot(self.t[i0], self.x[i0], 'r*')
        ax.plot(self.t[i1], self.x[i1], 'g*')
        ax.plot(self.t[iM], self.x[iM], 'm*')
        return ax


class VoltageSignal(AmpSignal):
    def __init__(self, x, t, amp=None, dtype=np.float64, **kw):
        kw.setdefault('type', 'Voltage')
        kw.setdefault('units', 'V')

        self.dt = kw.pop('dt', 0.1)
        self.min_ptp = kw.pop('min_ptp', 50)

        AmpSignal.__init__(self, x, t, amp, dtype, **kw)

    @memoized_property
    def PPF(self):
        return PeriodPhaseFinder(self.x)

    @memoized_property
    def iE(self): 
        return self.PPF.iE

    @memoized_property
    def is_swept(self):
        D = self.min_ptp
        cnd = self.t < self.t[0] + self.dt
        return self.x[cnd].ptp() > D and self.PPF.D > D

    def plot_sweeps(self, ax=None, **kw):
        ax = Signal.plot(self, ax, **kw)
        ax.plot(self.t[self.iE], self.x[self.iE], 'r+')
        return ax


class CurrentSignal(AmpSignal):
    def __init__(self, x, t, amp=None, dtype=np.float64, **kw):
        kw.setdefault('type', 'Current')
        kw.setdefault('units', 'A')

        AmpSignal.__init__(self, x, t, amp, dtype, **kw)

        self.V = kw.get('V', None)
        self.C = kw.get('C', None)

    def __getitem__(self, index):
        s = Signal.__getitem__(self, index)
        if self.V:
            s.V = self.V[index]
        return s

    def __add__(self, other):
        return self.__class__(self.x + other.x, self.t, V=self.V,
                              name=self.name + '+' + other.name)

    def masked_Isat(self, Vmax=-150.):
        return self.masked(self.V > Vmax)

    def capa_pickup(self):
        self.dV_dt = self.V.smooth(10).deriv()

        cnd = self.t - self.t[0] < 0.05
        dV_dtc = self.dV_dt.x[cnd]
        N = (dV_dtc*dV_dtc).sum()

        dI = self.copy()
        dI.norm_to_region(cnd)
        self.C = (dI.x[cnd]*dV_dtc).sum()/N

    @property
    def I_capa(self):
        if self.C is None:
            self.capa_pickup()
        I_capa = self.C*self.dV_dt.x
        return self.__class__(I_capa, self.t, V=self.V, name=self.name + '_capa')

    @property
    def I_corr(self):
        return self - self.I_capa


class PhysicalResults:
    def __init__(self, probe, i, meas):
        self.probe, self.i, self.meas = probe, i, meas
        
        self.shn = probe.digitizer.shn
        self.R = probe['Rs']

        self.keys = ('n_cs', 'Mach', 'nv', 'mnv', 'j', 'jsat', 'Vf', 'Te', 'Vp', 
                'cs', 'n', 'v', 'pe', 'pe_tot', 'R', 't', 'Dt')

        def sup(x): return r'$^{\text{%s}}$' % x
        def sub(x): return r'$_{\text{%s}}$' % x

        self.units = dict(
                Mach   = None,
                n_cs   = r'm%s s%s' % (sup('-2'), sup('-1')),
                nv     = r'm%s s%s' % (sup('-2'), sup('-1')),
                mnv    = r'g m%s s%s' % (sup('-2'), sup('-1')),
                j      = r'kA m%s' % sup('-2'),
                jsat   = r'kA m%s' % sup('-2'),
                Vf     = r'V',
                Te     = r'eV',
                Vp     = r'V',
                cs     = r'km s%s' % sup('-1'),
                n      = r'm%s' % sup('-3'),
                v      = r'km s%s' % sup('-1'),
                pe     = r'Pa',
                pe_tot = r'Pa',
                R      = r'cm',
                t      = r's',
                Dt     = r'ms')

        self.texlabels = dict(
                Mach   = r'Mach',
                n_cs   = r'$n c_s$',
                nv     = r'$nv$',
                mnv    = r'$mnv$',
                j      = r'$j$',
                jsat   = r'$j_{sat}$',
                Vf     = r'$V_f$',
                Te     = r'$T_e$',
                Vp     = r'$V_p$',
                cs     = r'$c_s$',
                n      = r'$n$',
                v      = r'$v$',
                pe     = r'$p_e$',
                pe_tot = r'$p_{e,tot}$',
                R      = r'$R$',
                t      = r'$t$',
                Dt     = r'$\Delta t$')

        self.fact = dict.fromkeys(self.keys, 1)
        self.fact['j'] = self.fact['jsat'] = self.fact['cs'] = self.fact['v'] = 1e-3
        self.fact['R'] = 100
        self.fact['mnv'] = self.fact['Dt'] = 1e3

        self.lim = dict.fromkeys(self.keys, (None, None))
        self.lim['n'] = (0, 2e20)
        self.lim['Mach'] = (-2, 2)
        self.lim['j'] = self.lim['jsat'] = (0, None)
        self.lim['Te'] = (0, 100)
        self.lim['R'] = (0, None)

    @memoized_property
    def res(self):
        qe = 1.6022e-19
        mi = 2*1.67e-27

        Gp = self.meas['jp']/qe
        Gm = self.meas['jm']/qe

        Gp[Gp <= 0] = np.nan
        Gm[Gm <= 0] = np.nan

        res = dict()

        R = self.probe.pos[self.i, 0]
        res['t'] = R.t
        res['R'] = R.x

        res['Mach'] = Mach = 0.5*np.log(Gm/Gp)
        res['n_cs'] = n_cs = np.e*np.sqrt(Gp*Gm)
        res['nv']   = nv = n_cs*Mach
        res['mnv']  = mi*nv
        res['j']    = qe*nv
        res['jsat'] = self.meas['jt']
        res['Vf']   = Vf = self.meas['Vf']
        res['Te']   = Te = self.meas['Te']
        res['Vp']   = Vf + 2.8*Te

        Ti = Te
        res['cs'] = cs = np.sqrt(qe/mi*(Te+Ti))
        res['n']  = n  = n_cs/cs
        res['v']  = v  = Mach*cs
        res['pe'] = pe = n*qe*Te
        res['pe_tot']  = pe*(1. + Mach*Mach)
        return res

    def eval(self, plunge=None, inout=None):
        mask = self.R.plunge_mask(self.i, plunge, inout)
        
        y = {k: v[mask] for k, v in self.res.iteritems()}

        y['Dt'] = y['t'].copy()
        tM = self.R.tM(plunge)
        if len(tM) == 1:
            y['Dt'] -= tM[0]
        return y

    def make_name(self, plunge=None, inout=None):        
        if plunge is None:
            name = "XPRres_%d" % self.shn
        else:
            name = "XPRres_%d_%d" % (self.shn, plunge)
        if inout is not None:
            name += "_" + inout
        return name

    def save(self, plunge=None, inout=None):
        name = self.make_name(plunge, inout)
        y = self.eval(plunge, inout)        
        IOH5(name + '.h5').save(y)

    def load(self, plunge=None, inout=None):
        name = self.make_name(plunge, inout)
        return IOH5(name + '.h5').load()

    def make_label(self, key, usetex=False):
        if usetex:
            lab = self.texlabels[key]
        else:
            lab = key
        if self.units[key] is not None:
            lab += r' (' + self.units[key] + r')'
        return lab

    def clip(self, y, lim):
        if lim == (None, None):
            return y
        else:
            mask = np.zeros_like(y, dtype=bool)
            if lim[0] is not None:
                mask[y < lim[0]] = True
            if lim[1] is not None:
                mask[y > lim[1]] = True
            return ma.masked_array(y, mask)

    def plot_key(self, key, x, y, ax=None, label=None):
        ax = get_axes(ax)

        ylab = self.make_label(key, ax.figure.usetex)
        ax.set_ylabel(ylab)

        yc = self.clip(y[key], self.lim[key])

        ax.plot(x, self.fact[key]*yc, label=label)

    def plot(self, fig=None, keys=None, xkey='t', plunge=None, inout=None, 
            mirror=False, figsize=(10, 10), usetex=False, **kw):
        y = self.eval(plunge, inout)
                
        x = self.fact[xkey]*self.clip(y[xkey], self.lim[xkey])
        if mirror:
            x = -x
        xlab = self.make_label(xkey, usetex)

        label = "%d" % self.shn
        tM = self.R.tM(plunge)
        if len(tM) == 1:
            label += ".%d" % (1e3*tM[0])
        if inout == 'out':
            label += " (out)"

        if keys is None:
            keys = (('Dt', 'R'     ), 
                    ('n' , 'Mach'  ), 
                    ('Vf', 'mnv'   ), 
                    ('Te', 'pe'    ), 
                    ('Vp', 'pe_tot'))
        keys = np.array(keys, ndmin=2)

        if xkey == 't':
            viewers = self.probe.viewers
        else:
            viewers = ()

        fig = get_tfig(fig, keys.shape, xlab=xlab, 
                       figsize=figsize, usetex=usetex, viewers=viewers)

        for key, ax in zip(keys.ravel(), fig.axes):
            self.plot_key(key, x, y, ax=ax, label=label)

        self.probe.plot_separatrix_crossings(fig.axes, xkey=xkey, 
                fact=self.fact[xkey], plunge=plunge, inout=inout, **kw)

        fig.axes[0].legend(loc='upper left')
        fig.canvas.draw()
        return fig

    def plot_R(self, **kw):
        kw['fig'] = self.plot_R_in(**kw)
        return self.plot_R_out(**kw)

    def plot_R_in(self, **kw):
        return self.plot(xkey='R', inout='in', **kw)

    def plot_R_out(self, **kw):
        return self.plot(xkey='R', inout='out', **kw)

    def _plot_all_factory(name):
        def plot_all(self, **kw):
            plotfun = getattr(self, name)
            kw['fig'] = plotfun(plunge=0, **kw)
            for plunge in xrange(1, self.R.Nplunges):
                plotfun(plunge=plunge, **kw)
            return kw['fig']
        return plot_all

    plot_R_all     = _plot_all_factory('plot_R')
    plot_R_all_in  = _plot_all_factory('plot_R_in')
    plot_R_all_out = _plot_all_factory('plot_R_out')


class Probe:
    def __init__(self, head, digitizer, R0, z0, eqi=None, viewers=()):
        self.head, self.digitizer, self.R0, self.z0 = head, digitizer, R0, z0
        self.eqi, self.viewers = eqi, viewers

    def __getitem__(self, index):
        return self.S[index]
   
    @memoized_property
    def S(self):
        self.S = self.mapsig()
        self.calib()
        return self.S

    @property
    def pos(self):
        raise NotImplementedError

    @memoized_property
    def psi(self):
        x, t = self.pos.x, self.pos.t
        return self.eqi(t, *x.T)

    def get_keys(self, name):
        raise NotImplementedError

    def get_mapping(self, key):
        raise NotImplementedError

    def get_amp(self, key):
        raise NotImplementedError

    def get_sig(self, key):
        return self.digitizer[self.get_mapping(key)]

    def mapsig(self):
        x = self.get_sig('ampR')
        nans = np.zeros(x.shape, np.float32)
        nans.fill(np.nan)
        Nans = AmpSignal(nans, x._t)

        def get_sig_amp(key):
            try:
                return self.get_sig(key), self.get_amp(key)
            except KeyError:
                return Nans, amp_unity
        
        x, amp = get_sig_amp('ampR')
        R = PositionSignal(x._x, x._t, amp*x.amp, name='R')

        S = OrderedDict(R=R)

        for tip in self.head.tips:
            i = tip.number
            keys = self.get_keys(tip.name)

            x, amp = get_sig_amp(keys['V'])
            V = VoltageSignal(x._x, x._t, amp*x.amp, 
                              number=i, name='V%d' % i, label=tip.label)
            
            x, amp = get_sig_amp(keys['I'])
            I = CurrentSignal(x._x, x._t, amp*x.amp, V=V, 
                              number=i, name='I%d' % i, label=tip.label)

            S[tip.name] = I

        return S

    def calib(self):
        pass

    def norm_to_region(self, s):
        for tip in self.head.tips:
            S = self.S[tip.name]
            S.norm_to_region(s)

            # I_keys is None: floating potential
            if self.get_keys(tip.name)['I'] is None:
                S.V.norm_to_region(s)
    
    def get_type(self, type):
        keys = [k for k, v in self.S.iteritems() if v.type == type]
        return DictView(self.S, keys)

    @memoized_property
    def I(self):
        return self.get_type('Current')

    @memoized_property
    def V(self):
        return OrderedDict(zip(self.I.keys(), [I.V for I in self.I.itervalues()]))
        
    @memoized_property
    def I_swept(self):
        keys = [k for k, I in self.I.iteritems() if I.V.is_swept]
        return DictView(self.I, keys)

    def _plot(self, S_list, ax=None, unmasked=True, sepmode='lines', 
              legend_loc='upper right', **kw):
        S = S_list[0]
        if ax is None:
            fig = get_tfig(xlab=S.xlab, viewers=self.viewers, **kw)
            ax = fig.axes[0]
        ax.set_ylabel(S.ylab)

        if unmasked:
            for S in S_list:
                S.unmasked().plot(ax)
        else:
            for S in S_list:
                S.plot(ax)

        if sepmode is not None:
            self.plot_separatrix_crossings([ax], sepmode=sepmode, color='k')
        if legend_loc is not None:
            ax.legend(loc=legend_loc)
        return ax

    def plot_R(self, *args, **kw):
        return self._plot([self.S['Rs'].to_cm()], *args, **kw)

    def plot_I(self, *args, **kw):
        return self._plot(self.I.values(), *args, **kw)

    def plot_V(self, *args, **kw):
        return self._plot(self.V.values(), *args, **kw)

    def plot(self, fig=None, **kw):
        kw2 = dict_pop(kw, unmasked=True, sepmode='lines', )
        R = self.S['Rs']
        fig = get_fig(fig, shape=(3, 1), xlab=R.xlab, viewers=self.viewers, **kw)
        
        self.plot_R(ax=fig.axes[0], legend_loc=None, **kw2)
        self.plot_I(ax=fig.axes[1], legend_loc='upper left', **kw2)
        self.plot_V(ax=fig.axes[2], legend_loc='upper left', **kw2)
        return fig

    def get_dwell_params(self):
        R = self.S['Rs']
        iM = R.t_ind[1]
        tM = R.t[iM]
        RM = R.x[iM]
        return tM, RM

    def trim(self, plunge='all'):
        R = self.S['Rs']
        i0, iM, i1 = R.t_ind

        if plunge == 'all':
            s = np.concatenate(map(np.arange, i0, i1))
        else:
            s = slice(i0[plunge], i1[plunge])

        for S in self.S.itervalues():
            S.trim(s)

    @memoized_property
    def t0(self):
        return self.digitizer.x['t'][0]

    def shift_t(self, plunge=None, phase=1):
        t0 = self.t0
        R = self.S['Rs']
        if plunge is None:
            R.t -= R.t[0] - t0
        else:
            R.t -= R.t[R.t_ind[phase][plunge]]
        return self

    def corr_capa(self):
        for I in self.I.itervalues():
            I.x[:] -= I.I_capa

    def calc_IV(self):
        return IVContainer(self.I_swept, self.S['Rs'])

    @memoized_property
    def IV(self):
        return self.calc_IV()
    
    def calc_res(self, ID='IV'):
        pass

    @memoized_property
    def res(self):
        return self.calc_res()

    def plot_head(self, ti, ax=None):
        ax = get_axes(ax)
        R, z = self.pos(ti, masked=True).x[0]
        self.head.plot(ax, R, z)
        return ax

    def _get_xsep(self, xkey='t', fact=1., offs=0., plunge=None, inout=None, **kw):
        isep = self.psi.separatrix_crossings
        mask = self.S['Rs'].plunge_mask(isep, plunge, inout)
        psep = self.pos[isep[mask]]
        if xkey == 't':
            xsep = fact*psep.t
        else:
            xsep = fact*psep.x[:,0]
        return xsep + offs

    def plot_separatrix_crossings(self, axes, color='last', linewidth=1, **kw):
        try:
            xsep = self._get_xsep(**kw)
            for ax in axes:
                vlines(ax, xsep, color=color, linewidth=linewidth)
        except:
            warn("Could not generate lines")


