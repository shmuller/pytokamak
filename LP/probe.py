import numpy as np
import numpy.ma as ma

import logging
reload(logging)
logging.basicConfig(level=logging.WARN)
logger = logging

from collections import OrderedDict

import fitter_IV
reload(fitter_IV)

FitterIV = fitter_IV.FitterIV
IVContainer = fitter_IV.IVContainer

from sm_pyplot.tight_figure import get_fig, get_tfig, get_axes

from sig import memoized_property, DictView, math_sel, usetex, \
        PositionSignal, VoltageSignal, CurrentSignal


class PhysicalResults:
    def __init__(self, shn, R, i, meas, usetex=usetex):
        self.shn, self.R, self.i, self.meas, self.usetex = shn, R, i, meas, usetex

        self.keys = ('n_cs', 'Mach', 'nv', 'mnv', 'j', 'Vf', 'Te', 'Vp', 
                'cs', 'n', 'v', 'pe', 'pe_tot', 'R', 't', 'Dt')

        def sup(x): return r'$^{\mathdefault{%s}}$' % x
        def sub(x): return r'$_{\mathdefault{%s}}$' % x

        self.units = dict(
                n_cs   = r'm%s s%s' % (sup('-2'), sup('-1')),
                Mach   = None,
                nv     = r'm%s s%s' % (sup('-2'), sup('-1')),
                mnv    = r'g m%s s%s' % (sup('-2'), sup('-1')),
                j      = r'kA m%s' % sup('-2'),
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
                n_cs   = math_sel.wrap(r'n c_s'),
                Mach   = r'Mach',
                nv     = math_sel.wrap(r'nv'),
                mnv    = math_sel.wrap(r'mnv'),
                j      = math_sel.wrap(r'j'),
                Vf     = math_sel.wrap(r'V_f'),
                Te     = math_sel.wrap(r'T_e'),
                Vp     = math_sel.wrap(r'V_p'),
                cs     = math_sel.wrap(r'c_s'),
                n      = math_sel.wrap(r'n'),
                v      = math_sel.wrap(r'v'),
                pe     = math_sel.wrap(r'p_e'),
                pe_tot = math_sel.wrap(r'p_{e,tot}'),
                R      = math_sel.wrap(r'R'),
                t      = math_sel.wrap(r't'),
                Dt     = math_sel.wrap(r'\Delta t'))

        self.fact = dict.fromkeys(self.keys, 1)
        self.fact['j'] = self.fact['cs'] = self.fact['v'] = 1e-3
        self.fact['R'] = 100
        self.fact['mnv'] = self.fact['Dt'] = 1e3

        self.lim = dict.fromkeys(self.keys, (None, None))
        self.lim['n'] = (0, 2e20)
        self.lim['Mach'] = (-2, 2)
        self.lim['Te'] = (0, 100)
        self.lim['R'] = (0, None)

    @memoized_property
    def res(self):
        qe = 1.6022e-19
        mi = 2*1.67e-27
        R0 = 1.645

        Gp = self.meas.jp/qe
        Gm = self.meas.jm/qe

        Gp[Gp <= 0] = np.nan
        Gm[Gm <= 0] = np.nan

        dtype = zip(self.keys, [np.double]*len(self.keys))
        res = np.empty(Gp.size, dtype).view(np.recarray)

        res.t = self.R.t[self.i]
        res.R = R0 - self.R.x[self.i]

        res.Mach = Mach = 0.5*np.log(Gm/Gp)
        res.n_cs = n_cs = np.e*np.sqrt(Gp*Gm)
        res.nv  = nv = n_cs*Mach
        res.mnv = mi*nv
        res.j   = qe*nv
        res.Vf  = Vf = self.meas.Vf
        res.Te  = Te = self.meas.Te
        res.Vp  = Vf + 2.8*Te

        Ti = Te
        res.cs = cs = np.sqrt(qe/mi*(Te+Ti))
        res.n  = n  = n_cs/cs
        res.v  = v  = Mach*cs
        res.pe = pe = n*qe*Te
        res.pe_tot  = pe*(1. + Mach*Mach)
        return res 

    def _mask(self, w):
        ind0, ind1 = np.searchsorted(self.i, w)
        return np.concatenate(map(np.arange, ind0, ind1))

    def eval(self, plunge=None, inout=None):
        w = self.R.plunges(plunge, inout)
        
        y = self.res[self._mask(w)]

        y.Dt = y.t
        tM = self.R.tM(plunge)
        if len(tM) == 1:
            y.Dt -= tM[0]
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

        tM = self.R.tM()
        if plunge is not None:
            tM = tM[plunge]
        
        f = h5py.File(name + '.h5', "w")
        for key in y.dtype.names:
            f.create_dataset(key, data=y[key], compression="gzip")
        f.close()

    def load(self, plunge=None, inout=None):
        name = self.make_name(plunge, inout)
        
        f = h5py.File(name + '.h5', "r")
        keys = [key.encode('ascii') for key in f.keys()]
        type = [f[key].dtype for key in keys]
        size = [f[key].len() for key in keys]
        
        y = np.empty(size[0], zip(keys, type))
        for key in keys:
            y[key] = f[key][:]
        f.close()
        return y.view(np.recarray)

    def make_label(self, key):
        if self.usetex:
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

        ylab = self.make_label(key)
        ax.set_ylabel(ylab)

        yc = self.clip(y[key], self.lim[key])

        ax.plot(x, self.fact[key]*yc, label=label)

    def plot(self, fig=None, keys=None, xkey='t', plunge=None, inout=None, 
            mirror=False, figsize=(10, 10)):
        y = self.eval(plunge, inout)
                
        x = self.fact[xkey]*self.clip(y[xkey], self.lim[xkey])
        if mirror:
            x = -x
        xlab = self.make_label(xkey)

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

        fig = get_tfig(fig, keys.shape, xlab=xlab, figsize=figsize)

        ax = fig.axes
        for i in xrange(keys.size):
            self.plot_key(keys.flat[i], x, y, ax=ax[i], label=label)

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
    def __init__(self, head, digitizer):
        self.head, self.digitizer = head, digitizer

        self.xlab = "t (s)"
        self.ylab = ("Isat (A)", "Vf (V)", "Te (eV)")

    def __getitem__(self, index):
        return self.S[index]
   
    @memoized_property
    def x(self):
        return self.digitizer.x

    @memoized_property
    def S(self):
        self.S = self.mapsig()
        self.calib()
        return self.S

    def get_mapping(self, key):
        raise NotImplementedError

    def get_amp(self, key):
        raise NotImplementedError

    @memoized_property
    def unique_sigs(self):
        return {k: self.x[self.get_mapping(k)].astype('d') 
                for k in self.head.unique_keys}

    def mapsig(self):
        t = self.x['t'].astype('d')
        R = self.unique_sigs[self.head.R_keys]
        S = OrderedDict(R=PositionSignal(R, t, name='R'))

        nans = np.zeros_like(t)
        nans.fill(np.nan)

        def get_sig(key):
            try:
                return self.unique_sigs[key]
            except KeyError:
                return nans

        for tip in self.head.tips:
            i = tip.number
            x_V = get_sig(tip.V_keys)
            x_I = get_sig(tip.I_keys)

            V = VoltageSignal(x_V, t, number=i, name='V%d' % i)
            I = CurrentSignal(x_I, t, V=V, number=i, name='I%d' % i)
            S[tip.name] = I

        return S

    def calib(self):
        for k, v in self.unique_sigs.iteritems():
            self.get_amp(k).apply(v)

    def norm_to_region(self, s):
        for tip in self.head.tips:
            S = self.S[tip.name]
            S.norm_to_region(s)

            # I_keys is None: floating potential
            if tip.I_keys is None:
                S.V.norm_to_region(s)

    def load_raw(self, loadfun='load', plunge=None, calib=True, corr_capa=False, **kw):
        self.x = getattr(self.digitizer, loadfun)(**kw)

        self.mapsig()
        if calib:
            self.calib()
        if corr_capa:
            self.corr_capa()
        if plunge is not None:
            self.trim(plunge)
    
    def get_type(self, type):
        keys = [k for k, v in self.S.iteritems() if v.type == type]
        return DictView(self.S, keys)

    @memoized_property
    def I(self):
        return self.get_type('Current')

    @memoized_property
    def V(self):
        return self.get_type('Voltage')

    @memoized_property
    def R(self):
        return self.get_type('Position')

    @memoized_property
    def I_swept(self):
        keys = [k for k, I in self.I.iteritems() if I.V.is_swept]
        return DictView(self.I, keys)

    def plot_raw(self, fig=None,
            keys = (('Position',), ('Current',), ('Voltage',))):
        keys = np.array(keys, ndmin=2)

        fig = get_fig(fig, keys.shape, xlab=self.xlab, ylab=keys)
        axes = np.array(fig.axes).reshape(keys.shape)

        for ax in axes[keys == 'Voltage']:
            for I in self.I.itervalues():
                I.V.plot(ax)

        for key, ax in zip(keys, fig.axes):
            for S in self.get_type(key[0]).itervalues():
                S.plot(ax)

        fig.canvas.draw()
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
        return self.x['t'][0]

    def shift_t(self, plunge=None, phase=1):
        t0 = self.t0
        R = self.S['Rs']
        if plunge is None:
            R.t -= R.t[0] - t0
        else:
            R.t -= R.t[R.t_ind[phase][plunge]]
        return self

    def smooth_I(self, w=10):
        for I in self.I.itervalues():
            I.smooth(w)

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

