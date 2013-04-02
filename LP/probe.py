import numpy as np

import scipy.interpolate as interp

import warnings

import logging
reload(logging)
logging.basicConfig(level=logging.WARN)
logger = logging

from pdb import set_trace

from sm_pyplot.contextmenupicker import ContextMenuPicker
from sm_pyplot.observer_viewer import ToggleViewer, ToggleViewerIntegrated

from sig import *

import fitter_IV
reload(fitter_IV)

FitterIV = fitter_IV.FitterIV
FitterError = fitter_IV.FitterError
IVContainer = fitter_IV.IVContainer


class PhysicalResults:
    def __init__(self, shn, R, i, meas, usetex=usetex):
        self.shn, self.R, self.i, self.meas, self.usetex = shn, R, i, meas, usetex

        self.keys = ('n_cs', 'Mach', 'nv', 'mnv', 'j', 'Vf', 'Te', 'Vp', 
                'cs', 'n', 'v', 'pe', 'R', 't', 'Dt')

        def sup(x): return r'$^{\mathdefault{%s}}$' % x
        def sub(x): return r'$_{\mathdefault{%s}}$' % x

        self.units = dict(
                n_cs = r'm%s s%s' % (sup('-2'), sup('-1')),
                Mach = None,
                nv   = r'm%s s%s' % (sup('-2'), sup('-1')),
                mnv  = r'g m%s s%s' % (sup('-2'), sup('-1')),
                j    = r'kA m%s' % sup('-2'),
                Vf   = r'V',
                Te   = r'eV',
                Vp   = r'V',
                cs   = r'km s%s' % sup('-1'),
                n    = r'm%s' % sup('-3'),
                v    = r'km s%s' % sup('-1'),
                pe   = r'Pa',
                R    = r'cm',
                t    = r's',
                Dt   = r'ms')

        self.texlabels = dict(
                n_cs = math_sel.wrap(r'n c_s'),
                Mach = r'Mach',
                nv   = math_sel.wrap(r'nv'),
                mnv  = math_sel.wrap(r'mnv'),
                j    = math_sel.wrap(r'j'),
                Vf   = math_sel.wrap(r'V_f'),
                Te   = math_sel.wrap(r'T_e'),
                Vp   = math_sel.wrap(r'V_p'),
                cs   = math_sel.wrap(r'c_s'),
                n    = math_sel.wrap(r'n'),
                v    = math_sel.wrap(r'v'),
                pe   = math_sel.wrap(r'p_e'),
                R    = math_sel.wrap(r'R'),
                t    = math_sel.wrap(r't'),
                Dt   = math_sel.wrap(r'\Delta t'))

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
        res.n  = n = n_cs/cs
        res.v  = v = Mach*cs
        res.pe = n*qe*Te
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
        ax = get_axes(ax, figure=tfigure)

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
            keys = ('Dt', 'R'), ('n', 'Mach'), ('Vf', 'v'), ('Te', 'mnv'), ('Vp', 'pe')
        keys = np.array(keys, ndmin=2)

        fig = get_tfig(fig, keys.shape, xlab=xlab, figsize=figsize)

        ax = fig.axes
        for i in xrange(keys.size):
            self.plot_key(keys.flat[i], x, y, ax=ax[i], label=label)

        fig.axes[0].legend(loc='upper left')
        fig.canvas.draw()
        return fig

    def plot_R_in(self, **kw):
        return self.plot(xkey='R', inout='in', **kw)

    def plot_R_out(self, **kw):
        return self.plot(xkey='R', inout='out', **kw)

    def plot_R(self, **kw):
        kw['fig'] = self.plot_R_in(**kw)
        return self.plot_R_out(**kw)


class Probe:
    def __init__(self, digitizer=None):
        self.digitizer = digitizer
        self.PP = None

        self.xlab = "t (s)"
        self.ylab = ("Isat (A)", "Vf (V)", "Te (eV)")

    def __getitem__(self, index):
        return self.S[index]
        
    def mapsig(self):
        pass

    def calib(self):
        pass
    
    @memoized_property
    def x(self):
        return self.digitizer.x

    @memoized_property
    def S(self):
        self.mapsig()
        self.calib()
        return self.S

    def load_raw(self, loadfun='load', plunge=None, calib=True, corr_capa=False):
        self.x = getattr(self.digitizer, loadfun)()

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

