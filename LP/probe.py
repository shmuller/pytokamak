import numpy as np

import scipy.interpolate as interp

import warnings

import logging
reload(logging)
logging.basicConfig(level=logging.WARN)
logger = logging

from pdb import set_trace

from sig import *

from fitter_IV import FitterIV, FitterError, IVSeriesSimpleGroup

from sm_pyplot.contextmenupicker import ContextMenuPicker
from sm_pyplot.observer_viewer import ToggleViewer, ToggleViewerIntegrated


class ArrayView(np.ndarray):
    def __new__(subtype, x, fields):
        dtype = {f: x.dtype.fields[f] for f in fields}
        return np.ndarray.__new__(subtype, x.shape, dtype,
                                  buffer=x, strides=x.strides)


class IVChar:
    def __init__(self, V, I, mask=None, **kw):
        self.V, self.I = V, I
        self.fitter_IV = FitterIV(self.V.x, self.I.x, **kw)

    def fit(self, out=None):
        if out is None:
            out = np.empty(3)
        try:
            out[:] = self.fitter_IV.fit()
        except FitterError:
            out[:] = np.nan
        return out
    
    def mask(self, out=None):
        if out is None:
            out = np.empty(self.I.size, bool)
        out.fill(True)
        if self.fitter_IV.p is not None:
            out[self.fitter_IV.get_ind()] = False
        return out

    def plot(self, ax=None, fun='get_xy', lines=None):
        return self.fitter_IV.plot(ax, fun, lines)


class IVGroup:
    def __init__(self, V, II, s=slice(None), **kw):
        self.IV_char = np.empty(len(II), object)
        for j, I in enumerate(II):
            self.IV_char[j] = IVChar(V[s], I[s], **kw)

    def __getitem__(self, index):
        return self.IV_char[index]

    def fit(self, out=None):
        if out is None:
            out = np.empty((self.IV_char.size, 3))
        for p, IV_char in zip(out, self.IV_char):
            IV_char.fit(p)
        return out

    def mask(self, out=None):
        if out is None:
            out = np.empty((self.IV_char.size, self.IV_char[0].I.size), bool)
        for I, IV_char in zip(out, self.IV_char):
            IV_char.mask(I)
        return out

    def plot(self, ax=None, fun='get_xy'):
        ax = get_axes(ax)
        cache = getattr(ax, 'lines_cache', None)
        if cache is None:
            ax.lines_cache = [x.plot(ax, fun) for x in self.IV_char]
        else:
            ax.lines_cache = [x.plot(ax, fun, l) for x, l in zip(self.IV_char, cache)]


class IVSeriesViewer(ToggleViewer):
    def __init__(self, IV_series):
        self.IV_series = IV_series

        ToggleViewer.__init__(self, 'IV viewer')

    def plotfun(self, event):
        t_event = event.xdata
        ti = self.IV_series.ti[:,1]
        i = min(np.searchsorted(ti, t_event), ti.size-1)
        self.IV_series.IV_group[i].plot(self.ax, 'get_xy')
        return self.ax.lines

    def viewer(self, event):
        fig = get_tfig(figsize=(6,5), xlab="V [V]")
        self.ax = fig.axes[0]
        self.ax.set_ylabel("I [A]")
        self.ax.set_xlim(self.IV_series.V_range)
        self.ax.set_ylim(self.IV_series.I_range)


class IVSeries:
    def __init__(self, V, II, iE, **kw):       
        self.V, self.II, self.iE = V, II, iE
        self.ti = V.t[iE]

        self.V_range = V.plot_range()
        self.I_range = self._plot_range(II)

        N = iE.shape[0]
        self.siz = (N, len(II))

        self.IV_group = np.empty(N, object)
        for j in xrange(N):
            s = self._slice(j)
            self.IV_group[j] = IVGroup(V, II, s, **kw)

    def __getitem__(self, index):
        return self.IV_group[index]

    def _plot_range(self, II):
        I_range = np.array([I.plot_range() for I in II])
        return np.nanmin(I_range[:,0]), np.nanmax(I_range[:,1])

    def _slice(self, j):
        return slice(self.iE[j,0], self.iE[j,1]+1)

    @memoized_property
    def PP(self):
        print "Calculating PP..."
        out = np.empty(self.siz + (3,))
        for p, IV_group in zip(out, self.IV_group):
            IV_group.fit(p)

        i0 = np.r_[self.iE[:,0], self.iE[-1,1]]
        return PiecewisePolynomial(out[None], self.V.t, i0=i0)
    
    def set_PP(self, PP):
        self.PP = PP

    def fit(self):
        del self.PP
        return self.PP

    def mask(self):
        out = np.empty((len(self.II), self.II[0].size), bool)
        for j, IV_group in enumerate(self.IV_group):
            s = self._slice(j)
            IV_group.mask(out[:,s])
        return out

    def eval(self, V_mask=False, I_mask=True):
        t, x = self.V.t, self.V.x
        
        if isinstance(V_mask, np.ndarray):
            x = ma.masked_array(x, V_mask)
       
        out = FitterIV.fitfun(self.PP(t).T, x)
        
        if I_mask is True:
            I_mask = self.mask()
        if isinstance(I_mask, np.ndarray):
            out = ma.masked_array(out, I_mask)
        return out

    @memoized_property
    def II_diff(self):
        II_fit = self.eval(I_mask=False)
        return [I - I_fit for I, I_fit in zip(self.II, II_fit)]

    def plot(self, fig=None, **kw):
        n = len(self.II)
        IIfit = self.eval(**kw)

        fig = get_fig(fig, (n, 1), xlab="t (s)")

        for ax, I, Ifit in zip(fig.axes, self.II, IIfit):
            ax.plot(I.t, I.x, I.t, Ifit)
            ax.set_ylabel(I.name)

        return fig
        
    def animate(self, fun='get_xy'):
        ion()
        fig = figure()
        ax = fig.gca()
        if fun == 'get_xy':
            ax.set_xlim(self.V_range)
            ax.set_ylim(self.I_range)
        else:
            ax.set_xlim(( 0.0, 1.0))
            ax.set_ylim((-1.2, 1.2))
        
        canvas = fig.canvas
        canvas.draw()
        background = canvas.copy_from_bbox(ax.bbox)

        for IV_group in self.IV_group:
            canvas.restore_region(background)
            IV_group.plot(ax, fun)
            for line in ax.lines:
                ax.draw_artist(line)
            canvas.blit(ax.bbox)


class PhysicalResults:
    def __init__(self, shn, R, i0, meas, usetex=usetex):
        self.shn, self.R, self.i0, self.meas, self.usetex = shn, R, i0, meas, usetex

        self.keys = ('n_cs', 'Mach', 'nv', 'mnv', 'j', 'Vf', 'Te', 'Vp', 
                'cs', 'n', 'v', 'pe', 'R', 't', 'Dt')

        sup = lambda x: r'$^{\mathdefault{%s}}$' % x
        sub = lambda x: r'$_{\mathdefault{%s}}$' % x

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
        self.lim['n'] = (0, 3e19)
        self.lim['Mach'] = (-2, 2)
        self.lim['Te'] = (0, 100)
        self.lim['R'] = (0, None)

    @memoized_property
    def PP(self):
        qe = 1.6022e-19
        mi = 2*1.67e-27

        Gp = self.meas.jp/qe
        Gm = self.meas.jm/qe

        Gp[Gp <= 0] = np.nan
        Gm[Gm <= 0] = np.nan

        dtype = zip(self.keys, [np.double]*len(self.keys))
        res = np.empty(Gp.size, dtype).view(np.recarray)

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

        return PiecewisePolynomial(res[None], self.R.t, i0=self.i0)

    def eval(self, plunge=None, inout=None):
        w = self.R.plunges(plunge, inout)
        
        i, y = self.PP.eval(w=w)
        
        R0 = 1.645
        y.t  = self.R.t[i]
        y.R  = R0 - self.R.x[i]
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

    def plot_R(self, **kw):
        return self.plot(xkey='R', **kw)

    def plot_R_in(self, **kw):
        return self.plot_R(inout='in', **kw)

    def plot_R_out(self, **kw):
        return self.plot_R(inout='out', **kw)


class ResultsIOError(Exception):
    pass

class ResultsIOWarning(Warning):
    pass

class Probe:
    def __init__(self, digitizer=None):
        self.digitizer = digitizer
        self.PP = None

        self.xlab = "t (s)"
        self.ylab = ("Isat (A)", "Vf (V)", "Te (eV)")

    def __getitem__(self, index):
        if not isinstance(index, tuple):
            return self.S[index]
        else:
            return tuple([self.S[i] for i in index])

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
        is_type = lambda x: x.type == type
        return np.array(filter(is_type, self.S.itervalues()))

    def _shortcut_factory(name):
        @memoized_property
        def wrapper(self):
            return self.get_type(name)
        return wrapper

    I = _shortcut_factory('Current')
    V = _shortcut_factory('Voltage')
    R = _shortcut_factory('Position')

    @memoized_property
    def is_swept(self):
        is_swept = lambda I: I.V.is_swept
        return np.array(map(is_swept, self.I))

    @memoized_property
    def I_swept(self):
        return self.I[self.is_swept]

    def plot_raw(self, fig=None,
            keys = (('Position',), ('Current',), ('Voltage',))):
        keys = np.array(keys, ndmin=2)

        fig = get_fig(fig, keys.shape, xlab=self.xlab, ylab=keys)

        ax = fig.axes[np.flatnonzero(keys == 'Voltage')]
        if ax:
            for I in self.I:
                I.V.plot(ax)

        for key, ax in zip(keys, fig.axes):
            for S in self.get_type(key[0]):
                S.plot(ax)

        fig.canvas.draw()
        return fig
    
    def get_dwell_params(self):
        R = self.R[0]
        iM = R.t_ind[1]
        tM = R.t[iM]
        RM = R.x[iM]
        return tM, RM

    def trim(self, plunge='all'):
        R = self.R[0]
        i0, iM, i1 = R.t_ind

        if plunge == 'all':
            s = np.concatenate(map(np.arange, i0, i1))
        else:
            s = slice(i0[plunge], i1[plunge])

        for S in self.S.itervalues():
            S.trim(s)

    def smooth_I(self, w=10):
        for I in self.I:
            I.smooth(w)

    def corr_capa(self):
        for I in self.I:
            I.x[:] -= I.I_capa()

    def calc_IV_series(self, n=1, **kw):
        II = self.I
        V = self.I_swept[0].V
        iE = np.c_[V.iE[:-n], V.iE[n:]]
        return IVSeries(V, II, iE, **kw)

    @memoized_property
    def IV_series(self):
        IV_series = self.calc_IV_series()
        if self.PP is not None:
            IV_series.set_PP(self.PP)
        return IV_series

    def analyze(self):
        self.PP = self.IV_series.PP

    def calc_IV_series_simple(self):
        return IVSeriesSimpleGroup(self.I_swept, self.S['Rs'])

    @memoized_property
    def IV_series_simple(self):
        return self.calc_IV_series_simple()

    def analyze_simple(self):
        self.PP_simple = [IV.fit() for IV in self.IV_series_simple]
    
    def PP_fluct(self, Vmax=-150):
        self.analyze()
        II_diff = self.IV_series.II_diff
        i0 = self.PP.i0

        out = np.empty((len(i0)-1, len(II_diff), 3))
        for i, I_diff in enumerate(II_diff):
            out[:,i,0] = I_diff.Isat(Vmax).apply_fun(np.std, i0[:-1], i0[1:])
            
        out[:,:,1] = self.PP.c[0,:,:,0]
        out[:,:,2] = out[:,:,0]/out[:,:,1]
        return PiecewisePolynomial(out[None], self.PP.x, i0=i0)

    def get_meas(self):
        pass

    @memoized_property
    def res(self):
        if self.PP is None:
            try:
                self.load()
            except ResultsIOError:
                self.analyze()
                self.save_res()

        Isat = self.PP.c[0,:,:,0].T
        Vf   = self.PP.c[0,:,:,1].T
        Te   = self.PP.c[0,:,:,2].T

        keys = ('jp', 'jm', 'Vf', 'Te')
        dtype = zip(keys, [np.double]*len(keys))
        meas = np.empty(Isat.shape[1], dtype).view(np.recarray)

        self.get_meas(Isat, Vf, Te, meas)
       
        shn = self.digitizer.shn
        return PhysicalResults(shn, self['Rs'], self.PP.i0, meas)
        
    @memoized_property
    def h5name_res(self):
        return self.digitizer.IO_file.h5name[:-3] + "_res.h5"

    def save_res(self):
        IO = IOH5(self.h5name_res)
        d = self.PP.savefields
        try:
            IO.save(d)
        except IOError:
            warnings.warn("could not save results", ResultsIOWarning)

    def load_res(self):
        IO = IOH5(self.h5name_res)
        try:
            d = IO.load()
        except IOError:
            raise ResultsIOError("no results available")
        self.PP = PiecewisePolynomial(**d)

    def load(self, **kw):
        self.load_raw(**kw)
        self.load_res()

    def analyze2(self, n=2):
        self.IV_series2 = self.calc_IV_series(n=2)
        self.PP2 = self.IV_series2.fit()

        """
        II = self.get_type('Current')
        V = II[0].V
        PP = self.PP(V.t)
        mask = self.IV_series.mask()
        
        iE = self.IV_series.iE
        self.IV_series2 = IVSeries2(V, II, PP, mask, iE)
        """

    def plot(self, fig=None, PP='PP', x=None, plunge=None, inout=None):
        if self.PP is None:
            self.analyze()
        PP = getattr(self, PP)

        ylab = self.ylab

        if x is None:
            xlab = self.xlab
        elif x == 'R':
            x = 100*self['Rs'].x
            xlab = "R [cm]"

        w = self['Rs'].plunges(plunge, inout)

        if fig is None:
            IV_series_viewer = IVSeriesViewer(self.IV_series)
            fig = IV_series_viewer.get_fig(shape=(3,1), figsize=(10,10),
                                           xlab=xlab, ylab=ylab)

        ax = fig.axes
        for i, pp in enumerate(PP.T):
            pp.plot(ax[i], x=x, w=w)

        fig.canvas.draw()
        return fig

    def plot_R(self, **kw):
        return self.plot(x='R', **kw)


