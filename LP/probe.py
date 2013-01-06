import numpy as np

import scipy.interpolate as interp
import scipy.optimize as opt
import scipy.odr as odr

import warnings

import logging
reload(logging)
logging.basicConfig(level=logging.WARN)
logger = logging

from pdb import set_trace

from sig import *

import LP.fitfun

class ArrayView(np.ndarray):
    def __new__(subtype, x, fields):
        dtype = {f: x.dtype.fields[f] for f in fields}
        return np.ndarray.__new__(subtype, x.shape, dtype,
                                  buffer=x, strides=x.strides)


class MouseMotionObserverViewer:
    def __init__(self, observers, viewer, plotfun):
        self.observers, self.viewer = observers, viewer
        self.plotfun = plotfun

        self.observer_canvas = self.observers[0].figure.canvas
        self.viewer_canvas = self.viewer.figure.canvas
        
        self.cid = self.observer_canvas.mpl_connect('motion_notify_event', self.on_move)
        self.viewer_canvas.mpl_connect('resize_event', self.on_resize)

        self.on_resize(None)

    def on_resize(self, event):
        logger.debug("Resizing")
        save_vis = [line.get_visible() for line in self.viewer.lines]
        for line in self.viewer.lines:
            line.set_visible(False)
        
        self.viewer_canvas.draw()
        self.background = self.viewer_canvas.copy_from_bbox(self.viewer.bbox)
        
        for line, vis in zip(self.viewer.lines, save_vis):
            line.set_visible(vis)
            self.viewer.draw_artist(line)

    def on_move(self, event):
        logger.debug("Moving")
        if event.inaxes in self.observers:
            self.viewer_canvas.restore_region(self.background)
                
            self.plotfun(event)

            for line in self.viewer.lines:
                self.viewer.draw_artist(line)
            self.viewer_canvas.blit(self.viewer.bbox)

    def __del__(self):
        logger.debug("Destroying")
        self.observer_canvas.mpl_disconnect(self.cid)


class FitterError(Exception):
    pass

class Fitter:
    def __init__(self, x, y, args=(), engine='fmin'):
        self.x, self.y, self.args = x, y, args

        self.OK = None
        self.X = self.Y = None
        self.P = self.P0 = None
        self.p = self.p0 = None

        if hasattr(self, 'fitfun_diff'):
            self.f = self.fitfun_diff
            self.wrap_fmin = self.wrap_fmin_diff
        elif hasattr(self, 'fitfun_fast'):
            self.f = self.fitfun_fast
            self.wrap_fmin = self.wrap_fmin_prealloc
        else:
            self.f = self.fitfun
            self.wrap_fmin = self.wrap_fmin_noprealloc

        self.engines = dict(fmin=self.wrap_fmin, odr=self.wrap_odr)
        self.set_engine(engine)

    # overload
    def set_OK(self):
        self.OK = True
        return self.OK

    def set_norm(self):
        self.X, self.Y = self.x, self.y
        return self.X, self.Y

    def set_unnorm(self):
        self.p, self.p0 = self.P, self.P0
        return self.p, self.p0

    def set_guess(self):
        self.P0 = 0.
        return self.P0

    @staticmethod
    def fitfun(P, X):
        pass

    # static
    @staticmethod
    def wrap_fmin_noprealloc(fun, p0, x, y, *args):
        def dy2(p):
            dy = fun(p, x, *args) - y
            return dy.dot(dy)/dy.size

        return opt.fmin(dy2, p0, disp=False)

    @staticmethod
    def wrap_fmin_prealloc(fun, p0, x, y, *args):
        out = np.empty_like(x)

        def dy2(p):
            fun(p, x, out, *args)
            dy = out - y
            return dy.dot(dy)/dy.size

        return opt.fmin(dy2, p0, disp=False)

    @staticmethod
    def wrap_fmin_diff(fun, p0, *args):
        return opt.fmin(fun, p0, args=args, disp=False)

    @staticmethod
    def wrap_odr(fun, p0, x, y):
        mod = odr.Model(fun)
        dat = odr.Data(x, y)
        o = odr.ODR(dat, mod, p0)
        out = o.run()
        return out.beta

    def is_OK(self):
        if self.OK is None:
            self.set_OK()
        return self.OK

    def get_norm(self):
        if self.X is None:
            self.set_norm()
        return self.X, self.Y

    def get_unnorm(self):
        if self.p is None:
            self.set_unnorm()
        return self.p, self.p0

    def get_guess(self):
        if self.P0 is None:
            self.set_guess()
        return self.P0
    
    def set_engine(self, engine):
        self.engine = self.engines[engine]

    def fit(self, P0=None):
        if not self.is_OK():
            raise FitterError("Cannot fit data that failed is_OK() check")
        
        if P0 is None:
            P0 = self.get_guess()
        X, Y = self.get_norm()
        self.P = self.engine(self.f, P0, X, Y, *self.args)
        self.set_unnorm()

        return self.p

    def eval_guess_norm(self, X):
        return self.fitfun(self.P0, X, *self.args)

    def eval_guess(self, x):
        return self.fitfun(self.p0, x, *self.args)

    def eval_norm(self, X):
        return self.fitfun(self.P, X, *self.args)

    def eval(self, x):
        return self.fitfun(self.p, x, *self.args)

    def get_XY(self):
        if self.is_OK(): 
            if self.P is None:
                self.fit()
            return self.X, (self.Y, self.eval_norm(self.X))
        else:
            return self.X, (self.Y,)

    def get_xy(self):
        if self.is_OK(): 
            if self.p is None:
                self.fit()
            return self.x, (self.y, self.eval(self.x))
        else:
            return self.x, (self.y,)

    def plot(self, ax=None, fun='get_xy', lines=None):
        x, y = getattr(self, fun)()
        sty = (dict(linewidth=1.5),)
       
        ax = get_axes(ax)

        if lines is None: lines = []
        Nl, Ny = len(lines), len(y)

        for li, yi in zip(lines, y):
            li.set_data(x, yi)
            li.set_visible(True)

        for li in lines[Ny:]:
            li.set_visible(False)

        if Nl > 0:
            color = lines[0].get_color()
        else:
            color = ax._get_lines.color_cycle.next()
            sty = (dict(),) + sty

        for yi, styi in zip(y[Nl:], sty):
            li, = ax.plot(x, yi, color=color, **styi)
            lines.append(li)

        return lines


class FitterIV(Fitter):
    def __init__(self, V, I, mask=None, cut_at_min=True, **kw):
        self.mask_ind = np.arange(V.size)
        if mask is not None:
            self.mask_ind = self.mask_ind[mask]

        self.sort_ind = V[self.mask_ind].argsort()
        self.ind = self.mask_ind[self.sort_ind]

        self.V, self.I = V[self.ind], I[self.ind]

        self.im = self.I.argmin()
        self.Vm, self.VM = self.V[0], self.V[-1]
        self.Im, self.IM = self.I[self.im], np.median(self.I[:self.I.size/2])
        self.dV = self.VM - self.Vm
        self.dI = self.IM - self.Im

        self.cut_at_min = cut_at_min
        self.M = self.V.size

        Fitter.__init__(self, self.V, self.I, **kw)

    def get_ind(self):
        return self.ind[:self.M]

    def set_OK(self):
        def medianstd(x):
            xm = np.median(x)
            dx = x-xm
            xs = np.sqrt(dx.dot(dx)/dx.size)
            return xm, xs

        N = self.I.size
        Ilm, Ils = medianstd(self.I[:N/10])
        Irm, Irs = medianstd(self.I[N*9/10:])

        cnd1 = Ilm > Ils
        cnd2 = Irm < -2*Ils
        self.OK = cnd1 & cnd2
        return self.OK

    def set_norm(self):
        if self.cut_at_min:
            self.M = self.im+1
        self.X = (self.V[:self.M].astype('d') - self.Vm)/self.dV
        self.Y = (self.I[:self.M].astype('d') - self.Im)/self.dI*2 - 1
        return self.X, self.Y

    def set_unnorm(self):
        self.p0 = self.LP_unnormalize(self.P0, self.Vm, self.VM, self.Im, self.IM)
        self.p  = self.LP_unnormalize(self.P , self.Vm, self.VM, self.Im, self.IM)
        return self.p, self.p0

    def set_guess(self):
        def find_i0(x):
            return np.flatnonzero(np.diff(np.sign(x)))[-1]

        V, I = self.get_norm()
        I0 = 1.
        i0 = find_i0(I)
        Vf = (V[i0]+V[i0+1])/2
        Te = (V[self.im]-Vf)/np.log(1-I[self.im]/I0)

        self.P0 = np.array([I0, Vf, Te])
        return self.P0

    @staticmethod
    def LP_unnormalize(P, Vm, VM, Im, IM):
        dV, dI = VM-Vm, IM-Im
        a = 2./dI
        b = -(Im*a+1)

        p = np.empty_like(P)
        p[0] = (P[0]-b)/a
        p[1] = (P[1]+P[2]*np.log(1.-b/P[0]))*dV+Vm
        p[2] = P[2]*dV
        return p

    @staticmethod
    def LP_normalize(p, Vm, VM, Im, IM):
        dV, dI = VM-Vm, IM-Im
        a = 2./dI
        b = -(Im*a+1)

        P = np.empty_like(p)
        P[0] = a*p[0]+b
        P[1] = (p[1]-p[2]*np.log(1.-b/P[0])-Vm)/dV
        P[2] = p[2]/dV
        return P

    @staticmethod
    def fitfun(P, X):
        iP2 = 1./P[2]
        return P[0]*(1.-np.exp((X-P[1])*iP2))

    #fitfun_fast = LP.fitfun.IV3
    fitfun_diff = LP.fitfun.IV3_diff

    def fit(self):
        Fitter.fit(self)
        save = self.X, self.Y
        Y0 = 0.
        while True:
            self.M *= 0.95
            self.X, self.Y = self.X[:self.M], self.Y[:self.M]
            P_old = self.P
            Fitter.fit(self, P0=self.P)
            
            if np.any(self.P > P_old) or (self.eval_norm(self.X[-1]) > Y0):
                self.P = P_old
                break

        self.X, self.Y = save
        self.set_unnorm()
        self.check_Te()
        return self.p

    def check_Te(self):
        n, Vf, Te = self.p
        if n < 0.:
            raise FitterError("Negative n")
        if Te < 0. or Te > 0.5*self.dV:
            raise FitterError("Unrealistic Te")


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
            

class IVSeriesViewer:
    def __init__(self, IV_series):
        self.IV_series = IV_series
        self.MMOV = None

    def on_close(self, event):
        logger.debug("Closing")
        self.MMOV = None

    def plotfun(self, event):
        t_event = event.xdata
        ti = self.IV_series.ti[:,1]
        i = min(np.searchsorted(ti, t_event), ti.size-1)
        self.IV_series.IV_group[i].plot(self.ax, 'get_xy')

    def toggle(self, event):
        fig = event.inaxes.figure
        if self.MMOV is not None:
            self.MMOV.viewer_canvas.manager.destroy()
            return

        fig2 = get_tfig(figsize=(6,5), xlab="V [V]")
        self.ax, = fig2.axes
        self.ax.set_ylabel("I [A]")
        self.ax.set_xlim(self.IV_series.V_range)
        self.ax.set_ylim(self.IV_series.I_range)

        fig2.canvas.mpl_connect('close_event', self.on_close)

        self.MMOV = MouseMotionObserverViewer(fig.axes, self.ax, self.plotfun)
        fig2.show()


class FitterIV2(Fitter):
    def __init__(self, V, I, a, p0):
        Fitter.__init__(self, V, I, args=(a,))

        self.dV = dV = V.ptp()
        self.dI = dI = I.ptp()
        self.fact = np.array((dI, dV, dV, dI, dV, dV))

        self.p0 = p0

    def set_OK(self):
        self.OK = np.isfinite(self.p0).all()
        return self.OK

    def set_norm(self):
        self.X = self.x / self.dV
        self.Y = self.y / self.dI
        return self.X, self.Y

    def set_unnorm(self):
        self.p = self.P * self.fact
        return self.p, self.p0

    def set_guess(self):
        self.P0 = self.p0 / self.fact
        return self.P0

    @staticmethod
    def fitfun(p, V, a):
        Is = p[3] + a*(p[0]-p[3])
        Vf = p[4] + a*(p[1]-p[4])
        Te = p[5] + a*(p[2]-p[5])
        return Is*(1.-np.exp((V-Vf)/Te))

    fitfun_diff = LP.fitfun.IV6_diff


class IVSeries2:
    def __init__(self, V, II, PP, mask, iE=None, **kw):
        self.V, self.II, self.PP, self.mask, self.iE = V, II, PP, mask, iE

    def plot(self, i, j):
        s = slice(self.iE[i], self.iE[j]+1)
        V = self.V[s]
        I = self.II[0][s]
        pp = self.PP[s, 0]
        m = self.mask[0, s]

        t = V.t[m]
        V = V.x[m]
        I = I.x[m]
        pm = pp[m]

        print pp.shape

        def fitfun(p, V): 
            return p[0]*(1.-np.exp((V-p[1])/p[2]))

        a = (t-t[0])/(t[-1]-t[0])

        def fitfun2(p, V):
            Is = p[0] + a*(p[3]-p[0])
            Vf = p[1] + a*(p[4]-p[1])
            Te = p[2] + a*(p[5]-p[2])
            return Is*(1.-np.exp((V-Vf)/Te))

        Ifit = fitfun(pm.T, V)

        """
        p0 = pm.mean(0)
        p = Fitter.wrap_fmin(fitfun, p0, V, I)
        Ifit2 = fitfun(p, V)
        """

        p0 = np.r_[pm[0], pm[-1]]
        p = Fitter.wrap_fmin(fitfun2, p0, V, I)
        Ifit2 = fitfun2(p, V)

        ax = get_axes()
        ax.plot(t, I, t, Ifit, t, Ifit2)
        return p


class IVSeriesSimple:
    def __init__(self, S):
        self.S = S

    def fit(self, n=1):
        iE = self.S.V.iE
        i0, i1 = iE[:-n], iE[n:]
        N = len(i0)
        out = np.empty((N, 3))
        out.fill(np.nan)
        self.mask = np.zeros(self.S.size, bool)

        for j in xrange(N):
            s = slice(i0[j], i1[j])
            S = self.S[s]
            fitter_IV = FitterIV(S.V.x, S.x)
            try:
                out[j] = fitter_IV.fit()
                self.mask[s][fitter_IV.get_ind()] = True
            except FitterError:
                pass
        self.PP = PiecewisePolynomial(out[None], self.S.t, i0=i0, i1=i1)
        return self.PP

    def fit2(self, n=5):
        Sfit = self.get_Sfit()
        #mask = ~Sfit.x.mask

        iE = self.S.V.iE
        i0, i1 = iE[:-n], iE[n:]

        N = len(i0)
        out = np.empty((N, 6))
        out.fill(np.nan)

        p = np.concatenate((self.PP.c[0,n-1:], self.PP.c[0,:N]), axis=1)
        
        t0, t1 = self.S.t[i0], self.S.t[i1]
        dt = t1 - t0

        for j in xrange(N):
            s = slice(i0[j], i1[j])
            mask = self.mask[s]
            if not np.any(mask):
                continue

            S = self.S[s]
            S = S[mask]

            a = (S.t - t0[j]) / dt[j]

            fitter_IV = FitterIV2(S.V.x, S.x, a, p[j])
            try:
                out[j] = fitter_IV.fit()
            except FitterError:
                pass

        c = np.swapaxes(out.reshape((N, 2, 3)), 0, 1).copy()
        c[0] -= c[1]
        c[0] /= dt[:, None]

        self.PP2 = PiecewisePolynomial(c, self.S.t, i0=i0, i1=i1)
        return self.PP2

    def get_Sfit(self):
        t, V = self.S.t, self.S.V
        p = self.PP(t).T
        Ifit = FitterIV.fitfun(p, V.x)
        Ifit_masked = ma.masked_array(Ifit, V.x > p[1] + p[2])
        return CurrentSignal(Ifit_masked, t, V=V)

    def plot(self, ax=None):
        ax = get_axes(ax)
        self.S.plot(ax=ax)
        self.get_Sfit().plot(ax=ax)
        return ax


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

        self.xlab = "t [s]"
        self.ylab = ("Isat [A]", "Vf [V]", "Te [eV]")

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
        return map(IVSeriesSimple, self.I_swept)

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
            menu_entries_ax = (('IV viewer', IV_series_viewer.toggle),)

            fig = get_fig(shape=(3,1), xlab=xlab, 
                          figsize=(10,10), menu_entries_ax=menu_entries_ax)

            fig.IV_series_viewer = IV_series_viewer

        ax = fig.axes
        for i, pp in enumerate(PP.T):
            pp.plot(ax[i], x=x, w=w)
            ax[i].set_ylabel(ylab[i])

        fig.canvas.draw()
        return fig

    def plot_R(self, **kw):
        return self.plot(x='R', **kw)


