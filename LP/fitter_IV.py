import numpy as np
import numpy.ma as ma

from itertools import cycle

from sig import median, memoized_property, get_fig, get_tfig
from sig import DictView, Container, CurrentSignal, PiecewisePolynomialEndpoints

from fitter import Fitter, FitterError

try:
    import LP.fitfun as ff
except ImportError:
    import fitfun as ff

import LP.mag_fit as mag_fit

from sm_pyplot.contextmenupicker import ContextMenuPicker
from sm_pyplot.observer_viewer import ToggleViewer, ToggleViewerIntegrated


class FitterIV(Fitter):
    def __init__(self, V, I, mask=None, **kw):
        self.mask_ind = np.arange(V.size)
        if mask is not None:
            self.mask_ind = self.mask_ind[mask]

        self.sort_ind = V[self.mask_ind].argsort()
        self.ind = self.mask_ind[self.sort_ind]

        self.V, self.I = V[self.ind], I[self.ind]
        
        self.im = self.I.argmin()
        self.Vm, self.VM = self.V[0], self.V[-1]
        self.Im, self.IM = self.I[self.im], median(self.I[:self.I.size/2])
        self.dV = self.VM - self.Vm
        self.dI = self.IM - self.Im

        self.M = self.V.size

        self.cut_at_min = kw.pop('cut_at_min', True)
        self.r = kw.pop('r', 0.95)
        
        Fitter.__init__(self, self.V, self.I, **kw)

    def get_ind(self):
        return self.ind[:self.M]

    def set_OK(self):
        def medianstd(x):
            xm = median(x)
            dx = x-xm
            xs = np.sqrt(dx.dot(dx)/dx.size)
            return xm, xs

        N = self.I.size
        Ilm, Ils = medianstd(self.I[:N/10])
        Irm, Irs = medianstd(self.I[N*9/10:])

        cnd1 = Ilm > Ils
        cnd2 = Irm < -2*Ils
        self.OK = cnd1 & cnd2

    def set_norm(self):
        if self.cut_at_min:
            self.M = self.im+1
        self.X = (self.V[:self.M].astype('d') - self.Vm)/self.dV
        self.Y = (self.I[:self.M].astype('d') - self.Im)/self.dI*2 - 1

    def set_unnorm(self):
        '''self.p0 commented out for now (should only be calculated on-demand)'''
        #self.p0 = self.LP_unnormalize(self.P0)
        self.p  = self.LP_unnormalize(self.P )

    def set_guess(self):
        def zero_crossings(x):
            return np.flatnonzero(np.diff(np.sign(x)))

        V, I = self.get_norm()
        I0 = 1.
        i0 = zero_crossings(I)
        Vf = np.mean((V[i0] + V[i0+1])/2)
        Te = (V[self.im]-Vf)/np.log(1-I[self.im]/I0)

        self.P0 = np.array([I0, Vf, Te])

    def LP_unnormalize(self, P):
        a = 2./self.dI
        b = -(self.Im*a+1)

        p = np.empty_like(P)
        p[0] = (P[0]-b)/a
        p[1] = (P[1]+P[2]*np.log(1.-b/P[0]))*self.dV+self.Vm
        p[2] = P[2]*self.dV
        return p

    def LP_normalize(self, p):
        a = 2./self.dI
        b = -(self.Im*a+1)

        P = np.empty_like(p)
        P[0] = a*p[0]+b
        P[1] = (p[1]-p[2]*np.log(1.-b/P[0])-self.Vm)/self.dV
        P[2] = p[2]/self.dV
        return P

    @staticmethod
    def fitfun(P, X):
        iP2 = 1./P[2]
        return P[0]*(1.-np.exp((X-P[1])*iP2))

    custom_engine = ff.IV3_fit
    fitfun_fast = ff.IV3
    fitfun_diff = ff.IV3_diff
    fitfun_rms = ff.IV3_rms

    def fit(self):
        Fitter.fit(self)

        if self.r < 1:
            save = self.X, self.Y
            Y0 = 0.
            while True:
                self.M *= self.r
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


class FitterIVMag(FitterIV):
    def __init__(self, V, I, **kw):
        FitterIV.__init__(self, V, I, **kw)

        self.Ifit = np.zeros_like(self.I)

        self.c_params = np.array([5e18, 18., 0., 12., 18., 1e-8, 
            np.pi/2, 0., 0., 0., 0.0009, -1., 2., 1., 2., 0., 1., 0.0016/2, 0., 0.])

        self.do_var = np.array([1, 1, 0, 1, 1, 1, 
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 'i')

    def fit(self):
        c_params = self.c_params.copy()

        mag_fit.magfit(self.V, -self.I, self.Ifit, c_params, self.do_var)

        #Ifit = np.zeros_like(self.V)
        #mag_fit.mag_doppel(self.V, Ifit, c_params)

        #print Ifit - self.Ifit

        np.negative(self.Ifit, self.Ifit)
        return c_params[[0, 4, 1]]


class FitterIV6(Fitter):
    def __init__(self, V, I, a, p0, **kw):
        Fitter.__init__(self, V, I, args=(a,), **kw)

        self.dV = dV = V.ptp()
        self.dI = dI = I.ptp()
        self.fact = np.array((dI, dV, dV, dI, dV, dV))

        self.p0 = p0

    def set_OK(self):
        self.OK = np.isfinite(self.p0).all()

    def set_norm(self):
        self.X = self.x / self.dV
        self.Y = self.y / self.dI

    def set_unnorm(self):
        self.p = self.P * self.fact

    def set_guess(self):
        self.P0 = self.p0 / self.fact

    @staticmethod
    def fitfun(p, V, a):
        Is = p[0] + a*(p[3]-p[0])
        Vf = p[1] + a*(p[4]-p[1])
        Te = p[2] + a*(p[5]-p[2])
        return Is*(1.-np.exp((V-Vf)/Te))

    custom_engine = ff.IV6_fit
    fitfun_fast = ff.IV6
    fitfun_diff = ff.IV6_diff
    fitfun_rms = ff.IV6_rms


class FitterIV6i(FitterIV6):
    @staticmethod
    def fitfun(p, V, a):
        Is = p[0] + a*(p[3]-p[0])
        Vf = p[1] + a*(p[4]-p[1])
        iTe = 1./p[2] + a*(1./p[5]-1./p[2])
        return Is*(1.-np.exp((V-Vf)*iTe))

    custom_engine = ff.IV6i_fit
    fitfun_fast = ff.IV6i
    fitfun_diff = ff.IV6i_diff
    fitfun_rms = ff.IV6i_rms


class FitterIV6Perm(FitterIV6):
    def __init__(self, *args, **kw):
        # derived class needs to assign self.perm and self.iperm here

        FitterIV6.__init__(self, *args, **kw)
        self.fact = self.fact[self.perm]

    def set_unnorm(self):
        self.p = (self.P * self.fact)[self.iperm]
        return self.p, self.p0

    def set_guess(self):
        self.P0 = self.p0[self.perm] / self.fact


class FitterIV5(FitterIV6Perm):
    def __init__(self, *args, **kw):
        self.perm = np.array((0,1,2,3,4))
        self.iperm = np.array((0,1,2,3,4,2))
        FitterIV6Perm.__init__(self, *args, **kw)

    custom_engine = ff.IV5_fit
    fitfun_fast = ff.IV5
    fitfun_diff = ff.IV5_diff
    fitfun_rms = ff.IV5_rms


class FitterIV4(FitterIV6Perm):
    def __init__(self, *args, **kw):
        self.perm = np.array((0,1,2,3))
        self.iperm = np.array((0,1,2,3,1,2))
        FitterIV6Perm.__init__(self, *args, **kw)

    custom_engine = ff.IV4_fit
    fitfun_fast = ff.IV4
    fitfun_diff = ff.IV4_diff
    fitfun_rms = ff.IV4_rms


class IVViewer(ToggleViewer):
    def __init__(self, IV, PP='PP'):
        self.IV, self.PP = IV, PP

        ToggleViewer.__init__(self, 'IV viewer')

    def plotfun(self, event):
        t_event = event.xdata
        res = self.IV.get_Sfit_at_event(t_event, PP=self.PP)
        if not isinstance(res, list):
            res = [res]

        self.clear()
        lines = []
        colorcycle = cycle(self.colors)
        for V, I, Ifit, t in res:
            color = colorcycle.next()
            lines += self.ax.plot(V, I, color=color)
            lines += self.ax.plot(V, Ifit, color=color, linewidth=1.0)
        return lines

    def viewer(self, event):
        fig = get_tfig(figsize=(6,5), xlab="V (V)")
        self.ax = fig.axes[0]
        self.ax.set_ylabel("I (A)")
        self.ax.set_xlim(self.IV.plot_range_V())
        self.ax.set_ylim(self.IV.plot_range_I())
        self.colors = ('b', 'g', 'r', 'c')


class IVViewerIt(ToggleViewer):
    def __init__(self, IV, PP='PP'):
        self.IV, self.PP = IV, PP

        ToggleViewer.__init__(self, 'I(t) viewer')

    def plotfun(self, event):
        t_event = event.xdata
        res = self.IV.get_Sfit_at_event(t_event, PP=self.PP)
        if not isinstance(res, list):
            res = [res]
        
        self.clear()
        lines = []
        colorcycle = cycle(self.colors)
        for V, I, Ifit, t in res:
            dt = t - t[0]
            color = colorcycle.next()
            lines += self.ax.plot(dt, I, color=color)
            lines += self.ax.plot(dt, Ifit, color=color, linewidth=1.5)
        return lines

    def viewer(self, event):
        fig = get_tfig(figsize=(6,5), xlab="$\Delta$t (s)")
        self.ax = fig.axes[0]
        self.ax.set_ylabel("I (A)")

        self.ax.set_xlim(self.IV.plot_range_dt(self.PP))
        self.ax.set_ylim(self.IV.plot_range_I())
        self.colors = ('b', 'g', 'r', 'c')


class IVViewerItIntegrated(ToggleViewerIntegrated):
    def __init__(self, IV, PP='PP'):
        self.IV, self.PP = IV, PP

        ToggleViewerIntegrated.__init__(self, 'I(t) integrated viewer')

    def plotfun(self, event):
        t_event = event.xdata
        V, I, Ifit, t = self.IV.get_Sfit_at_event(t_event, PP=self.PP)

        self.clear()
        lines = self.ax.plot(t, I, 'k-')
        lines_fit = self.ax.plot(t, Ifit, 'r-', linewidth=1.5)
        return lines + lines_fit


# method factories for IV and IVContainer
def _plot_factory(plotfun, PP):
    def plot(self, **kw):
        return getattr(self, plotfun)(PP=PP, **kw)
    return plot

def _plot_range_factory(rangefun):
    def plot_range(self, *args):
        r = np.array([getattr(x, rangefun)(*args) for x in self])
        return r[:,0].min(), r[:,1].max()
    return plot_range

def _forall_factory(method):
    def forall(self, *args, **kw):
        for x in self:
            getattr(x, method)(*args, **kw)
        return self
    return forall


class IV:
    def __init__(self, S, R):
        self.S, self.R = S, R

    @staticmethod
    def _slices(N, n, incr):
        sl = slice(0, N - n, incr)
        sr = slice(n, N, incr)
        return sl, sr

    def _get_valid_indices(self, i0, i1):
        w = self.R.plunges()
        isin = np.zeros_like(i0, bool)
        for a, b in w.T:
            isin[(a <= i0) & (i1 <= b)] = True
        return np.flatnonzero(isin)

    def _prepare(self, n, incr, nvars):
        iE = self.S.V.iE
        sl, sr = self._slices(len(iE), n, incr)
        i0, i1 = iE[sl], iE[sr]
        N = len(i0)

        ind = self._get_valid_indices(i0, i1)
        #ind = xrange(N)

        out = np.empty((N, nvars))
        out.fill(np.nan)

        shift = -(n // (2*incr))
        return sl, sr, i0, i1, ind, out, shift

    def _fit_const(self, FitterIVClass, n=1, incr=1, mask=None, **kw):
        sl, sr, i0, i1, ind, out, shift = self._prepare(n, incr, 3)
        
        S = self.S
        V, I, t = S.V.x, S.x, S.t
        self.mask = np.zeros(V.size, bool)

        if mask is None:
            def get_mask(s): return None
        else:
            def get_mask(s): return mask[s]

        for j in ind:
            s = slice(i0[j], i1[j] + 1)
            fitter_IV = FitterIVClass(V[s], I[s], mask=get_mask(s), **kw)
            try:
                out[j] = fitter_IV.fit()
                self.mask[s][fitter_IV.get_ind()] = True
            except FitterError:
                pass
        return PiecewisePolynomialEndpoints(out[None], t, i0=i0, i1=i1, shift=shift)

    def _fit_linear(self, FitterIVClass, n=5, incr=1, use_mask=True, **kw):
        sl, sr, i0, i1, ind, out, shift = self._prepare(n, incr, 6)

        S = self.S
        V, I, t = S.V.x, S.x, S.t
        t0, t1 = t[i0], t[i1]
        dt = t1 - t0

        try:
            # try fast version first
            c = self.PP.c[0]
            p_knots = np.concatenate((c[:1], 0.5*(c[:-1] + c[1:]), c[-1:]), axis=0)
            p = np.concatenate((p_knots[sl], p_knots[sr]), axis=1)
        except:
            # fallback to full evaluation
            p = np.concatenate((self.PP(t0), self.PP(t1)), axis=1)

        #Sfit = self.get_Sfit()
        #mask = ~Sfit.x.mask
        if use_mask:
            mask = self.mask
        else:
            mask = np.ones_like(self.mask)

        for j in ind:
            s = slice(i0[j], i1[j] + 1)
            m = mask[s]
            if not np.any(m):
                continue

            Vj, Ij, tj = V[s][m], I[s][m], t[s][m]

            aj = (tj - t0[j]) / dt[j]

            fitter_IV = FitterIVClass(Vj, Ij, aj, p[j], **kw)
            try:
                out[j] = fitter_IV.fit()
            except FitterError:
                pass

        c = np.swapaxes(out.reshape((-1, 2, 3)), 0, 1)[::-1].copy()
        c[0] -= c[1]
        c[0] /= dt[:, None]

        return PiecewisePolynomialEndpoints(c, t, i0=i0, i1=i1, shift=shift)

    def fit(self, **kw):
        self.PP = self._fit_const(FitterIV, **kw)
        return self

    def fitmag(self, **kw):
        self.PPmag = self._fit_const(FitterIVMag, **kw)
        return self

    def fit6(self, **kw):
        self.PP6 = self._fit_linear(FitterIV6, **kw)
        return self

    def fit6i(self, **kw):
        self.PP6i = self._fit_linear(FitterIV6i, **kw)
        return self

    def fit5(self, **kw):
        self.PP5 = self._fit_linear(FitterIV5, **kw)
        return self

    def fit4(self, **kw):
        self.PP4 = self._fit_linear(FitterIV4, **kw)
        return self

    @memoized_property
    def PP(self):
        print "Calculating PP..."
        return self._fit_const(FitterIV)

    @memoized_property
    def PPmag(self):
        print "Calculating PPmag..."
        return self._fit_const(FitterIVMag)
    
    @memoized_property
    def PP6(self):
        print "Calculating PP6..."
        return self._fit_linear(FitterIV6)

    @memoized_property
    def PP6i(self):
        print "Calculating PP6i..."
        return self._fit_linear(FitterIV6i)

    @memoized_property
    def PP5(self):
        print "Calculating PP5..."
        return self._fit_linear(FitterIV5)

    @memoized_property
    def PP4(self):
        print "Calculating PP4..."
        return self._fit_linear(FitterIV4)

    def get_PP(self, PP='PP'):
        return getattr(self, PP)

    def get_Sfit(self, PP='PP'):
        t, V = self.S.t, self.S.V
        p = getattr(self, PP)(t).T
        Ifit = FitterIV.fitfun(p, V.x)
        #mask = V.x > p[1] + p[2]
        Ifit_masked = ma.masked_array(Ifit, ~self.mask)
        return CurrentSignal(Ifit_masked, t, V=V)

    def get_Sfit_at_event(self, t_event, PP='PP'):
        PP = getattr(self, PP)
        s, p = PP.eval_at_event(t_event)

        S = self.S
        V, I, t = S.V.x[s], S.x[s], S.t[s]
        
        Ifit = FitterIV.fitfun(p.T, V)
        Ifit = ma.masked_array(Ifit, ~self.mask[s])
        return V, I, Ifit, t

    def plot_range_dt(self, PP='PP'):
        PP = getattr(self, PP)
        t = PP.x[[PP.i0[0], PP.i1[0]]]
        return (0, t[1] - t[0])

    def plot_range_V(self, r=0):
        return self.S.V.plot_range(r)

    def plot_range_I(self, r=0):
        return self.S.plot_range(r)

    def plot_raw(self, ax=None, PP='PP'):
        if ax is None:
            fig = get_tfig(xlab=self.S.xlab, ylab=self.S.ylab)
            ax = fig.axes[0]

            self.viewers = (IVViewer(self, PP=PP),
                            IVViewerIt(self, PP=PP),
                            IVViewerItIntegrated(self, PP=PP))
            
            menu_entries_ax = []
            for v in self.viewers:
                menu_entries_ax += v.menu_entries_ax

            fig.context_menu_picker = ContextMenuPicker(
                    fig, menu_entries_ax=menu_entries_ax)

        self.S.plot(ax=ax)
        self.get_Sfit(PP=PP).plot(ax=ax)
        return ax

    def plot(self, fig=None, PP='PP', **kw):
        if fig is None:
            xlab = 't (s)'
            ylab = ('Isat (A)', 'Vf (V)', 'Te (eV)')
            fig = get_tfig(shape=(3, 1), figsize=(10, 10),
                           xlab=xlab, ylab=ylab)

            self.viewers = (IVViewer(self, PP=PP),
                            IVViewerIt(self, PP=PP))
            
            menu_entries_ax = []
            for v in self.viewers:
                menu_entries_ax += v.menu_entries_ax

            fig.context_menu_picker = ContextMenuPicker(
                    fig, menu_entries_ax=menu_entries_ax)

        for ax, p in zip(fig.axes, getattr(self, PP)):
             p.plot(ax=ax, **kw)
        return fig

    plot_raw6 = _plot_factory('plot_raw', 'PP6')
    plot_raw6i = _plot_factory('plot_raw', 'PP6i')
    plot_raw5 = _plot_factory('plot_raw', 'PP5')
    plot_raw4 = _plot_factory('plot_raw', 'PP4')
    plot6 = _plot_factory('plot', 'PP6')
    plot6i = _plot_factory('plot', 'PP6i')
    plot5 = _plot_factory('plot', 'PP5')
    plot4 = _plot_factory('plot', 'PP4')


class IVContainer(Container):
    def __init__(self, S=None, R=None, x=None):
        if x is not None:
            self.x = x
        else:
            Container.__init__(self)
            for k, s in S.iteritems():
                self.x[k] = IV(s, R)

    def __getitem__(self, indx):
        try:
            return self.x[indx]
        except (KeyError, TypeError):
            try:
                k = np.array(self.x.keys())[indx]
                if isinstance(k, np.ndarray):
                    return self.__class__(x=DictView(self.x, k))
                else:
                    return self.x[k]
            except (ValueError, IndexError):
                raise KeyError(indx)

    def get_Sfit_at_event(self, t_event, PP='PP'):
        return [x.get_Sfit_at_event(t_event, PP) for x in self]

    plot_range_dt = _plot_range_factory('plot_range_dt')
    plot_range_V  = _plot_range_factory('plot_range_V')
    plot_range_I  = _plot_range_factory('plot_range_I')

    fit   = _forall_factory('fit')
    fit6  = _forall_factory('fit6')
    fit6i = _forall_factory('fit6i')
    fit5  = _forall_factory('fit5')
    fit4  = _forall_factory('fit4')

    def plot_raw(self, fig=None, PP='PP'):
        if fig is None:
            xlab = self.x.values()[0].S.xlab
            ylab = [x.S.ylab for x in self]
            fig = get_tfig(shape=(len(self.x), 1), figsize=(10, 10), 
                           xlab=xlab, ylab=ylab)

            self.viewers = (IVViewer(self, PP=PP),
                            IVViewerIt(self, PP=PP))

            menu_entries_ax = []
            for v in self.viewers:
                menu_entries_ax += v.menu_entries_ax

            fig.context_menu_picker = ContextMenuPicker(
                    fig, menu_entries_ax=menu_entries_ax)

        for ax, x in zip(fig.axes, self.x.values()):
            x.plot_raw(ax=ax, PP=PP)
        return fig

    def plot(self, fig=None, PP='PP', **kw):
        if fig is None:
            xlab = "t (s)"
            ylab = ('Isat (A)', 'Vf (V)', 'Te (eV)')
            fig = get_tfig(shape=(3, 1), figsize=(10, 10), 
                           xlab=xlab, ylab=ylab)

            self.viewers = (IVViewer(self, PP=PP),
                            IVViewerIt(self, PP=PP))

            menu_entries_ax = []
            for v in self.viewers:
                menu_entries_ax += v.menu_entries_ax

            fig.context_menu_picker = ContextMenuPicker(
                    fig, menu_entries_ax=menu_entries_ax)

        for x in self:
            for ax, p in zip(fig.axes, getattr(x, PP)):
                p.plot(ax=ax, **kw)
        return fig

    plot_raw6 = _plot_factory('plot_raw', 'PP6')
    plot_raw6i = _plot_factory('plot_raw', 'PP6i')
    plot_raw5 = _plot_factory('plot_raw', 'PP5')
    plot_raw4 = _plot_factory('plot_raw', 'PP4')
    plot6 = _plot_factory('plot', 'PP6')
    plot6i = _plot_factory('plot', 'PP6i')
    plot5 = _plot_factory('plot', 'PP5')
    plot4 = _plot_factory('plot', 'PP4')


