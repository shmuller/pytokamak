import numpy as np
import numpy.ma as ma

from itertools import cycle

from sig import median, memoized_property, get_fig, get_tfig
from sig import CurrentSignal, PiecewisePolynomialEndpoints

from fitter import Fitter, FitterError

try:
    import LP.fitfun as ff
except ImportError:
    import fitfun as ff

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
        self.red_fact = kw.pop('red_fact', 0.95)
        
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

        if self.red_fact < 1:
            save = self.X, self.Y
            Y0 = 0.
            while True:
                self.M *= self.red_fact
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


class FitterIV2(Fitter):
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


class FitterIV3(FitterIV2):
    def __init__(self, *args, **kw):
        FitterIV2.__init__(self, *args, **kw)

        self.perm = np.array((0,1,2,3))
        self.iperm = np.array((0,1,2,3,1,2))
        self.fact = self.fact[self.perm]

    def set_unnorm(self):
        self.p = (self.P * self.fact)[self.iperm]
        return self.p, self.p0

    def set_guess(self):
        self.P0 = self.p0[self.perm] / self.fact

    custom_engine = ff.IV4_fit
    fitfun_fast = ff.IV4
    fitfun_diff = ff.IV4_diff
    fitfun_rms = ff.IV4_rms


class IVSeriesSimpleViewerIV(ToggleViewer):
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


class IVSeriesSimpleViewerIt(ToggleViewer):
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


class IVSeriesSimpleViewerItIntegrated(ToggleViewerIntegrated):
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


class IVSeriesSimple:
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

    def _fit_const(self, n=1, incr=1, **kw):
        sl, sr, i0, i1, ind, out, shift = self._prepare(n, incr, 3)
        
        S = self.S
        V, I, t = S.V.x, S.x, S.t
        self.mask = np.zeros(V.size, bool)

        for j in ind:
            s = slice(i0[j], i1[j])
            fitter_IV = FitterIV(V[s], I[s], **kw)
            try:
                out[j] = fitter_IV.fit()
                self.mask[s][fitter_IV.get_ind()] = True
            except FitterError:
                pass
        return PiecewisePolynomialEndpoints(out[None], t, i0=i0, i1=i1, shift=shift)

    def _fit_linear(self, FitterIVClass, n=5, incr=1, **kw):
        sl, sr, i0, i1, ind, out, shift = self._prepare(n, incr, 6)

        c = self.PP.c[0]
        p_knots = np.concatenate((c[:1], 0.5*(c[:-1] + c[1:]), c[-1:]), axis=0)

        p = np.concatenate((p_knots[sl], p_knots[sr]), axis=1)
       
        S = self.S
        V, I, t = S.V.x, S.x, S.t
        t0, t1 = t[i0], t[i1]
        dt = t1 - t0

        #Sfit = self.get_Sfit()
        #mask = ~Sfit.x.mask
        mask = self.mask
        #mask = np.ones_like(self.mask)

        for j in ind:
            s = slice(i0[j], i1[j])
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
        self.PP = self._fit_const(**kw)

    def fit2(self, **kw):
        self.PP2 = self._fit_linear(FitterIV2, **kw)

    def fit3(self, **kw):
        self.PP3 = self._fit_linear(FitterIV3, **kw)

    @memoized_property
    def PP(self):
        print "Calculating PP..."
        return self._fit_const()

    @memoized_property
    def PP2(self):
        print "Calculating PP2..."
        return self._fit_linear(FitterIV2)

    @memoized_property
    def PP3(self):
        print "Calculating PP3..."
        return self._fit_linear(FitterIV3)

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
        dtM = np.diff(PP.x[np.array((PP.i0[0], PP.i1[0]))])
        return (0, dtM)

    def plot_range_V(self, r=0):
        return self.S.V.plot_range(r)

    def plot_range_I(self, r=0):
        return self.S.plot_range(r)

    def plot_raw(self, ax=None, PP='PP'):
        if ax is None:
            fig = get_tfig()
            ax = fig.axes[0]

            self.viewers = (IVSeriesSimpleViewerIV(self, PP=PP),
                            IVSeriesSimpleViewerIt(self, PP=PP),
                            IVSeriesSimpleViewerItIntegrated(self, PP=PP))
            
            menu_entries_ax = []
            for v in self.viewers:
                menu_entries_ax += v.menu_entries_ax

            fig.context_menu_picker = ContextMenuPicker(
                    fig, menu_entries_ax=menu_entries_ax)

        self.S.plot(ax=ax)
        self.get_Sfit(PP=PP).plot(ax=ax)
        return ax

    def plot_raw2(self, ax=None):
        return self.plot_raw(ax=ax, PP='PP2')

    def plot(self, fig=None, PP='PP'):
        if fig is None:
            xlab = 't (s)'
            ylab = ('Isat (A)', 'Vf (V)', 'Te (eV)')
            fig = get_tfig(shape=(3, 1), figsize=(10, 10),
                           xlab=xlab, ylab=ylab)

            self.viewers = (IVSeriesSimpleViewerIV(self, PP=PP),
                            IVSeriesSimpleViewerIt(self, PP=PP))
            
            menu_entries_ax = []
            for v in self.viewers:
                menu_entries_ax += v.menu_entries_ax

            fig.context_menu_picker = ContextMenuPicker(
                    fig, menu_entries_ax=menu_entries_ax)

        for ax, p in zip(fig.axes, getattr(self, PP)):
             p.plot(ax=ax)
        return fig

    def plot2(self, fig=None):
        return self.plot(fig=fig, PP='PP2')


class IVSeriesSimpleGroup:
    def __init__(self, S, R):
        self.x = [IVSeriesSimple(s, R) for s in S]

    def get_Sfit_at_event(self, t_event, PP='PP'):
        return [x.get_Sfit_at_event(t_event, PP) for x in self.x]

    def _plot_range_factory(range_fun):
        def plot_range(self, *args):
            r = np.array([getattr(x, range_fun)(*args) for x in self.x])
            return r[:,0].min(), r[:,1].max()
        return plot_range

    plot_range_dt = _plot_range_factory('plot_range_dt')
    plot_range_V  = _plot_range_factory('plot_range_V')
    plot_range_I  = _plot_range_factory('plot_range_I')

    def plot_raw(self, fig=None, PP='PP'):
        n = len(self.x)
        xlab = "t (%s)" % self.x[0].S.tunits
        ylab = ["%s (%s)" % (x.S.name, x.S.units) for x in self.x]
        fig = get_tfig(shape=(n, 1), figsize=(10, 10), xlab=xlab, ylab=ylab)

        self.viewers = (IVSeriesSimpleViewerIV(self, PP=PP),
                        IVSeriesSimpleViewerIt(self, PP=PP))

        menu_entries_ax = []
        for v in self.viewers:
            menu_entries_ax += v.menu_entries_ax

        fig.context_menu_picker = ContextMenuPicker(
                fig, menu_entries_ax=menu_entries_ax)

        for ax, x in zip(fig.axes, self.x):
            x.plot_raw(ax=ax, PP=PP)
        return fig

    def plot_raw2(self, fig=None):
        return self.plot_raw(fig=fig, PP='PP2')

    def plot(self, fig=None, PP='PP'):
        xlab = "t (s)"
        ylab = ('Isat (A)', 'Vf (V)', 'Te (eV)')
        fig = get_tfig(shape=(3, 1), figsize=(10, 10), xlab=xlab, ylab=ylab)

        self.viewers = (IVSeriesSimpleViewerIV(self, PP=PP),
                        IVSeriesSimpleViewerIt(self, PP=PP))

        menu_entries_ax = []
        for v in self.viewers:
            menu_entries_ax += v.menu_entries_ax

        fig.context_menu_picker = ContextMenuPicker(
                fig, menu_entries_ax=menu_entries_ax)

        for x in self.x:
            for ax, p in zip(fig.axes, getattr(x, PP)):
                p.plot(ax=ax)
        return fig

    def plot2(self, fig=None):
        return self.plot(fig=fig, PP='PP2')


