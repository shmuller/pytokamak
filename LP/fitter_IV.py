import numpy as np
import numpy.ma as ma

from itertools import cycle

from sig import memoized_property, median
from sig import DictView, GeneratorDict, Container, \
        CurrentSignal, PiecewisePolynomialEndpoints

from fitter import Fitter, FitterError

import fitfun as ff
#import fitfun_ctypes as ff
#import fitfun_cython as ff
#import fitfun_cffi as ff

from LP import mag_fit

from sm_pyplot.tight_figure import get_tfig
from sm_pyplot.contextmenupicker import ContextMenuPicker
from sm_pyplot.observer_viewer import ToggleViewer, ToggleViewerIntegrated


class FitterIVBase(Fitter):
    linear = False
    nvars = 3
    def __init__(self, V, I, mask=None, cut_at_min=False, r=1, **kw):
        self.cut_at_min, self.r = cut_at_min, r

        self.mask_ind = np.arange(V.size)
        if mask is not None:
            self.mask_ind = self.mask_ind[mask]

        self.sort_ind = V[self.mask_ind].argsort()
        self.ind = self.mask_ind[self.sort_ind]

        self.V, self.I = V[self.ind], I[self.ind]
        self.M = self.V.size
        
        self.im = self.I.argmin()
        self.Vm = Vm = self.V[0]
        self.VM = VM = self.V[-1]
        self.Im = Im = self.I[self.im]
        self.IM = IM = median(self.I[:self.M/2])
        self.dV = dV = VM - Vm
        self.dI = dI = IM - Im

        self.abcd = (2./dI, -(2.*Im/dI + 1.), 1./dV, -Vm/dV)

        Fitter.__init__(self, self.V, self.I, **kw)

    @classmethod
    def zero_crossings(cls, x):
        return np.flatnonzero(np.diff(np.sign(x)))

    def get_ind(self):
        return self.ind[:self.M]

    def normalize(self, V, I):
        if self.cut_at_min:
            self.M = self.im+1
        a, b, c, d = self.abcd
        X = c*V[:self.M].astype('d') + d
        Y = a*I[:self.M].astype('d') + b
        return X, Y

    def unnormalize(self, P):
        return self.LP_unnormalize(P)

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

    def is_old_better(self, P_old):
        return np.any(self.P > P_old)
        #return np.any(self.P > P_old) or (self.eval_norm(self.X[-1]) > 0.)

    def fit(self):
        Fitter.fit(self)

        if self.r < 1:
            Y0 = 0.
            save = self.X, self.Y
            while True:
                P_old, M_old = self.P, self.M
                self.M *= self.r
                self.X, self.Y = self.X[:self.M], self.Y[:self.M]
                Fitter.fit(self)

                if self.is_old_better(P_old) or (self.Y[-1] > Y0):
                    self.P, self.M = P_old, M_old
                    break
            self.X, self.Y = save

        self.check()

    def check(self):
        n, Vf, Te = self.p[:3]
        if n < 0.:
            raise FitterError("Negative n")
        if Te < 0. or Te > 0.5*self.dV:
            raise FitterError("Unrealistic Te")


class FitterIV(FitterIVBase):
    def __init__(self, V, I, **kw):
        kw.setdefault('cut_at_min', True)
        kw.setdefault('r', 0.95)
        FitterIVBase.__init__(self, V, I, **kw)

    def set_guess(self):
        if not self.OK:
            raise FitterError("Cannot compute guess for data that failed OK check")

        V, I = self.X, self.Y
        I0 = 1.
        i0 = self.zero_crossings(I)
        Vf = np.mean((V[i0] + V[i0+1])/2)
        
        #I1 = 0.5*I0
        #i1 = self.zero_crossings(I - I1)
        #V1 = np.mean((V[i1] + V[i1+1])/2)
        V1, I1 = V[self.im], I[self.im]

        Te = (V1-Vf)/np.log(1-I1/I0)

        self.P0 = np.array((I0, Vf, Te))

    def LP_unnormalize(self, P):
        a, b, c, d = self.abcd

        p = P.copy()
        p[0] = (P[0]-b)/a
        p[1] = (P[1]+P[2]*np.log(1.-b/P[0]) - d)/c
        p[2] = P[2]/c
        return p

    def LP_normalize(self, p):
        a, b, c, d = self.abcd

        P = p.copy()
        P[0] = a*p[0]+b
        P[1] = c*(p[1]-p[2]*np.log(1.-b/P[0])) + d
        P[2] = c*p[2]
        return P

    def is_old_better(self, P_old):
        dP = self.P - P_old
        return np.any(dP[[0, 2]] > 0)

    @classmethod
    def fitfun(cls, P, X):
        iP2 = 1./P[2]
        return P[0]*(1.-np.exp((X-P[1])*iP2))

    custom_engine = staticmethod(ff.IV3_fit)
    fitfun_fast = staticmethod(ff.IV3)
    fitfun_diff = staticmethod(ff.IV3_diff)
    fitfun_rms = staticmethod(ff.IV3_rms)


class FitterIVDbl(FitterIVBase):
    nvars = 5
    def __init__(self, V, I, **kw):
        FitterIVBase.__init__(self, V, I, **kw)

        self.do_var = np.array((1, 1, 1, 0, 2), 'i')

    def LP_unnormalize(self, P):
        a, b, c, d = self.abcd

        p = P.copy()
        p[0] = (P[0] - b)/a
        p[4] = (P[0] - b)*P[4]/(P[0] + b*P[4])
        p[1] = (P[1] + P[2]*np.log(1. - b*(1.+p[4])/P[0]) - d)/c
        p[2] = P[2]/c
        return p

    def set_guess(self):
        if not self.OK:
            raise FitterError("Cannot compute guess for data that failed OK check")

        V, I = self.X, self.Y
        I0, B = 1., 1.
        
        i0 = self.zero_crossings(I)
        Vf = np.mean((V[i0] + V[i0+1])/2)
        
        I1 = 0.5*I0
        i1 = self.zero_crossings(I - I1)
        V1 = np.mean((V[i1] + V[i1+1])/2)
        
        Te = (V1-Vf)/np.log((1-I1/I0)/(1+B))

        self.P0 = np.array((I0, Vf, Te, 0., B))
    
    @classmethod
    def pow_075(cls, x):
        cnd = x < 0.
        y = np.zeros_like(x)
        tmp = np.sqrt(1. - 2.*x[cnd])
        y[cnd] = (tmp + 2.)*np.sqrt(tmp - 1.)
        return y

    @classmethod
    def fitfun(cls, P, X):
        iP2 = 1./P[2]
        arg = (P[1] - X)*iP2
        exp_arg = np.exp(arg)
        return P[0]*(exp_arg - 1. - P[3]*cls.pow_075(arg)) / (exp_arg + P[4])

    custom_engine = staticmethod(ff.IVdbl_fit)
    fitfun_fast = staticmethod(ff.IVdbl)
    fitfun_diff = staticmethod(ff.IVdbl_diff)
    fitfun_rms = staticmethod(ff.IVdbl_rms)


class FitterIVMag(Fitter):
    linear = False
    nvars = 20
    def __init__(self, V, I, mask=None, **kw):
        Fitter.__init__(self, V, I, **kw)
        
        self.engine = self.magfit

        self.do_var = np.array([1, 1, 0, 1, 1, 1, 
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 'i')

        self.yfit = np.zeros_like(self.y)

    def get_ind(self):
        return np.arange(self.x.size)

    def set_guess(self):
        self.P0 = np.array([5e18, 18., 0., 12., 18., 1e-8, 
            np.pi/2, 0., 0., 0., 0.0009, -1., 2., 1., 2., 0., 1., 0.0016/2, 0., 0.])

    @classmethod
    def fitfun(cls, P, X):
        do_var = np.zeros(P.shape[0], 'i')
        Y = np.zeros_like(X)

        mag_fit.magfit(X, Y.copy(), Y, P[:,0].copy(), do_var)
        return -Y

    def magfit(self, p0, x, y):
        p = p0.copy()
        mag_fit.magfit(x, -y, self.yfit, p, self.do_var)
        return p


class FitterIVLinear(Fitter):
    linear = True
    def __init__(self, V, I, a, p0, **kw):
        Fitter.__init__(self, V, I, args=(a,), **kw)

        self.dV = V.ptp()
        self.dI = I.ptp()
        self.p0 = p0

    def normalize(self, V, I):
        return V/self.dV, I/self.dI

    def unnormalize(self, P):
        return P * self.fact

    def set_OK(self):
        self.OK = np.isfinite(self.p0).all()

    def set_guess(self):
        self.P0 = self.p0 / self.fact


class FitterIV6(FitterIVLinear):
    nvars, base_ID = 6, 'IV'
    def __init__(self, *args, **kw):
        FitterIVLinear.__init__(self, *args, **kw)

        dV, dI = self.dV, self.dI
        self.fact = np.array((dI, dV, dV, dI, dV, dV))

    @classmethod
    def fitfun(cls, p, V, a):
        n = cls.nvars/2
        pp = p[:,None]
        pj = pp[:n] + a[None]*(pp[n:] - pp[:n])
        return FitterIV.fitfun(pj, V)

    custom_engine = staticmethod(ff.IV6_fit)
    fitfun_fast = staticmethod(ff.IV6)
    fitfun_diff = staticmethod(ff.IV6_diff)
    fitfun_rms = staticmethod(ff.IV6_rms)


class FitterIV6Perm(FitterIV6):
    def __init__(self, *args, **kw):
        # derived class needs to assign self.perm and self.iperm here

        FitterIV6.__init__(self, *args, **kw)
        self.fact = self.fact[self.perm]

    def unnormalize(self, P):
        return (P * self.fact)[self.iperm]

    def set_guess(self):
        self.P0 = self.p0[self.perm] / self.fact


class FitterIV5(FitterIV6Perm):
    def __init__(self, *args, **kw):
        self.perm = np.array((0,1,2,3,4))
        self.iperm = np.array((0,1,2,3,4,2))
        FitterIV6Perm.__init__(self, *args, **kw)

    custom_engine = staticmethod(ff.IV5_fit)
    fitfun_fast = staticmethod(ff.IV5)
    fitfun_diff = staticmethod(ff.IV5_diff)
    fitfun_rms = staticmethod(ff.IV5_rms)


class FitterIV4(FitterIV6Perm):
    def __init__(self, *args, **kw):
        self.perm = np.array((0,1,2,3))
        self.iperm = np.array((0,1,2,3,1,2))
        FitterIV6Perm.__init__(self, *args, **kw)

    custom_engine = staticmethod(ff.IV4_fit)
    fitfun_fast = staticmethod(ff.IV4)
    fitfun_diff = staticmethod(ff.IV4_diff)
    fitfun_rms = staticmethod(ff.IV4_rms)


class FitterIVDbl2(FitterIVLinear):
    nvars, base_ID = 10, 'IVdbl'
    def __init__(self, *args, **kw):
        FitterIVLinear.__init__(self, *args, **kw)

        dV, dI = self.dV, self.dI
        self.fact = np.array((dI, dV, dV, 1., 1., dI, dV, dV, 1., 1.))

        self.do_var = np.array((1, 1, 1, 0, 2) * 2, 'i')
    
    @classmethod
    def fitfun(cls, p, V, a):
        n = cls.nvars/2
        pp = p[:,None]
        pj = pp[:n] + a[None]*(pp[n:] - pp[:n])
        return FitterIVDbl.fitfun(pj, V)

    custom_engine = staticmethod(ff.IVdbl2_fit)
    fitfun_fast = staticmethod(ff.IVdbl2)
    fitfun_diff = staticmethod(ff.IVdbl2_diff)
    fitfun_rms = staticmethod(ff.IVdbl2_rms)


FitterIVClasses = dict(
        IV     = FitterIV,
        IVdbl  = FitterIVDbl,
        IVmag  = FitterIVMag,
        IV6    = FitterIV6,
        IV5    = FitterIV5,
        IV4    = FitterIV4,
        IVdbl2 = FitterIVDbl2)


class IVViewer(ToggleViewer):
    def __init__(self, IV, ID='IV'):
        self.IV, self.ID = IV, ID

        ToggleViewer.__init__(self, 'IV viewer')

    def plotfun(self, event):
        t_event = event.xdata
        res = self.IV.get_Sfit_at_event(t_event, ID=self.ID)
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
        self.colors = ('b', 'g', 'r', 'c', 'm')


class IVViewerIt(ToggleViewer):
    def __init__(self, IV, ID='IV'):
        self.IV, self.ID = IV, ID

        ToggleViewer.__init__(self, 'I(t) viewer')

    def plotfun(self, event):
        t_event = event.xdata
        res = self.IV.get_Sfit_at_event(t_event, ID=self.ID)
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

        self.ax.set_xlim(self.IV.plot_range_dt(self.ID))
        self.ax.set_ylim(self.IV.plot_range_I())
        self.colors = ('b', 'g', 'r', 'c', 'm')


class IVViewerItIntegrated(ToggleViewerIntegrated):
    def __init__(self, IV, ID='IV'):
        self.IV, self.ID = IV, ID

        ToggleViewerIntegrated.__init__(self, 'I(t) integrated viewer')

    def plotfun(self, event):
        t_event = event.xdata
        V, I, Ifit, t = self.IV.get_Sfit_at_event(t_event, ID=self.ID)

        self.clear()
        lines = self.ax.plot(t, I, 'k-')
        lines_fit = self.ax.plot(t, Ifit, 'r-', linewidth=1.5)
        return lines + lines_fit


class IV:
    def __init__(self, S, R):
        self.S, self.R = S, R
        self.PP = GeneratorDict(generator=self._fit)

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
        nvars = FitterIVClass.nvars
        sl, sr, i0, i1, ind, out, shift = self._prepare(n, incr, nvars)
        
        S = self.S
        V, I, t = S.V.x, S.x, S.t
        out_mask = np.zeros(V.size, bool)

        if mask is None:
            def get_mask(s): return None
        else:
            def get_mask(s): return mask[s]

        for j in ind:
            s = slice(i0[j], i1[j] + 1)
            fitter_IV = FitterIVClass(V[s], I[s], mask=get_mask(s), **kw)
            try:
                out[j] = fitter_IV.p
                out_mask[s][fitter_IV.get_ind()] = True
            except FitterError:
                pass

        PP = PiecewisePolynomialEndpoints(out[None], t, i0=i0, i1=i1, shift=shift)
        PP.mask = out_mask
        return PP

    def _fit_linear(self, FitterIVClass, n=5, incr=1, use_mask=True, **kw):
        nvars, base_ID = FitterIVClass.nvars, FitterIVClass.base_ID
        sl, sr, i0, i1, ind, out, shift = self._prepare(n, incr, nvars)

        S = self.S
        V, I, t = S.V.x, S.x, S.t
        t0, t1 = t[i0], t[i1]
        dt = t1 - t0

        base_PP = self.PP[base_ID]
        try:
            # try fast version first
            c = base_PP.c[0]
            p_knots = np.concatenate((c[:1], 0.5*(c[:-1] + c[1:]), c[-1:]), axis=0)
            p = np.concatenate((p_knots[sl], p_knots[sr]), axis=1)
        except:
            # fallback to full evaluation
            p = np.concatenate((base_PP(t0), base_PP(t1)), axis=1)

        if use_mask:
            mask = base_PP.mask
        else:
            mask = np.ones(V.size, bool)

        for j in ind:
            s = slice(i0[j], i1[j] + 1)
            m = mask[s]
            if not np.any(m):
                continue

            Vj, Ij, tj = V[s][m], I[s][m], t[s][m]

            aj = (tj - t0[j]) / dt[j]

            fitter_IV = FitterIVClass(Vj, Ij, aj, p[j], **kw)
            try:
                out[j] = fitter_IV.p
            except FitterError:
                pass

        c = np.swapaxes(out.reshape((-1, 2, nvars/2)), 0, 1)[::-1].copy()
        c[0] -= c[1]
        c[0] /= dt[:, None]

        PP = PiecewisePolynomialEndpoints(c, t, i0=i0, i1=i1, shift=shift)
        PP.mask = mask
        return PP

    def _merge(self, ID1='IVdbl2', ID2='IV6', Te_c=10.):
        self.fit('IVdbl', engine='odr')
        self.fit(ID1, engine='odr')
        self.fit('IV', r=1)
        self.fit(ID2, engine='odr')

        PP = self.PP[ID1].copy()
        cnd = PP[2] > Te_c
        PP.c[:,cnd,:3] = self.PP[ID2].c[:,cnd]
        PP.c[:,cnd,3:] = 0.
        return PP

    def _fit(self, ID='IV', **kw):
        if ID == 'IVhyb':
            return self._merge(**kw)
        else:
            FitterIVClass = FitterIVClasses[ID]
            if FitterIVClass.linear:
                fit = self._fit_linear
            else:
                fit = self._fit_const
            return fit(FitterIVClass, **kw)

    def fit(self, ID='IV', **kw):
        self.PP[ID] = self._fit(ID=ID, **kw)
        return self

    @classmethod
    def fitfun(cls, P, X):
        nvars = P.shape[0]
        if nvars == 3:
            return FitterIV.fitfun(P, X)
        elif nvars == 5:
            return FitterIVDbl.fitfun(P, X)
        elif nvars == 20:
            return FitterIVMag.fitfun(P, X)
        else:
            raise Exception("No matching fit function")

    def get_Sfit(self, ID='IV'):
        t, V = self.S.t, self.S.V
        PP = self.PP[ID]
        p = PP(t).T
        Ifit = self.fitfun(p, V.x)
        #mask = V.x > p[1] + p[2]
        Ifit_masked = ma.masked_array(Ifit, ~PP.mask)
        return CurrentSignal(Ifit_masked, t, V=V)

    def get_Sfit_at_event(self, t_event, ID='IV'):
        PP = self.PP[ID]
        s, p = PP.eval_at_event(t_event)

        S = self.S
        V, I, t = S.V.x[s], S.x[s], S.t[s]
        
        Ifit = self.fitfun(p.T, V)
        Ifit = ma.masked_array(Ifit, ~PP.mask[s])
        return V, I, Ifit, t

    def plot_range_dt(self, ID='IV'):
        PP = self.PP[ID]
        t = PP.x[[PP.i0[0], PP.i1[0]]]
        return (0, t[1] - t[0])

    def plot_range_V(self, r=0):
        return self.S.V.plot_range(r)

    def plot_range_I(self, r=0):
        return self.S.plot_range(r)

    def plot_raw(self, ID='IV', ax=None):
        if ax is None:
            fig = get_tfig(xlab=self.S.xlab, ylab=self.S.ylab)
            ax = fig.axes[0]

            self.viewers = (IVViewer(self, ID=ID),
                            IVViewerIt(self, ID=ID),
                            IVViewerItIntegrated(self, ID=ID))
            
            menu_entries_ax = []
            for v in self.viewers:
                menu_entries_ax += v.menu_entries_ax

            fig.context_menu_picker = ContextMenuPicker(
                    fig, menu_entries_ax=menu_entries_ax)

        self.S.plot(ax=ax)
        self.get_Sfit(ID=ID).plot(ax=ax)
        return ax

    def plot(self, ID='IV', fig=None, **kw):
        if fig is None:
            xlab = 't (s)'
            ylab = ('Isat (A)', 'Vf (V)', 'Te (eV)')
            fig = get_tfig(shape=(3, 1), figsize=(10, 10),
                           xlab=xlab, ylab=ylab)

            self.viewers = (IVViewer(self, ID=ID),
                            IVViewerIt(self, ID=ID))
            
            menu_entries_ax = []
            for v in self.viewers:
                menu_entries_ax += v.menu_entries_ax

            fig.context_menu_picker = ContextMenuPicker(
                    fig, menu_entries_ax=menu_entries_ax)

        for ax, p in zip(fig.axes, self.PP[ID]):
             p.plot(ax=ax, **kw)
        return fig


# method factories for IVContainer
def _plot_range_factory(rangefun):
    def plot_range(self, *args):
        r = np.array([getattr(x, rangefun)(*args) for x in self])
        return r[:,0].min(), r[:,1].max()
    return plot_range


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

    def get_Sfit_at_event(self, t_event, ID='IV'):
        return [x.get_Sfit_at_event(t_event, ID) for x in self]

    plot_range_dt = _plot_range_factory('plot_range_dt')
    plot_range_V  = _plot_range_factory('plot_range_V')
    plot_range_I  = _plot_range_factory('plot_range_I')

    def fit(self, *args, **kw):
        for x in self:
            x.fit(*args, **kw)
        return self
    
    def plot_raw(self, ID='IV', fig=None):
        if fig is None:
            xlab = self.x.values()[0].S.xlab
            ylab = [x.S.ylab for x in self]
            fig = get_tfig(shape=(len(self.x), 1), figsize=(10, 10), 
                           xlab=xlab, ylab=ylab)

            self.viewers = (IVViewer(self, ID=ID),
                            IVViewerIt(self, ID=ID))

            menu_entries_ax = []
            for v in self.viewers:
                menu_entries_ax += v.menu_entries_ax

            fig.context_menu_picker = ContextMenuPicker(
                    fig, menu_entries_ax=menu_entries_ax)

        for ax, x in zip(fig.axes, self.x.values()):
            x.plot_raw(ax=ax, ID=ID)
        return fig

    def plot(self, ID='IV', fig=None, **kw):
        if fig is None:
            xlab = "t (s)"
            ylab = ('Isat (A)', 'Vf (V)', 'Te (eV)')
            fig = get_tfig(shape=(3, 1), figsize=(10, 10), 
                           xlab=xlab, ylab=ylab)

            self.viewers = (IVViewer(self, ID=ID),
                            IVViewerIt(self, ID=ID))

            menu_entries_ax = []
            for v in self.viewers:
                menu_entries_ax += v.menu_entries_ax

            fig.context_menu_picker = ContextMenuPicker(
                    fig, menu_entries_ax=menu_entries_ax)

        for x in self:
            for ax, p in zip(fig.axes, x.PP[ID]):
                p.plot(ax=ax, **kw)
        return fig

    
