import matplotlib
matplotlib.use('qt4agg')
import matplotlib.pyplot as plt

from LP import probe_xpr

#shn = [27688, 27689, 27690, 27691, 27692]
shn = [28444, 28445, 28446]

fig = None

for s in shn:
    XPR = probe_xpr.ProbeXPR(shn=s)

    for i in xrange(XPR['R'].Nplunges):
        fig = XPR.res.plot(xkey='R', fig=fig, plunge=i, inout='in')

plt.show()
    

