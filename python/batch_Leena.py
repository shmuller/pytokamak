import matplotlib.pyplot as plt

import probe_xpr

#shn = [27688, 27689, 27690, 27691, 27692]
shn = [28444, 28445, 28446]

fig = None

for s in shn:
    XPR = probe_xpr.ProbeXPR(shn=s)
    res = XPR.results()
    res.save()

    for i in xrange(XPR.iplunges.size):
        kw = dict(plunge=i, inout='in')
        res.save(**kw)
        fig = res.plot(xkey='R', fig=fig, **kw)

plt.show()
    

