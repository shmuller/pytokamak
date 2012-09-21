import matplotlib.pyplot as plt

import probe_xpr

#shn = [27688, 27689, 27690, 27691, 27692]
shn = [28444, 28445, 28446]

fig = None
plt.ion()

for s in shn:
    XPR = probe_xpr.ProbeXPR(shn=s)
    XPR.load()
    XPR.trim(plunge=0)
    XPR.analyze()
    fig = XPR.plot_R(fig=fig)
    plt.show()
    plt.draw()

