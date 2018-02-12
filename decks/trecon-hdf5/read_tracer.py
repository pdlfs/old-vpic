import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]

font = {
    'family': 'serif',
    # 'color':'darkred',
    'color': 'black',
    'weight': 'normal',
    'size': 24,
}

tracer = 'tracer/ion_tracer.40108035'
fdata = np.fromfile(tracer, dtype=np.float32)
sz, = fdata.shape
nt  = sz / 16

ptl = fdata.reshape((nt, 16))
t = ptl[:, 0]
x = ptl[:, 1]
y = ptl[:, 2]
z = ptl[:, 3]
px = ptl[:, 4]
py = ptl[:, 5]
pz = ptl[:, 6]
ex = ptl[:, 7]
ey = ptl[:, 8]
ez = ptl[:, 9]
bx = ptl[:, 10]
by = ptl[:, 11]
bz = ptl[:, 12]
ux = ptl[:, 13]
uy = ptl[:, 14]
uz = ptl[:, 15]

ee = np.sqrt(1.+px*px+py*py+pz*pz)-1.
emax = max(ee)
print "Emax", emax

fig = plt.figure(figsize=[7, 5])
xs, ys = 0.15, 0.15
w1, h1 = 0.8, 0.8
ax = fig.add_axes([xs, ys, w1, h1])
ax.plot(t, ee, linewidth=2)
#ax.plot(x, z, linewidth=2)
ax.tick_params(labelsize=16)

plt.show()
