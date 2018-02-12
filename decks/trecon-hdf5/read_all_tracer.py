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

tracer = 'tracer/ion_tracer.9961473'
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

gama = np.sqrt(1.+px*px+py*py+pz*pz)
ee = gama-1.
emax = max(ee)
print "Emax", emax
dt = t[1]-t[0]
print "dt", dt
vx = px/gama
vy = py/gama
vz = pz/gama
exi = - (uy*bz - uz*by)
eyi = - (uz*bx - ux*bz)
ezi = - (ux*by - uy*bx)
exn = ex - exi
eyn = ey - eyi
ezn = ez - ezi

Egain = np.zeros(nt)

Egain[0] = ee[0]
for i in range (1, nt-1):
  Egain[i] = Egain[i-1] + dt*(ex[i]*vx[i] + ey[i]*vy[i] + ez[i]*vz[i])



fig = plt.figure(figsize=[7, 5])
xs, ys = 0.15, 0.15
w1, h1 = 0.8, 0.8
ax = fig.add_axes([xs, ys, w1, h1])
ax.plot(t, ee, linewidth=2, color='b')
ax.plot(t, Egain+1, linewidth=2, color='g')
#ax.plot(x, z, linewidth=2)
ax.tick_params(labelsize=16)

plt.show()
