import math
import os
import random
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

def calc_jdote_single(fname):
    """
    """
    fdata = np.fromfile(fname, dtype=np.float32)
    sz, = fdata.shape
    nt  = sz / 16

    print nt

    jdote1 = np.zeros(nt)
    jdote2 = np.zeros(nt)
    gamma_t = np.zeros(nt)

    dt = fdata[16] - fdata[0]
    
    q = 1.0 if 'ion' in fname else -1.0

    for ct in range(nt):
        ptl = fdata[ct*16:(ct+1)*16]
        ux = ptl[4]
        uy = ptl[5]
        uz = ptl[6]
        ex = ptl[7]
        ey = ptl[8]
        ez = ptl[9]
        bx = ptl[10]
        by = ptl[11]
        bz = ptl[12]
        vx = ptl[13]
        vy = ptl[14]
        vz = ptl[15]

        einx = by*vz - bz*vy
        einy = bz*vx - bx*vz
        einz = bx*vy - by*vx
        gamma = math.sqrt(1 + ux*ux + uy*uy + uz*uz)
        gamma_t[ct] = gamma
        ux /= gamma
        uy /= gamma
        uz /= gamma
        jdote1[ct] = ex * ux + ey * uy + ez * uz
        jdote2[ct] = einx * ux + einy * uy + einz * uz

    ene = gamma_t - 1
    dene_dt = np.gradient(ene, dt)
    dgamma_dt = np.gradient(gamma_t, dt)
    jdote1_cum = np.cumsum(jdote1) * dt * q
    jdote2_cum = np.cumsum(jdote2) * dt * q

    return (ene, jdote1_cum, jdote2_cum)


def plot_jdote_single(fname):
    """
    """
    print(fname)
    ene, jdote1_cum, jdote2_cum = calc_jdote_single(fname)
    fig = plt.figure(figsize=[7, 5])
    xs, ys = 0.15, 0.15
    w1, h1 = 0.8, 0.8
    ax = fig.add_axes([xs, ys, w1, h1])
    ax.plot(jdote1_cum, linewidth=2, color='r')
    ax.plot(jdote2_cum, linewidth=2, color='b')
    ax.plot(ene-ene[0], linewidth=2, color='k')
    ax.tick_params(labelsize=16)

    plt.show()


def calc_total_jdote():
    """
    Loop over all particles and calculate the total energy gain
    """
    tracer_dir = 'tracer/'
    fname = random.choice(os.listdir(tracer_dir))
    fname = tracer_dir + fname
    ene, jdote1_cum, jdote2_cum = calc_jdote_single(fname)
    nt, = ene.shape
    ene = np.zeros(nt)
    qve = np.zeros(nt)
    qv_uxb = np.zeros(nt)
    for fname in os.listdir(tracer_dir):
        fname = tracer_dir + fname
        ene_s, data1, data2 = calc_jdote_single(fname)
        qve += data1
        qv_uxb += data2
        ene += ene_s

    fig = plt.figure(figsize=[7, 5])
    xs, ys = 0.15, 0.15
    w1, h1 = 0.8, 0.8
    ax = fig.add_axes([xs, ys, w1, h1])
    ax.plot(qve, linewidth=2, color='r')
    ax.plot(qv_uxb, linewidth=2, color='b')
    ax.plot(ene-ene[0], linewidth=2, color='k')
    ax.tick_params(labelsize=16)

    plt.show()
    

if __name__ == "__main__":
    tracer_dir = 'tracer/'
    fname = random.choice(os.listdir(tracer_dir))
    fname = tracer_dir + fname
    # plot_jdote_single(fname)
    calc_total_jdote()
