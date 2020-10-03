import os
import matplotlib as mpl
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

code = 'deterministic_shallow_demo'
data = Dataset('results/'+code+'/data.nc', 'r')

# plot details
ms = 15
lw = 2
fs = 40
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = {'size':fs}
plt.rc('font',**font)

fig, axs = plt.subplots(6, 1, figsize=(15, 15), sharex=True, gridspec_kw = {'hspace':0.06})

tmax = 20
dt = 0.001
dt_dump = tmax / 2000
timesteps = [0, 300, 600, 900, 1200, 1500]
x = data['x'][:]

for i, (timestep, ax) in enumerate(zip(timesteps, axs)):

    u = data['u_field'][timestep,:]
    xmin = data['min_du_loc'][timestep]
    xmax = data['max_du_loc'][timestep]

    for j, x_value in enumerate(x):
        if xmin > x_value:
            xmin_idx = j
        if xmax > x_value:
            xmax_idx = j

    # for plotting turn them into arrays
    xmin = [data['min_du_loc'][timestep], data['min_du_loc'][timestep]]
    xmax = [data['max_du_loc'][timestep], data['max_du_loc'][timestep]]
    u_xmin = [0, u[xmin_idx]]
    u_xmax = [0, u[xmax_idx]]

    time = timestep * dt_dump

    ax.plot(x, u, color='black', linestyle='-', marker='', ms=ms, lw=lw)
    ax.plot(xmin, u_xmin, color='black', linestyle='--', marker='', ms=ms, lw=lw)
    ax.plot(xmax, u_xmax, color='black', linestyle='--', marker='', ms=ms, lw=lw)

    ax.set_ylim([-0.05,0.6])
    ax.set_xlim([5, 25])
    ax.text(0.98, 0.9, r'$t$ = %.1f s' % time,
            horizontalalignment='right',
            verticalalignment='top',
            transform=ax.transAxes)

    if i == len(timesteps) - 1:
        ax.set_xlabel(r'$x \ /$ m', fontsize=48)

    fig.text(0.04, 0.5, r'$u \ / $ m s$^{-1}$', fontsize=48,
             ha='center', va='center', rotation='vertical')

fig.subplots_adjust(wspace=0, hspace=0)
plt.tight_layout()

plt.savefig('figures/shallow_demo.png', bbox_inches='tight')
