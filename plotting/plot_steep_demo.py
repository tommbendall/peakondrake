import os
import matplotlib as mpl
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

code = 'deterministic_steep_demo'
data = Dataset('results/'+code+'/data.nc', 'r')

# plot details
ms = 15
lw = 2
fs = 40
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = {'size':fs}
plt.rc('font',**font)

fig, axs = plt.subplots(5, 2, figsize=(15, 15), sharex=True, sharey=True, gridspec_kw = {'wspace':0.08, 'hspace':0.02})

tmax = 20
dt = 0.001
dt_dump = tmax / 2000
timesteps = [0, 300, 50, 400, 100, 500, 150, 600, 200, 800] # the funny ordering is to match the ordering of the axs
x = data['x'][:]

for i, (timestep, ax) in enumerate(zip(timesteps, axs.flatten())):

    u = data['u_field'][timestep,:]

    time = timestep * dt_dump

    ax.plot(x, u, color='black', linestyle='-', marker='', ms=ms, lw=lw)

    ax.set_ylim([-0.3,0.7])
    ax.set_xlim([5, 20])
    if (i % 2 == 0):
        ax.text(0.98, 0.9, r'$t$ = %.1f s' % time,
                horizontalalignment='right',
                verticalalignment='top',
                transform=ax.transAxes)
    else:
        ax.text(0.05, 0.9, r'$t$ = %.1f s' % time,
                horizontalalignment='left',
                verticalalignment='top',
                transform=ax.transAxes)

    fig.text(0.5, 0.05, r'$x \ /$ m', fontsize=48,
             ha='center', va='center')

    if i == 4:
        ax.set_ylabel(r'$u \ / $ m s$^{-1}$', fontsize=48)

fig.subplots_adjust(wspace=0, hspace=0)
plt.tight_layout()

plt.savefig('figures/steep_demo.png', bbox_inches='tight')
