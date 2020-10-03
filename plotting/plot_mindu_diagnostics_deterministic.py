import os
import math
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from datetime import datetime


# plot details
ms = 15
lw = 3
fs = 32
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = {'size':fs}
plt.rc('font',**font)

tmax = 20
Ld = 40
resolutions = [1000, 1500, 2000, 2500, 3000]
dxs = [Ld / res for res in resolutions]
res_idx_to_plot = [0, 2, 4]
colors = ['red', 'purple', 'blue', 'black']
soft_colors = ['salmon', 'plum', 'lightsteelblue', 'gainsboro']
markers = ['s', '^', 'o', 'x']


# Define function for string formatting of scientific notation
def sci_notation(num, decimal_digits=1, precision=0, exponent=None):
    """
    Returns a string representation of the scientific
    notation of the given number formatted for use with
    LaTeX or Mathtext, with specified number of significant
    decimal digits and precision (number of decimal digits
    to show). The exponent to be used can also be specified
    explicitly.
    """
    if exponent is None:
        exponent = int(math.floor(math.log10(abs(num))))
    coeff = round(num / float(10**exponent), decimal_digits)
    if precision is None:
        precision = decimal_digits

    return r"${0:.{2}f}\times 10^{{{1:d}}}$".format(coeff, exponent, precision)


for mode in ['shallow', 'steep']:

    code = mode+'_deterministic_diagnostics'
    data = Dataset('results/'+code+'/data.nc', 'r')
    times = data['time'][:]

    for diagnostic in ['min_du']:


        fig = plt.figure(figsize=(15,5))
        # Use grid spec to get subplots within a subplot
        outer_grid = gridspec.GridSpec(1, 2, wspace=0.35)
        inner_grid = gridspec.GridSpecFromSubplotSpec(3, 1,
                                                      subplot_spec=outer_grid[0],
                                                      wspace=0.05, hspace=0.05)


        # First concentrate on plotting the time series
        # Ensure that axes are shared
        ax2 = fig.add_subplot(inner_grid[2])
        ax0 = fig.add_subplot(inner_grid[0], sharex=ax2)
        ax1 = fig.add_subplot(inner_grid[1], sharex=ax2)

        # Hide inner tick labels
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.setp(ax0.get_xticklabels(), visible=False)

        axs = [ax0, ax1, ax2]

        max_min = np.nanmax(data[diagnostic][:,:]) - np.nanmin(data[diagnostic][:,:])
        ylim_upper = np.nanmax(data[diagnostic][:,:]) + 0.02 * max_min
        ylim_lower = np.nanmin(data[diagnostic][:,:]) - 0.02 * max_min


        for i, res_idx in enumerate(res_idx_to_plot):

            ax = axs[i]

            dx = dxs[res_idx]
            values = data[diagnostic][:,res_idx]

            ax.plot(times, values, color='black', linestyle='-')
            if mode == 'steep':
                ax.text(0.95, 0.05, r'$\Delta x = \ $'+sci_notation(dx),
                        horizontalalignment='right',
                        verticalalignment='bottom',
                        transform=ax.transAxes,
                        fontsize=20)
            else:
                ax.text(0.95, 0.95, r'$\Delta x = \ $'+sci_notation(dx),
                        horizontalalignment='right',
                        verticalalignment='top',
                        transform=ax.transAxes,
                        fontsize=20)
            ax.set_ylim([ylim_lower, ylim_upper])
            ax.set_xlim([times[0], times[-1]])

            if i == 1:
                ax.set_ylabel(r'$\min(u_x) \ / $ s$^{-1}$')

            if i == 2:
                ax.set_xlabel(r'$t \ / $ s')

        # Plot the time series on the same subplot
        ax = plt.Subplot(fig, outer_grid[1])

        dmindu = [np.poly1d(np.polyfit(dxs, data[diagnostic][t,:], deg=1))[1] for t in range(len(times))]
        ax.plot(times, dmindu, color='black')
        max_dmindu = np.nanmax(dmindu)
        ax.plot([times[0], times[-1]], [max_dmindu, max_dmindu], color='black', linestyle='--')
        ax.set_xlabel(r'$t \ / $ s')
        if mode == 'shallow':
            ax.set_ylim([-5, 25])
        else:
            ax.set_ylim([-5, 85])
        ax.set_xlim([times[0], times[-1]])

        fig.add_subplot(ax)
        fig.text(0.48, 0.5, r'$d\min(u_x) / d(\Delta x) \ / $ m$^{-1}$ s$^{-1}$', va='center', rotation='vertical')


        plt.savefig('figures/'+mode+'_deterministic_'+diagnostic+'.png', bbox_inches='tight')
        plt.close()
