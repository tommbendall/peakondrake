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


for mode in ['shallow']:

    code = mode+'_deterministic_diagnostics'
    data = Dataset('results/'+code+'/data.nc', 'r')
    times = data['time'][:]

    for diagnostic in ['peakon_nu']:


        fig = plt.figure(figsize=(15,5))
        # Use grid spec to get subplots within a subplot
        outer_grid = gridspec.GridSpec(1, 2, wspace=0.36, width_ratios=[1,2])
        left_inner_grid = gridspec.GridSpecFromSubplotSpec(3, 1,
                                                           subplot_spec=outer_grid[0],
                                                           wspace=0.05, hspace=0.05)


        # First concentrate on plotting the time series
        # Ensure that axes are shared
        ax2 = fig.add_subplot(left_inner_grid[2])
        ax0 = fig.add_subplot(left_inner_grid[0], sharex=ax2)
        ax1 = fig.add_subplot(left_inner_grid[1], sharex=ax2)

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
            if diagnostic == 'peakon_nu':
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

            if i == 1:
                ax.set_ylabel(r'$\nu \ / $ m$^{-1}$ s$^{-1}$')

            if i == 2:
                ax.set_xlabel(r'$t \ / $ s')

        # Plot the time series on the same subplot
        right_inner_grid = gridspec.GridSpecFromSubplotSpec(1, 2,
                                                            subplot_spec=outer_grid[1],
                                                            wspace=0.08, hspace=0.05)
        ax0 = fig.add_subplot(right_inner_grid[0])
        ax1 = fig.add_subplot(right_inner_grid[1], sharey=ax0)



        dnu = [np.poly1d(np.polyfit(dxs, data[diagnostic][t,:], deg=1))[1] for t in range(len(times))]
        ax0.axvspan(15, 20, color='lightgrey')
        ax0.plot(times, dnu, color='black')
        ax0.set_xlabel(r'$t \ / $ s')

        ax1.hist(dnu[749:], orientation='horizontal', color='dimgrey')
        mean_dnu = np.mean(dnu[749:])
        ax1.plot([0,100], [mean_dnu,mean_dnu], color='black', linestyle='--')
        ax1.set_xlim([0, 100])
        ax1.set_xlabel('Num. points')

        # Hide inner tick labels
        plt.setp(ax1.get_yticklabels(), visible=False)



        fig.text(0.355, 0.5, r'$d\nu / d(\Delta x) \ / $ m$^{-2}$ s$^{-1}$', va='center', rotation='vertical')


        plt.savefig('figures/'+mode+'_deterministic_'+diagnostic+'.png', bbox_inches='tight')
        plt.close()
