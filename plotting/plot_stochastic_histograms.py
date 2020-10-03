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

sigmas = [0.05, 0.2, 0.5, 1.0, 2.0]
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

    data = Dataset('results/isambard_peakon_formation/'+mode+'_flat_formation_times.nc', 'r')

    fig, axs = plt.subplots(5, 2, figsize=(15, 15), sharex='col')

    for m, diagnostic in enumerate(['peakon_measure', 'wave_breaking_measure']):

        # max_min = np.nanmax(data[diagnostic][:,:]) - np.nanmin(data[diagnostic][:,:])
        # ylim_upper = np.nanmax(data[diagnostic][:,:]) + 0.02 * max_min
        # ylim_lower = np.nanmin(data[diagnostic][:,:]) - 0.02 * max_min


        for i, sigma in enumerate(sigmas):

            ax = axs[i,m]

            values = data[diagnostic][i,:]

            # make bins wider so we can actually see results
            if (i in [0, 1]) and m == 1 and mode == 'shallow':
                ax.hist(values, bins=range(-10,10), color='grey')
            else:
                ax.hist(values, color='grey')

            # add text
            if m == 0:
                ax.text(0.05, 0.95, r'$\xi = %.2f $' % sigma,
                        horizontalalignment='left',
                        verticalalignment='top',
                        transform=ax.transAxes,
                        fontsize=20)
            else:
                ax.text(0.95, 0.95, r'$\xi = %.2f $' % sigma,
                        horizontalalignment='right',
                        verticalalignment='top',
                        transform=ax.transAxes,
                        fontsize=20)

            # add labels
            if i == 4:
                if m == 0:
                    ax.set_xlabel(r'$\Pi \ / $ m$^{-2}$ s$^{-1}$')
                elif m == 1:
                    ax.set_xlabel(r'$\omega \ / $ m$^{-1}$ s$^{-1}$')
            if i == 2:
                if m == 0:
                    ax.set_ylabel(r'Number of realisations')

            # set limits so that plots can be compared
            if m == 1:
                ax.set_xlim([-25, 150])
            elif m == 0:
                ax.set_xlim([-300,100])

    plt.savefig('figures/'+mode+'_stochastic_histograms.png', bbox_inches='tight')
    plt.close()
