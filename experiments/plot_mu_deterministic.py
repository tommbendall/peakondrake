import os
import math
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from datetime import datetime

code = 'mu_deterministic_standard_dt'
data = Dataset('results/'+code+'/data.nc', 'r')

def smallest_mus(data, dx_idx):
    """
    Pick the lowest mu for each of the recorded mus.

    :arg data: the netCDF4 Dataset.
    :arg idx: the deltax index.
    :arg sigma: the strength of the noise.
    :arg seed: the seed.
    """

    time = data['time'][:]
    out_mus = []
    for tidx, t in enumerate(time):
        now_mus = [data['mu_'+str(i)][tidx, dx_idx] for i in range(3)]
        out_mus.append(np.nanmin(np.array(now_mus)))

    return out_mus

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


# plot details
ms = 15
lw = 3
fs = 32
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = {'size':fs}
plt.rc('font',**font)

fig = plt.figure(figsize=(15,8))
# Use grid spec to get subplots within a subplot
outer_grid = gridspec.GridSpec(1, 2, wspace=0.25)
inner_grid = gridspec.GridSpecFromSubplotSpec(2, 2,
                                              subplot_spec=outer_grid[0],
                                              wspace=0.05, hspace=0.05)

tmax = 20
Ld = 40
times = data['time'][:]
resolutions = [20000, 16000, 12500, 10000, 7500, 5000, 3000, 2000, 1000]
res_idx_to_plot = [8, 5, 3, 0]
colors = ['red', 'purple', 'blue', 'black']
soft_colors = ['salmon', 'plum', 'lightsteelblue', 'gainsboro']
markers = ['s', '^', 'o', 'x']

# First concentrate on plotting the time series
# Ensure that axes are shared
ax2 = fig.add_subplot(inner_grid[2])
ax0 = fig.add_subplot(inner_grid[0], sharex=ax2, sharey=ax2)
ax3 = fig.add_subplot(inner_grid[3], sharex=ax2, sharey=ax2)
ax1 = fig.add_subplot(inner_grid[1], sharex=ax3, sharey=ax0)

# Hide inner tick labels
plt.setp(ax0.get_xticklabels(), visible=False)
plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax1.get_yticklabels(), visible=False)
plt.setp(ax3.get_yticklabels(), visible=False)

axs = [ax0, ax1, ax2, ax3]

for i, res_idx in enumerate(res_idx_to_plot):

    ax = axs[i]

    dx = Ld / resolutions[res_idx]

    best_mu = smallest_mus(data, res_idx)

    nu = (data['max_du'][:,res_idx] - data['min_du'][:,res_idx]) / best_mu[:]

    ax.plot(times, nu, color=colors[i], linestyle='-')
    ax.text(0.05, 0.95, r'$\Delta x = \ $'+sci_notation(dx),
            horizontalalignment='left',
            verticalalignment='top',
            transform=ax.transAxes,
            fontsize=20)


    # ax.set_xlabel(r'$t \ / $ s')
    # ax.set_ylabel(r'$\nu \ / $ m$^{-1}$ s$^{-1}$')
    # fig.add_subplot(ax)


# Plot the time series on the same subplot
ax = plt.Subplot(fig, outer_grid[1])
for i, res_idx in enumerate(res_idx_to_plot):

    dx = Ld / resolutions[res_idx]

    best_mu = smallest_mus(data, res_idx)

    nu = (data['max_du'][:,res_idx] - data['min_du'][:,res_idx]) / best_mu[:]

    ax.plot(times, nu, color=soft_colors[i], linestyle='-')
    # ax.set_xlabel(r'$t \ / $ s')
    # ax.set_ylabel(r'$\nu \ / $ m$^{-1}$ s$^{-1}$')

    # Determine peakon formation time
    nu_jump = nu[1:] - nu[:-1]
    biggest_jump_time = np.nan

    # Find formation time
    for n in range(2,len(times)-1):

        if (nu_jump[n-1]) > (nu_jump[n-2]):
            # is the jump bigger than the next average jumps?
            if nu_jump[n-1] > np.mean(np.sqrt(nu_jump[n:]**2)):
                # time of jump is mu_jump[j] = mu[j] - mu[j-1]
                biggest_jump_time = 0.5*(times[n] + times[n-1])
                biggest_jump_index = n
                break

    # Best fit before peakon jump
    before_times = np.linspace(0, biggest_jump_time, 5)
    best_fit_before = np.poly1d(np.polyfit(times[0:n], nu[0:n], deg=3))
    ax.plot(times[0:n], best_fit_before(times[0:n]), color=colors[i], linestyle='-', marker='')
    ax.plot(before_times, best_fit_before(before_times), color=colors[i], linestyle='', marker=markers[i], label=r'$\Delta x = \ $'+sci_notation(dx), ms=ms)

    # Best fit after peakon jump
    after_times = np.linspace(biggest_jump_time, times[-1], 5)
    best_fit_after = np.poly1d(np.polyfit(times[n+1:], nu[n+1:], deg=4))
    ax.plot(times[n+1:], best_fit_after(times[n+1:]), color=colors[i], linestyle='-', marker='')
    ax.plot(after_times, best_fit_after(after_times), color=colors[i], linestyle='', marker=markers[i], ms=ms)

    # Join up the before and after
    ax.plot(times[n-1:n+1], [best_fit_before(times[n-1]), best_fit_after(times[n])], color=colors[i], marker='')

    handles, labels = ax.get_legend_handles_labels()

fig.add_subplot(ax)
fig.text(0.5, 0.0, r'$t \ / $ s', ha='center')
fig.text(0.04, 0.5, r'$\nu \ / $ m$^{-1}$ s$^{-1}$', va='center', rotation='vertical')
leg = fig.legend(handles, labels, ncol=4, loc='upper center', fontsize=24)

# Do lines of best fit before and after peakon formation
# Use colors and symbols (rather than soft_colors)

plt.savefig('figures/nu_deterministic_grid.png', bbox_inches='tight')
plt.close()

# Now plot all peakon formation times
dxs = data['deltax'][:-1]
peakon_formation_times = []
for dx_idx, dx in enumerate(dxs):

    best_mu = smallest_mus(data, dx_idx)

    nu = (data['max_du'][:,dx_idx] - data['min_du'][:,dx_idx]) / best_mu[:]

    # Determine peakon formation time
    nu_jump = nu[1:] - nu[:-1]
    biggest_jump_time = np.nan

    # Find formation time
    for n in range(2,len(times)-1):

        if (nu_jump[n-1]) > (nu_jump[n-2]):
            # is the jump bigger than the next average jumps?
            if nu_jump[n-1] > np.mean(np.sqrt(nu_jump[n:]**2)):
                # time of jump is mu_jump[j] = mu[j] - mu[j-1]
                biggest_jump_time = 0.5*(times[n] + times[n-1])
                break

    peakon_formation_times.append(biggest_jump_time)

fig = plt.figure(1, figsize=(8,6))
ax = fig.add_subplot(111)
ax.plot(dxs, peakon_formation_times, marker='x', color='black', linestyle='', ms=ms)
fit_dxs = np.linspace(0.0, np.max(dxs), 20)
best_fit = np.poly1d(np.polyfit(dxs, peakon_formation_times, deg=3))
ax.plot(fit_dxs, best_fit(fit_dxs), color='black', linestyle='-', marker='')

ax.set_xlabel(r'$\Delta x \ / $ m')
ax.set_ylabel(r'$\tau \ / $ s')
ax.set_xlim([0, np.max(dxs)*1.05])
plt.savefig('figures/deterministic_peakon_formation_times.png', bbox_inches='tight')
