import os
import matplotlib as mpl
if os.environ.get('DISPLAY', '') == '':
    print('no display found. Using the non-interactive Agg backend')
    mpl.use('Agg')
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

generate_data = False

# plot details
ms = 15
lw = 3
fs = 40
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = {'size':fs}
plt.rc('font',**font)
colors = ['coral', 'maroon', 'purple', 'navy', 'royalblue']
linestyles = ['solid', 'dashed', 'dashdot', (0, (3, 1, 1, 1, 1, 1)), 'dotted']


base_code = 'final_mu_experiment'
codes = [base_code+'_'+str(i) for i in range(10)]

# get times to output mu histograms at
# to do this, open a first data file
# define the time indices to get the t_outs
first_data = Dataset('src/peakondrake/results/'+codes[0]+'/data.nc', 'r')
times = first_data['time'][:]
t_out_idxs = [int(tidx) for tidx in np.linspace(0, len(times)/2, num=5)]
t_outs = [times[tidx] for tidx in t_out_idxs]

# make an array for storing mus for histograms
dxs = first_data['deltax'][:]
seeds = first_data['seed'][:]
total_seeds = len(seeds) * len(codes)
mu_hist_data = np.zeros((len(dxs), len(t_outs), total_seeds))

# make things for peakon formation time
threshold = 0.4
peakon_formation_times = np.full(total_seeds, np.nan)
did_not_form = []
total_seed_counter = 0
individual_seeds_to_plot = [3, 15, 20]
unfiltered = False

def smallest_mus(data, idx, seed=None):
    """
    Pick the lowest mu for each of the recorded mus.

    :arg data: the netCDF4 Dataset.
    :arg idx: the deltax index.
    :arg seed: the seed.
    """

    time = data['time'][:]
    out_mus = []
    for tidx, t in enumerate(time):
        if seed is not None:
            now_mus = [data['mu_'+str(i)][tidx, idx, seed] for i in range(3)]
        else:
            now_mus = [data['mu_'+str(i)][tidx, idx] for i in range(3)]
        out_mus.append(np.nanmin(np.array(now_mus)))

    return out_mus


def low_pass_filter(t, f, tau):
    """
    Performs a low pass filter of f, where tau is the time scale below which to discard.

    solves g' = (1/tau)(g - f)
    g_n1 = g_n + (dt / tau) * (f_n1 - g_n1)

    :arg t: array of times.
    :arg f: array of function f.
    :arg tau: the time scale to filter below.
    """

    g = []
    g.append(f[0])
    for i, (tn, tn1, fn1) in enumerate(zip(t, t[1:], f[1:])):
        dt = tn1 - tn
        gn = g[-1]
        gn1 = (gn + (dt/tau)*fn1) / (1 + (dt/tau))
        g.append(gn1)

    return g

def make_plot_for_individual_seeds(filtered_mus, times, dxs, base_code, seed=None, peakon_time=None, best_fits=None, colors=None, unfiltered=False, linestyles=None):

    """
    Makes a plot of mu against time, with lines for mu at each dx.

    :arg filtered_mus: a list of the mu values.
    :arg times: a list of the times at which values are recorded.
    :arg dxs: a list of the different dx values for each run.
    :arg base_code: the base code for the folder to save the image to.
    :arg seed: for labelling the image. Otherwise it is assumed to be deterministic.
    :arg peakon_time: the time at which a peakon was decided to have formed. None if not formed.
    :arg best_fits: the gradients of best fit lines used to determine peakon formation.
    :arg colors: colors for the mu lines to have.
    """

    fig = plt.figure(1, figsize=(12,8))
    ax = fig.add_subplot(111)
    plt.subplots_adjust(left=0.15, right=0.9, bottom=0.17)

    for filtered_mu, color, dx, linestyle in zip(filtered_mus, colors, dxs, linestyles):
        ax.plot(times, filtered_mu, color=color, label=r'$\Delta x = %.1e$' % dx, linestyle=linestyle)

    if peakon_formed:
        ax.plot([peakon_time, peakon_time], [0, np.max(filtered_mus)], color='black', linestyle='--', marker='')

    if best_fits is not None:
        ax.plot(times, best_fits, color='grey', label='best fit gradient', marker='')

    ax.set_xlabel(r'$t$')
    ax.set_xlim([0, 60])
    ax.set_ylabel(r"$\mu$")
    ax.grid('on')
    handles, labels = ax.get_legend_handles_labels()
    #lgd = fig.legend(handles, labels, fontsize=20, loc='upper right', edgecolor='black', bbox_to_anchor=(1.0, 1.0))
    lgd = fig.legend(handles, labels, fontsize=24, loc='upper right', bbox_to_anchor=(1.0, 1.0))

    figcode = 'figures/final/mu_with_time_'
    if seed is not None:
        figcode += str(seed)
    else:
        figcode += 'deterministic'
    if unfiltered:
        figcode += '_unfiltered'

    plt.savefig(figcode+'.png', bbox_extra_artists=(lgd,))
    plt.close()

starttime = datetime.now()

# deal with deterministic code first
data = Dataset('src/peakondrake/results/final_mu_deterministic/data.nc', 'r')

deterministic_times = data['time'][:]
deterministic_t_out_idxs = [int(tidx) for tidx in np.linspace(0, len(deterministic_times)/2, num=5)]
deterministic_dxs = data['deltax'][:]
mu_hist_deterministic = np.zeros((len(deterministic_dxs), len(deterministic_t_out_idxs)))
peakon_formed = False
peakon_time = None
mus_for_this_seed = []
deterministic_peakon_time = None

for i, dx in enumerate(deterministic_dxs):

    best_mu = smallest_mus(data, i)
    filtered_mu = low_pass_filter(deterministic_times, best_mu, 0.3)
    chosen_data = best_mu if unfiltered else filtered_mu

    for j, t_out_idx in enumerate(deterministic_t_out_idxs):
        mu_hist_deterministic[i, j] = chosen_data[t_out_idx]

    mus_for_this_seed.append(chosen_data)

# go through each time finding gradient
best_fits = []
for t_idx, time in enumerate(deterministic_times):
    # extract values at different dxs for this time
    filtered_mus_at_t = [filtered_mu[t_idx] for filtered_mu in mus_for_this_seed]
    best_fit = np.poly1d(np.polyfit(np.log(deterministic_dxs), np.log(filtered_mus_at_t), deg=1))
    best_fits.append(best_fit[1])

    if best_fit[1] > threshold and not peakon_formed:
        peakon_formed = True
        deterministic_peakon_time = time

    make_plot_for_individual_seeds(mus_for_this_seed, deterministic_times, deterministic_dxs, base_code, seed=None, peakon_time=deterministic_peakon_time, best_fits=None, colors=colors, unfiltered=unfiltered, linestyles=linestyles)


# now loop through actual files
if generate_data:
    for code in codes:

        data = Dataset('src/peakondrake/results/'+code+'/data.nc', 'r')

        for seed in seeds:
            print('SEED', total_seed_counter)
            peakon_formed = False
            peakon_time = None
            mus_for_this_seed = []

            for i, dx in enumerate(dxs):

                best_mu = smallest_mus(data, i, seed)
                filtered_mu = low_pass_filter(times, best_mu, 0.3)
                chosen_data = best_mu if unfiltered else filtered_mu

                for j, t_out_idx in enumerate(t_out_idxs):
                    mu_hist_data[i, j, total_seed_counter] = chosen_data[t_out_idx]

                mus_for_this_seed.append(chosen_data)

                # go through each time finding gradient
                best_fits = []
            for t_idx, time in enumerate(times):
                # extract values at different dxs for this time
                filtered_mus_at_t = [filtered_mu[t_idx] for filtered_mu in mus_for_this_seed]
                best_fit = np.poly1d(np.polyfit(np.log(dxs), np.log(filtered_mus_at_t), deg=1))
                best_fits.append(best_fit[1])

                if best_fit[1] > threshold and not peakon_formed:
                    peakon_formation_times[total_seed_counter] = time
                    peakon_formed = True
                    peakon_time = time

            if total_seed_counter in individual_seeds_to_plot:
                make_plot_for_individual_seeds(mus_for_this_seed, times, dxs, base_code, seed=total_seed_counter, peakon_time=peakon_time, best_fits=best_fits, colors=colors, unfiltered=unfiltered, linestyles=linestyles)

            if not peakon_formed:
                did_not_form.append(total_seed_counter)
            total_seed_counter += 1

        data.close()

else:
    # now we are picking up data
    hist_file = Dataset('hist_data.nc', 'r')
    peakon_formation_times = hist_file['peakon_time'][:]
    mu_hist_data = hist_file['mu'][:,:,:]

# make histogram plot of peakon formation time
fig = plt.figure(1, figsize=(8,8))
ax = fig.add_subplot(111)
plt.subplots_adjust(left=0.1, right=0.9, bottom=0.15)
heights, bins, patches = ax.hist(peakon_formation_times, bins=20, color='grey')
ax.plot([deterministic_peakon_time, deterministic_peakon_time], [0, np.max(heights)], color='black', linestyle='--')
ax.set_xlabel(r'Time of peakon formation')
ax.grid('on')
figcode = 'figures/final/peakon_formation_hist'
if unfiltered:
    figcode += '_unfiltered'
plt.savefig(figcode+'.png')


font = {'size':30}
plt.rc('font',**font)
# make grid of histograms of mu
fig, axs = plt.subplots(len(t_outs), len(dxs), figsize=(4*len(dxs), 15))
fig.subplots_adjust(wspace=0.1, hspace=0.4)

dx_offset = len(deterministic_dxs) - len(dxs)

for t_idx, t_out in enumerate(t_outs):
    for i, dx in enumerate(dxs):
        # make mu hist list so as to cut out any nans
        mu_hist_list = []
        for mu in mu_hist_data[i,t_idx,:]:
            if not np.isnan(mu):
                mu_hist_list.append(mu)

        axs[t_idx, i].hist(mu_hist_list, bins=np.linspace(0, 2, 40), color='grey')
        axs[t_idx, i].plot([mu_hist_deterministic[i+dx_offset, t_idx], mu_hist_deterministic[i+dx_offset,t_idx]], [0, total_seeds], color='black', linestyle='--')
        if i == 0:
            axs[t_idx,i].set_ylabel(r'$t = $ %0.f' % t_out)
        if t_idx == 0:
             axs[t_idx,i].set_title(r'$\Delta x = $ %.1e' % dx)
        if t_idx == len(t_outs) - 1:
             axs[t_idx,i].set_xlabel(r"$\mu$")

for ax in axs.flat:
    ax.set_xlim([0, 2])
    ax.set_ylim([0, total_seeds])
    ax.label_outer()

figcode = 'figures/final/histograms'
if unfiltered:
    figcode += '_unfiltered'
plt.savefig(figcode+'.png')

print('Runs that did not form peakons were', did_not_form)

# now write out data to a file so that this data can quickly be picked up
if generate_data:
    hist_file = Dataset('hist_data.nc', 'w')
    hist_file.createDimension('seed', total_seeds)
    hist_file.createDimension('deltax', len(dxs))
    hist_file.createDimension('time', len(t_outs))
    hist_file.createVariable('seed', float, ('seed',))
    hist_file.createVariable('deltax', float, ('deltax',))
    hist_file.createVariable('time', float, ('time',))
    hist_file.createVariable('mu', float, ('deltax', 'time', 'seed'))
    hist_file.createVariable('peakon_time', float, ('seed',))

    hist_file['seed'][:] = range(total_seeds)
    hist_file['deltax'][:] = dxs
    hist_file['time'][:] = t_outs
    hist_file['mu'][:,:,:] = mu_hist_data
    hist_file['peakon_time'][:] = peakon_formation_times
    hist_file.close()

endtime = datetime.now()
print('Total runtime was %i' % (endtime - starttime).seconds)
