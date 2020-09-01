import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

num_data_files = 20
base_data_name = 'results/new_mu_stochastic_flat_'
sigmas = [0.05, 0.1, 0.2, 0.5, 1.0]
resolutions = [1000, 1500, 2000, 3000, 5000]
num_seeds_per_file = 5

# Create a netcdf file for storing all of the peakon formation data
formation_data_name = 'results/peakon_formation/formation_times.nc'
formation_data = Dataset(formation_data_name, 'w')
formation_data.createDimension('resolution', len(resolutions))
formation_data.createDimension('sigma', len(sigmas))
formation_data.createDimension('seed', num_data_files*num_seeds_per_file)
formation_data.createVariable('resolution', int, ('resolution',))
formation_data.createVariable('sigma', float, ('sigma',))
formation_data.createVariable('seed', int, ('seed',))
formation_data.createVariable('formation_time', float, ('resolution', 'sigma', 'seed'))
formation_data.close()

def smallest_mus(data, dx_idx, sigma_idx, seed_idx):
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
        now_mus = [data['mu_'+str(i)][tidx, dx_idx, seed_idx, sigma_idx] for i in range(3)]
        out_mus.append(np.nanmin(np.array(now_mus)))

    return out_mus


# Loop through all data files
for i in range(num_data_files):
    data_name = base_data_name+str(i)+'/data.nc'
    data = Dataset(data_name, 'r')
    time = data['time'][:]

    for j, res in enumerate(resolutions):
        for k, sigma in enumerate(sigmas):
            for l in range(num_seeds_per_file):
                seed = 5 * i + l

                print(i, j, k, l)

                # Find the best mu data
                best_mu = smallest_mus(data, j, k, l)

                # Make actual diagnostic
                nu = (data['max_du'][:,j,l,k] - data['min_du'][:,j,l,k]) / best_mu[:]
                nu_jump = nu[1:] - nu[:-1]

                # Find formation time
                for n in range(1,len(time)-1):
                    # is the jump bigger than the previous jump?
                    biggest_jump_time = np.nan
                    if (nu_jump[n]) > (nu_jump[n-1]):
                        # is the jump bigger than the next average jumps?
                        if nu_jump[n] > np.mean(np.sqrt(nu_jump[n+1:]**2)):
                            # time of jump is mu_jump[j] = mu[j] - mu[j-1]
                            print(n)
                            biggest_jump_time = 0.5*(time[n] + time[n-1])
                            break

                # Put formation time into netcdf data
                formation_data = Dataset(formation_data_name, 'a')
                formation_data['formation_time'][j,k,l] = biggest_jump_time
                formation_data.close()

                # Make plot of diagnostic so we can check if it's good
                formation_times_to_plot = [biggest_jump_time, biggest_jump_time]
                jump_values = [0, np.nanmax(nu)]
                figcode = 'figures/peakon_diagnostic_flat_res'+str(j)+'_sigma'+str(k)+'_seed'+str(l)

                fig = plt.figure(1, figsize=(8,8))
                ax = fig.add_subplot(111)
                ax.plot(time, nu)
                ax.plot(formation_times_to_plot, jump_values, color='black', linestyle='--')
                plt.savefig(figcode+'.png')
                plt.close()

    data.close()
