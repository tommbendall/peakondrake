import os
import matplotlib as mpl
if os.environ.get('DISPLAY', '') == '':
    print('no display found. Using the non interactive Agg backend')
    mpl.use('Agg')
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

ms = 15
lw = 3
fs = 30

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
font = {'size':fs}
plt.rc('font',**font)
plt.locator_params(nbins=3)

code = 'periodic_peakon_convergence_dx'
data = Dataset('results/'+code+'/data.nc', 'r')
dxs = []
ncells = [500., 750., 1000., 1500., 2000.]
Ld = 40.0

fig = plt.figure(1, figsize=(8,8))
ax = fig.add_subplot(111)
plt.subplots_adjust(left=0.22, right=0.9, bottom=0.15)

error = []
dts = []
best_fit = True

for i, ncell in enumerate(ncells):
    dxs.append(Ld/ncell)
    error.append(data['u_error_with_sde'][-1,i])

ax.plot(np.log(dxs), np.log(error), color='black', linestyle='', marker='+', ms=ms)

if best_fit:
    p = np.poly1d(np.polyfit(np.log(dxs), np.log(error), deg=1))
    ax.plot(np.log(dxs), np.poly1d(np.polyfit(np.log(dxs), np.log(error), deg=1))(np.log(dxs)), linestyle='-', color='black', label=r'gradient = %.2f' % p[1])

handles, labels = ax.get_legend_handles_labels()
leg = ax.legend(handles, labels, loc='upper left')

ax.set_xlabel(r'$\log(\Delta x)$')
ax.set_ylabel(r'$\log\left(||u_{PDE}-u_{SDE}||\right)$')
plt.locator_params(nbins=3)

plt.savefig('figures/periodic_peakon_dx_convergence.png', bbox_extra_artists=(leg,), bbox_inches='tight')
