import os
import matplotlib as mpl
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

base_code = 'convergence_strong_dt_'
all_dts = [0.1, 0.05, 0.025, 0.02, 0.0125, 0.01, 0.00625, 0.005, 0.004, 0.0025, 0.002, 0.00125, 0.001, 0.0008, 0.0005, 0.0004, 0.00025, 0.0002, 0.0001, 0.00005]
colors = ['red', 'blue', 'purple']
Ld = 40

fig = plt.figure(1, figsize=(8,8))
ax = fig.add_subplot(111)
plt.subplots_adjust(left=0.22, right=0.9, bottom=0.25)

error = []
dts = []
best_fit = True

for i, dt in enumerate(all_dts):
    code = base_code+str(i)
    data = Dataset('results/'+code+'/data.nc', 'r')
    dts.append(dt)
    error.append(data['u_error_with_sde'][-1])

ax.plot(np.log(dts), np.log(error), color='black', linestyle='', marker='+', ms=ms)

if best_fit:
    p = np.poly1d(np.polyfit(np.log(dts), np.log(error), deg=1))
    ax.plot(np.log(dts), np.poly1d(np.polyfit(np.log(dts), np.log(error), deg=1))(np.log(dts)), linestyle='-', color='black', label=r'gradient = %.2f' % p[1])

handles, labels = ax.get_legend_handles_labels()
leg = ax.legend(handles, labels, loc='upper left')

ax.set_xlabel(r'$\log(\Delta t)$')
ax.set_ylabel(r'$\log\left(||u_{PDE}-u_{SDE}||\right)$')

plt.savefig('figures/dt_convergence.png')
