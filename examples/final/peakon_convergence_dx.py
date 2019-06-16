from peakondrake import *
from netCDF4 import Dataset
from datetime import datetime

code = 'final_peakon_convergence_dx'
Ld = 40.
tmax = 200
dt = 0.0001

starttime = datetime.now()

experiment(code, Ld, tmax,
           resolutions=[400, 1000, 2500, 5000, 10000],
           dts=dt,
           sigmas=[0.02],
           seeds=0,
           schemes='hydrodynamic',
           timesteppings='midpoint',
           ics='peakon',
           num_Xis=1,
           Xi_family='sines',
           alphasq=1.0,
           c0=0.,
           gamma=0.,
           diagnostics=['p_pde', 'q_pde'],
           fields_to_output=[],
           ndump=int(tmax / (1000 * dt)),
           field_ndump=int(tmax / (1 * dt)),
           allow_fail=True,
           peakon_equations=True)

endtime = datetime.now()
print(code)
print('Total runtime was %i' % (endtime - starttime).seconds)

data = Dataset('results/'+code+'/data.nc', 'r')
runtimes = []


for j, k in zip(initial_times, final_times):
    runtimes.append((datetime.fromtimestamp(k) - datetime.fromtimestamp(j)).seconds)

print('runtimes', runtimes)
print('failed times', failed_times)
