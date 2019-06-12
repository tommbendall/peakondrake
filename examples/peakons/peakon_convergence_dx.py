from peakondrake import *
from netCDF4 import Dataset
from datetime import datetime

code = 'peakon_dx'
Ld = 40.
tmax = 200
dt = 0.001

starttime = datetime.now()

experiment(code, Ld, tmax,
           resolutions=[400, 800, 1000, 1500, 2000],
           dts=dt,
           sigmas=[0.0, 0.2],
           seeds=0,
           schemes='hydrodynamic',
           timesteppings='midpoint',
           ics='one_peak',
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
           nXi_update=1,
           peakon_equations=True)

endtime = datetime.now()
print(code)
print('Total runtime was %i' % (endtime - starttime).seconds)

data = Dataset('results/'+code+'/data.nc', 'r')
runtimes = []

failed_times = data['failed_time'][:]
initial_times = [initial_time for sublist in data['wallclock_time'][0,:,:] for initial_time in sublist]
final_times = [final_time for sublist in data['wallclock_time'][1,:,:] for final_time in sublist]

for j, k in zip(initial_times, final_times):
    runtimes.append((datetime.fromtimestamp(k) - datetime.fromtimestamp(j)).seconds)

print('runtimes', runtimes)
print('failed times', failed_times)
