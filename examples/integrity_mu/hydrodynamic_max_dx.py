from peakondrake import *
from netCDF4 import Dataset
from datetime import datetime

code = 'integrity_hydrodynamic_max_dx'
Ld = 40.
tmax = 50
dt = 0.001

starttime = datetime.now()

experiment(code, Ld, tmax,
           resolutions=[2000, 800, 500, 200, 100],
           dts=dt,
           sigmas=[0.0, 0.2],
           seeds=0,
           schemes='hydrodynamic',
           timesteppings='midpoint',
           ics='one_peak',
           num_Xis=3,
           Xi_family='sines',
           alphasq=1.0,
           c0=0.,
           gamma=0.,
           diagnostics=['mu'],
           fields_to_output=['du'],
           ndump=int(tmax / (1000 * dt)),
           field_ndump=int(tmax / (10 * dt)),
           allow_fail=True,
           nXi_update=1)

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
