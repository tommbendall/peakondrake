from peakondrake import *
from netCDF4 import Dataset
from datetime import datetime

base_code = 'integrity_conforming_dx_with_dt'
Ld = 40.
tmax = 50

starttime = datetime.now()

dts = [0.0001, 0.002, 0.001, 0.05, 0.01]
dx = [10000, 2000, 1000, 500, 100]

for i, dt in enumerate(dts):

    code = base_code+'_'+str(i)

    experiment(code, Ld, tmax,
               resolutions=10000,
               dts=dt,
               sigmas=[0.0, 0.2],
               seeds=0,
               schemes='conforming',
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
               allow_fail=True)

endtime = datetime.now()
print(base_code)
print('Total runtime was %i' % (endtime - starttime).seconds)

for i, dt in enumerate(dts):
    code = base_code + '_'+str(i)
    data = Dataset('results/'+code+'/data.nc', 'r')
    runtimes = []
    failed_times = data['failed_time'][:]

    for j, k in zip(data['wallclock_time'][0,:], data['wallclock_time'][1,:]):
        runtimes.append((datetime.fromtimestamp(k) - datetime.fromtimestamp(j)).seconds)

    print('dt', i, 'runtimes', runtimes)
    print('dt', i, 'failed times', failed_times)
    print('')
