from peakondrake import *
from netCDF4 import Dataset
from datetime import datetime

base_code = 'integrity_quiet_loud'
Ld = 40.
tmax = 20

starttime = datetime.now()

dts = [0.0005, 0.001, 0.002]

for i, dt in enumerate(dts):

    code = base_code + '_'+str(i)

    experiment(code, Ld, tmax,
               resolutions=1000,
               dts=dt,
               sigmas=[0.002, 0.2],
               seeds=range(20),
               schemes='hydrodynamic',
               timesteppings='midpoint',
               ics='one_peak',
               num_Xis=1,
               Xi_family='sines',
               alphasq=1.0,
               c0=0.,
               gamma=0.,
               diagnostics=['mu'],
               fields_to_output=['du'],
               ndump=int(tmax / (1000 * dt)),
               field_ndump=int(tmax / (1 * dt)),
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