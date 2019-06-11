from peakondrake import *
from netCDF4 import Dataset
from datetime import datetime

base_code = 'dump_cost'
Ld = 40.
tmax = 0.02

starttime = datetime.now()

dumps = [1, 10, 100, 1000]
dt = 0.00001

for i, dump in enumerate(dumps):

    code = base_code + '_'+str(i)

    experiment(code, Ld, tmax,
               resolutions=10000,
               dts=dt,
               sigmas=0.0,
               seeds=0,
               schemes='upwind',
               timesteppings='midpoint',
               ics='one_peak',
               num_Xis=3,
               Xi_family='sines',
               alphasq=1.0,
               c0=0.,
               gamma=0.,
               diagnostics=['mu'],
               fields_to_output=['uscalar', 'du'],
               ndump=int(tmax / (dump * dt)),
               field_ndump=int(tmax / (10 * dt)),
               allow_fail=True,
               nXi_update=1)

endtime = datetime.now()
print(base_code)
print('Total runtime was %i' % (endtime - starttime).seconds)

for i, dump in enumerate(dumps):
    code = base_code + '_'+str(i)
    data = Dataset('results/'+code+'/data.nc', 'r')
    runtimes = []
    failed_times = data['failed_time'][:]

    j = data['wallclock_time'][0]
    k = data['wallclock_time'][1]
    runtimes.append((datetime.fromtimestamp(k) - datetime.fromtimestamp(j)).seconds)

    print('dt', i, 'runtimes', runtimes)
    print('dt', i, 'failed times', failed_times)
    print('')
