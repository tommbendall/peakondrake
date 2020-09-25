from peakondrake import *
from netCDF4 import Dataset
from datetime import datetime

code = 'mu_deterministic_standard_dt'
Ld = 40.
tmax = 20
dt = 0.001

starttime = datetime.now()

experiment(code, Ld, tmax,
           resolutions=[20000, 16000, 12500, 10000, 7500, 5000, 3000, 2000, 1500, 1000],
           dts=dt,
           sigmas=0.0,
           seeds=0,
           schemes='hydrodynamic',
           timesteppings='midpoint',
           ics='proper_peak',
           num_Xis=1,
           Xi_family='sines',
           alphasq=1.0,
           c0=0.,
           gamma=0.,
           diagnostics=['mu', 'max_du', 'min_du'],
           fields_to_output=['du'],
           ndump=int(tmax / (2000 * dt)),
           field_ndump=int(tmax / (1 * dt)),
           allow_fail=True,
           nXi_updates=1)

endtime = datetime.now()
print(code)
print('Total runtime was %i' % (endtime - starttime).seconds)
