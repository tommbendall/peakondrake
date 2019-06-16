from peakondrake import *
from netCDF4 import Dataset
from datetime import datetime

code = 'final_mu_deterministic'
Ld = 40.
tmax = 80
dt = 0.0001

starttime = datetime.now()

experiment(code, Ld, tmax,
           resolutions=[400, 1000, 2500, 5000, 10000],
           dts=dt,
           sigmas=0.0,
           seeds=0,
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
           ndump=int(tmax / (2000 * dt)),
           field_ndump=int(tmax / (1 * dt)),
           allow_fail=True,
           nXi_updates=1)

endtime = datetime.now()
print(code)
print('Total runtime was %i' % (endtime - starttime).seconds)
