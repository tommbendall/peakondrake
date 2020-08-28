from peakondrake import *
from netCDF4 import Dataset
from datetime import datetime

code = 'mu_sigma_example'
Ld = 40.
tmax = 0.1
dt = 0.001

starttime = datetime.now()

experiment(code, Ld, tmax,
           resolutions=500,
           dts=dt,
           sigmas=[0.04, 0.08, 0.12, 0.16, 0.2],
           seeds=0,
           schemes='hydrodynamic',
           timesteppings='midpoint',
           ics='one_peak',
           num_Xis=6,
           Xi_family='sines',
           alphasq=1.0,
           c0=0.,
           gamma=0.,
           diagnostics=['mu', 'l2_u'],
           fields_to_output=['du'],
           ndump=int(tmax / (2000 * dt)),
           field_ndump=int(tmax / (1 * dt)),
           allow_fail=True,
           nXi_updates=1)

endtime = datetime.now()
print(code)
print('Total runtime was %i' % (endtime - starttime).seconds)
