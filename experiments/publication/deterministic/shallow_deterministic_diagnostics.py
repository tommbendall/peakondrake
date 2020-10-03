from peakondrake import *
from netCDF4 import Dataset
from datetime import datetime

code = 'shallow_deterministic_diagnostics'
Ld = 40.
tmax = 20
dt = 0.0005

starttime = datetime.now()

experiment(code, Ld, tmax,
           resolutions=[1000, 1500, 2000, 2500, 3000],
           dts=dt,
           sigmas=0.0,
           seeds=0,
           schemes='hydrodynamic',
           timesteppings='midpoint',
           ics='new_peak',
           num_Xis=1,
           Xi_family='sines',
           alphasq=1.0,
           c0=0.,
           gamma=0.,
           diagnostics=['max_du', 'min_du', 'peakon_suite'],
           fields_to_output=['du'],
           ndump=int(tmax / (1000 * dt)),
           field_ndump=int(tmax / (1 * dt)),
           allow_fail=True,
           nXi_updates=1,
           peak_width=1.0)

endtime = datetime.now()
print(code)
print('Total runtime was %i' % (endtime - starttime).seconds)
