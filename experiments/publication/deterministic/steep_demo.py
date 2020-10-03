from peakondrake import *
from netCDF4 import Dataset
from datetime import datetime

code = 'deterministic_steep_demo'
Ld = 40.
tmax = 20
dt = 0.001

starttime = datetime.now()

experiment(code, Ld, tmax,
           resolutions=20000,
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
           diagnostics=['u_field','max_du_loc','min_du_loc'],
           fields_to_output=[],
           ndump=int(tmax / (2000 * dt)),
           field_ndump=int(tmax / (1 * dt)),
           allow_fail=True,
           nXi_updates=1,
           peak_width=(1/6))

endtime = datetime.now()
print(code)
print('Total runtime was %i' % (endtime - starttime).seconds)
