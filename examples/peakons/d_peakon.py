from peakondrake import *
from netCDF4 import Dataset
from datetime import datetime

code = 'd_peakon'
Ld = 40.
tmax = 200
dt = 0.001

starttime = datetime.now()

experiment(code, Ld, tmax,
           resolutions=[500, 1000, 1100, 1200, 2000],
           dts=dt,
           sigmas=[0.0],
           seeds=0,
           schemes='hydrodynamic',
           timesteppings='midpoint',
           ics='d_peakon',
           num_Xis=1,
           Xi_family='sines',
           alphasq=1.0,
           c0=0.,
           gamma=0.,
           diagnostics=['p_pde', 'q_pde'],
           fields_to_output=[],
           ndump=int(tmax / (1000 * dt)),
           field_ndump=int(tmax / (10 * dt)),
           allow_fail=True,
           peakon_equations=True)

endtime = datetime.now()
print(code)
print('Total runtime was %i' % (endtime - starttime).seconds)
