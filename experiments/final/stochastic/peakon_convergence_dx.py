from peakondrake import *
from netCDF4 import Dataset
from datetime import datetime

code = 'convergence_dx'
Ld = 40.
tmax = 0.1
dt = 1e-5

starttime = datetime.now()

experiment(code, Ld, tmax,
           resolutions=[500, 750, 1000, 1500, 2000],
           dts=dt,
           sigmas=1.0,
           seeds=0,
           schemes='hydrodynamic',
           timesteppings='midpoint',
           ics='periodic_peakon',
           num_Xis=1,
           Xi_family='double_sines',
           alphasq=1.0,
           c0=0.,
           gamma=0.,
           diagnostics=['u_error_with_sde', 'u_error_weak', 'p_pde', 'q_pde'],
           fields_to_output=['u_sde', 'u_sde_weak'],
           ndump=int(tmax / (1 * dt)),
           field_ndump=int(tmax / (4 * dt)),
           allow_fail=True,
           peakon_equations=False,
           true_peakon_data='true_peakon_data',
           nXi_updates=1,
           periodic=True)

endtime = datetime.now()
print(code)
print('Total runtime was %i' % (endtime - starttime).seconds)
