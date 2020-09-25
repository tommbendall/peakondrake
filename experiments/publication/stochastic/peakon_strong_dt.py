from peakondrake import *
from netCDF4 import Dataset
from datetime import datetime

base_code = 'convergence_strong_dt_'
Ld = 40.
tmax = 0.1
dt_true = 1e-6
dts = [0.1, 0.05, 0.025, 0.02, 0.0125, 0.01, 0.00625, 0.005, 0.004, 0.0025, 0.002, 0.00125, 0.001, 0.0008, 0.0005, 0.0004, 0.00025, 0.0002, 0.0001, 0.00005]


starttime = datetime.now()

for i, dt in enumerate(dts):
    code = base_code+str(i)
    experiment(code, Ld, tmax,
               resolutions=20000,
               dts=dt,
               sigmas=1.0,
               seeds=0,
               schemes='hydrodynamic',
               timesteppings='midpoint',
               ics='periodic_peakon',
               num_Xis=1,
               Xi_family='constant',
               alphasq=1.0,
               c0=0.,
               gamma=0.,
               diagnostics=['u_error_with_sde', 'u_error_weak', 'p_pde', 'q_pde'],
               fields_to_output=['u_sde', 'u_sde_weak'],
               ndump=int(tmax / (1 * dt)),
               field_ndump=int(tmax / (1 * dt)),
               allow_fail=True,
               peakon_equations=False,
               true_peakon_data='true_peakon_data',
               nXi_updates=int(dt/dt_true),
               periodic=True)


endtime = datetime.now()
print(code)
print('Total runtime was %i' % (endtime - starttime).seconds)
