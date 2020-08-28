from peakondrake import *
from netCDF4 import Dataset
from datetime import datetime

base_code = 'final_peakon_convergence_old_deterministic_strong_dt'
Ld = 40.
tmax = 10
dt_true = 1e-4
dts = [0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005]

starttime = datetime.now()

for i, dt in enumerate(dts):
    code = base_code+'_'+str(i)
    experiment(code, Ld, tmax,
               resolutions=10000,
               dts=dt,
               sigmas=0,
               seeds=0,
               schemes='hydrodynamic',
               timesteppings='midpoint',
               ics='peakon',
               num_Xis=0,
               Xi_family='double_sines',
               alphasq=1.0,
               c0=0.,
               gamma=0.,
               diagnostics=['u_error_with_sde', 'u_error_weak', 'p_pde', 'q_pde'],
               fields_to_output=['u_sde', 'u_sde_weak'],
               ndump=int(tmax / (100 * dt)),
               field_ndump=int(tmax / (1 * dt)),
               allow_fail=True,
               peakon_equations=False,
               true_peakon_data='true_old_deterministic_peakon_data',
               nXi_updates=int(dt/dt_true))

endtime = datetime.now()
print(code)
print('Total runtime was %i' % (endtime - starttime).seconds)