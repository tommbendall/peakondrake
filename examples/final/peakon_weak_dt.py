from peakondrake import *
from netCDF4 import Dataset
from datetime import datetime

Ld = 40.
tmax = 10
dt_true = 1e-6
dts = [1e-4, 2e-4, 4e-4]
seeds = range(20)

base_true_code = 'true_weak'
base_pde_code = 'final_peakon_convergence_weak_dt'

starttime = datetime.now()

for seed in seeds:
    true_code = base_true_code+'_'+str(seed)

    experiment(true_code, Ld, tmax,
               resolutions=10000,
               dts=dt,
               sigmas=0.02,
               seeds=0,
               schemes='hydrodynamic',
               timesteppings='midpoint',
               ics='peakon',
               num_Xis=1,
               Xi_family='double_sines',
               alphasq=1.0,
               c0=0.,
               gamma=0.,
               diagnostics=['l2_u'],
               fields_to_output=[],
               ndump=int(tmax / (100 * dt)),
               field_ndump=int(tmax / (1 * dt)),
               allow_fail=True,
               peakon_equations=True,
               only_peakons=True)

    for i, dt in enumerate(dts):
        code = base_code+'_'+str(seed)+'_'+str(i)
        experiment(code, Ld, tmax,
                   resolutions=1000,
                   dts=dt,
                   sigmas=0.02,
                   seeds=seed,
                   schemes='hydrodynamic',
                   timesteppings='midpoint',
                   ics='peakon',
                   num_Xis=1,
                   Xi_family='sines',
                   alphasq=1.0,
                   c0=0.,
                   gamma=0.,
                   diagnostics=['p_pde', 'q_pde', 'u_error_with_sde', 'u_error_weak'],
                   fields_to_output=['u_sde', 'u_sde_weak'],
                   ndump=int(tmax / (100 * dt)),
                   field_ndump=int(tmax / (4 * dt)),
                   allow_fail=True,
                   peakon_equations=False,
                   true_peakon_data=true_code,
                   nXi_updates=int(dt/dt_true))

endtime = datetime.now()
print(code)
print('Total runtime was %i' % (endtime - starttime).seconds)
